#!/usr/bin/env python3

"""
Core demultiplexing pipeline logic.

This module contains demultiplex functionality
extracted from the original specimux.py monolithic file.
"""

import argparse
import logging
from enum import Enum
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from Bio.Seq import reverse_complement
from Bio.SeqRecord import SeqRecord

from .constants import (
    AlignMode, Barcode, Primer, ResolutionType, SampleId, Orientation, TrimMode, MultipleMatchStrategy
)
from .databases import BarcodePrefilter, Specimens
from .models import (
    AlignmentResult, CandidateMatch, MatchParameters, PrimerInfo, WriteOperation
)
from .alignment import align_seq, get_quality_seq
from .trace import TraceLogger


def create_write_operation(sample_id, args, seq, match, resolution_type, trace_sequence_id=None,
                           trace_logger: Optional[TraceLogger] = None):
    formatted_seq = seq.seq
    quality_scores = get_quality_seq(seq)

    s = 0
    e = len(formatted_seq)

    if args.trim == TrimMode.PRIMERS:
        (s, e) = match.interprimer_extent()
    elif args.trim == TrimMode.BARCODES:
        (s, e) = match.interbarcode_extent()
    elif args.trim == TrimMode.TAILS:
        (s, e) = match.intertail_extent()

    if args.trim != TrimMode.NONE:
        # Check if trimming would result in empty sequence
        if s >= e:
            logging.debug(f"Sequence {seq.id}: trimming would produce empty sequence (trim region {s}:{e}), routing to fallback")
            # Log trace event
            if trace_logger:
                p1_name = match.get_p1().name if match.get_p1() else "unknown"
                p2_name = match.get_p2().name if match.get_p2() else "unknown"
                trace_logger.log_sequence_trim_empty(
                    trace_sequence_id, args.trim, s, e, len(formatted_seq), p1_name, p2_name)
            # Return fallback WriteOperation with untrimmed sequence to unknown/unknown/unknown-unknown/
            quality_seq = "".join(chr(q + 33) for q in quality_scores)
            return WriteOperation(
                sample_id=SampleId.UNKNOWN,
                seq_id=seq.id,
                distance_code=match.distance_code(),
                sequence=str(formatted_seq),
                quality_sequence=quality_seq,
                quality_scores=quality_scores,
                p1_location=match.get_p1_location(),
                p2_location=match.get_p2_location(),
                b1_location=match.get_barcode1_location(),
                b2_location=match.get_barcode2_location(),
                primer_pool="unknown",
                p1_name="unknown",
                p2_name="unknown",
                resolution_type=ResolutionType.UNKNOWN,
                trace_sequence_id=trace_sequence_id
            )
        formatted_seq = seq.seq[s:e]
        quality_scores = quality_scores[s:e]
        match.trim_locations(s)

    quality_seq = "".join(chr(q + 33) for q in quality_scores)

    # Get primer names, using "unknown" if not matched
    p1_name = match.get_p1().name if match.get_p1() else "unknown"
    p2_name = match.get_p2().name if match.get_p2() else "unknown"

    # Get pool name - could come from specimen or default to "unknown"
    primer_pool = match.get_pool() if match.get_pool() else "unknown"

    return WriteOperation(
        sample_id=sample_id,
        seq_id=seq.id,
        distance_code=match.distance_code(),
        sequence=str(formatted_seq),
        quality_sequence=quality_seq,
        quality_scores=quality_scores,
        p1_location=match.get_p1_location(),
        p2_location=match.get_p2_location(),
        b1_location=match.get_barcode1_location(),
        b2_location=match.get_barcode2_location(),
        primer_pool=primer_pool,
        p1_name=p1_name,
        p2_name=p2_name,
        resolution_type=resolution_type,
        trace_sequence_id=trace_sequence_id
    )




def process_sequences(seq_records: List[SeqRecord],
                      parameters: MatchParameters,
                      specimens: Specimens,
                      args: argparse.Namespace,
                      prefilter: Optional[BarcodePrefilter],
                      trace_logger: Optional[TraceLogger] = None,
                      record_offset: int = 0) -> Tuple[List[WriteOperation], int, int]:
    """Process sequences and track pipeline events.
    
    Returns:
        write_ops: List of write operations to perform
        total_count: Number of sequences processed
        matched_count: Number of sequences successfully matched
    """
    write_ops = []
    total_count = 0
    matched_count = 0

    for idx, seq in enumerate(seq_records):
        total_count += 1
        
        # Generate unique sequence ID for tracing
        sequence_id = None
        if trace_logger:
            sequence_id = trace_logger.get_sequence_id(seq, record_offset + idx)
            trace_logger.log_sequence_received(sequence_id, len(seq), seq.id)
        
        if args.min_length != -1 and len(seq) < args.min_length:
            if trace_logger:
                trace_logger.log_sequence_filtered(sequence_id, len(seq), 'too_short')
        elif args.max_length != -1 and len(seq) > args.max_length:
            if trace_logger:
                trace_logger.log_sequence_filtered(sequence_id, len(seq), 'too_long')
        else:
            rseq = seq.reverse_complement()
            rseq.id = seq.id
            rseq.description = seq.description

            matches = find_candidate_matches(prefilter, parameters, seq, rseq, specimens, trace_logger, sequence_id)

            # Handle match selection and processing
            if matches:
                best_matches = select_best_matches(matches, trace_logger, sequence_id)

                # Process matches based on strategy
                has_full_match = False

                if args.dereplicate == MultipleMatchStrategy.BEST:
                    # Dereplicate: group by specimen, select best match per specimen
                    derep_results = dereplicate_matches(best_matches, specimens, trace_logger, sequence_id)

                    for match, specimen_id, b1, b2 in derep_results:
                        if specimen_id is not None:
                            # Full match with specific specimen - use directly
                            resolution_type = ResolutionType.DEREPLICATED_FULL
                            pool = specimens.get_specimen_pool(specimen_id)
                            match.set_pool(pool)

                            op = create_write_operation(specimen_id, args, match.sequence, match,
                                                        resolution_type, sequence_id, trace_logger)
                            if op is not None:
                                write_ops.append(op)
                            has_full_match = True
                        else:
                            # Partial/unknown - handle normally via resolve_specimen
                            final_sample_id, resolution_type = resolve_specimen(
                                match, specimens, trace_logger, sequence_id)
                            op = create_write_operation(final_sample_id, args, match.sequence, match,
                                                        resolution_type, sequence_id, trace_logger)
                            if op is not None:
                                write_ops.append(op)
                            if resolution_type.is_full_match():
                                has_full_match = True
                else:
                    # No dereplication: process all equivalent matches
                    for match in best_matches:
                        # Determine final sample ID for this match (includes specimen resolution and fallback logic)
                        final_sample_id, resolution_type = resolve_specimen(
                            match, specimens, trace_logger, sequence_id)

                        # Create and add write operation directly - use the sequence in the matched orientation
                        # This ensures proper orientation normalization
                        op = create_write_operation(final_sample_id, args, match.sequence, match, resolution_type,
                                                     sequence_id, trace_logger)
                        if op is not None:
                            write_ops.append(op)

                        # Count successful matches
                        if resolution_type.is_full_match():
                            has_full_match = True

                # Only count sequences with at least one full match as successful
                if has_full_match:
                    matched_count += 1
            else:
                # No matches found - create minimal match object for output
                match = CandidateMatch(seq, specimens.b_length())
                if trace_logger:
                    trace_logger.log_no_match_found(sequence_id, 'primer_search', 'No primer matches found')
                op = create_write_operation(SampleId.UNKNOWN, args, match.sequence, match, ResolutionType.UNKNOWN,
                                            sequence_id, trace_logger)
                if op is not None:
                    write_ops.append(op)

    return write_ops, total_count, matched_count



def select_best_matches(matches: List[CandidateMatch],
                        trace_logger: Optional[TraceLogger] = None,
                        sequence_id: Optional[str] = None) -> List[CandidateMatch]:
    """Select all equivalent best matches based on match quality scoring"""
    if not matches:
        raise ValueError("No matches provided to choose_best_match")

    # Score all matches and group by score
    scored_matches = []
    for m in matches:
        score = 0
        if m.p1_match and m.p2_match and m.has_b1_match() and m.has_b2_match():
            score = 5
        elif m.p1_match and m.p2_match and (m.has_b1_match() or m.has_b2_match()):
            score = 4
        elif (m.p1_match or m.p2_match) and (m.has_b1_match() or m.has_b2_match()):
            score = 3
        elif m.p1_match and m.p2_match:
            score = 2
        elif m.p1_match or m.p2_match:
            score = 1
        
        # Log match scoring
        if trace_logger:
            trace_logger.log_match_scored(sequence_id, m, float(score))
        
        scored_matches.append((score, m))
    
    # Sort by score (descending)
    scored_matches.sort(key=lambda x: x[0], reverse=True)
    best_score = scored_matches[0][0]
    
    # Collect all matches with best score - these are all equivalent
    equivalent_matches = [m for score, m in scored_matches if score == best_score]
    
    # Log discarded matches (those with lower scores)
    if trace_logger:
        for score, m in scored_matches:
            if score < best_score:
                trace_logger.log_match_discarded(sequence_id, m, float(score), 'lower_score')
    
    # Note: Multiple matches are handled downstream in match_sample and logged via MATCH_SELECTED/MATCH_DISCARDED events

    return equivalent_matches


def dereplicate_matches(matches: List[CandidateMatch],
                        specimens: Specimens,
                        trace_logger: Optional[TraceLogger] = None,
                        sequence_id: Optional[str] = None) -> List[Tuple[CandidateMatch, str, str, str]]:
    """Dereplicate matches by grouping equivalent matches by specimen ID.

    For each unique specimen that could be matched, select the best match using:
    1. Total barcode edit distance (b1_distance + b2_distance) - lower wins
    2. Total primer edit distance (p1_distance + p2_distance) - lower wins
    3. File order of primers (p1.file_index + p2.file_index) - lower wins

    Args:
        matches: List of equivalent CandidateMatches (all with same score)
        specimens: Specimens database for lookups
        trace_logger: Optional trace logger
        sequence_id: Sequence ID for tracing

    Returns:
        List of (CandidateMatch, specimen_id, b1, b2) tuples - one per unique specimen.
        For partial matches (no specimen found), returns (match, None, None, None).
    """
    # Step 1: Expand matches to all possible (match, specimen_id, b1, b2, b1_dist, b2_dist) combinations
    expanded: List[Tuple[CandidateMatch, Optional[str], Optional[str], Optional[str], float, float]] = []

    for match in matches:
        if not match.has_full_match():
            # For partial matches, keep as-is (handled downstream)
            expanded.append((match, None, None, None, 999.0, 999.0))
            continue

        p1 = match.get_p1()
        p2 = match.get_p2()

        # Get all equally-good barcodes
        b1_candidates = match.best_b1()  # List of barcodes within tolerance
        b2_candidates = match.best_b2()

        # Get the distance dictionaries for specific barcode lookups
        b1_distances = match.get_barcode_distances(Barcode.B1)
        b2_distances = match.get_barcode_distances(Barcode.B2)

        # For each barcode combination, find the specimen
        found_any = False
        for b1 in b1_candidates:
            for b2 in b2_candidates:
                specimen_id = specimens.specimen_for_exact_match(b1, b2, p1, p2)
                if specimen_id:
                    b1_dist = b1_distances.get(b1, 999.0)
                    b2_dist = b2_distances.get(b2, 999.0)
                    expanded.append((match, specimen_id, b1, b2, b1_dist, b2_dist))
                    found_any = True

        # If no specimen found for any barcode combo, treat as partial
        if not found_any:
            expanded.append((match, None, None, None, 999.0, 999.0))

    # Log expansion
    if trace_logger:
        trace_logger.log_dereplicate_expanded(sequence_id, len(matches), len(expanded))

    # Step 2: Group by specimen_id
    specimen_groups: Dict[Optional[str], List[Tuple[CandidateMatch, Optional[str], Optional[str], Optional[str], float, float]]] = defaultdict(list)
    for entry in expanded:
        specimen_id = entry[1]
        specimen_groups[specimen_id].append(entry)

    # Step 3: For each specimen group, select the best match
    results: List[Tuple[CandidateMatch, str, str, str]] = []

    for specimen_id, group in specimen_groups.items():
        if specimen_id is None:
            # Separate matches by barcode status
            def has_single_barcode(m: CandidateMatch) -> bool:
                has_b1 = m.has_b1_match()
                has_b2 = m.has_b2_match()
                return (has_b1 and not has_b2) or (has_b2 and not has_b1)

            def has_no_barcodes(m: CandidateMatch) -> bool:
                return not m.has_b1_match() and not m.has_b2_match()

            single_barcode_matches = [entry[0] for entry in group if has_single_barcode(entry[0])]
            no_barcode_matches = [entry[0] for entry in group if has_no_barcodes(entry[0])]
            other_matches = [entry[0] for entry in group
                             if not has_single_barcode(entry[0]) and not has_no_barcodes(entry[0])]

            # Dereplicate single-barcode matches (partials)
            if single_barcode_matches:
                derep_partials = dereplicate_partial_matches(
                    single_barcode_matches, trace_logger, sequence_id)
                for match in derep_partials:
                    results.append((match, None, None, None))

            # Dereplicate no-barcode matches (unknowns)
            if no_barcode_matches:
                derep_unknowns = dereplicate_unknown_matches(
                    no_barcode_matches, trace_logger, sequence_id)
                for match in derep_unknowns:
                    results.append((match, None, None, None))

            # Keep other matches as-is (full match without specimen - rare edge case)
            for match in other_matches:
                results.append((match, None, None, None))

            continue

        # Sort by tiebreaker hierarchy:
        # 1. Total barcode distance (b1_dist + b2_dist) - lower wins
        # 2. Total primer distance - lower wins
        # 3. File index sum - lower wins
        def sort_key(entry: Tuple[CandidateMatch, Optional[str], Optional[str], Optional[str], float, float]) -> Tuple[float, int, int]:
            match, _, b1, b2, b1_dist, b2_dist = entry
            barcode_dist = b1_dist + b2_dist
            primer_dist = match.primer_distance(Primer.FWD) + match.primer_distance(Primer.REV)
            p1 = match.get_p1()
            p2 = match.get_p2()
            file_index = (p1.file_index if p1 else 999) + (p2.file_index if p2 else 999)
            return (barcode_dist, primer_dist, file_index)

        group.sort(key=sort_key)
        best = group[0]
        results.append((best[0], specimen_id, best[2], best[3]))

        # Log selection
        if trace_logger:
            scores = sort_key(best)
            trace_logger.log_dereplicate_selected(
                sequence_id, specimen_id,
                len(group),  # alternatives considered
                scores  # winning scores (barcode_dist, primer_dist, file_idx)
            )

    return results


def dereplicate_partial_matches(
        matches: List[CandidateMatch],
        trace_logger: Optional[TraceLogger] = None,
        sequence_id: Optional[str] = None) -> List[CandidateMatch]:
    """Dereplicate partial matches by grouping by matched barcode.

    For matches with exactly one barcode (regardless of primer count), group by
    the matched barcode and select the best match per barcode using:
    1. Barcode edit distance (only the matched one)
    2. Combined primer distance (p1_distance + p2_distance)
    3. File order sum (p1.file_index + p2.file_index)

    Args:
        matches: List of partial CandidateMatches
        trace_logger: Optional trace logger
        sequence_id: Sequence ID for tracing

    Returns:
        List of dereplicated CandidateMatches - one per unique barcode
    """
    # Group by (direction, barcode)
    barcode_groups: Dict[Tuple[str, str], List[CandidateMatch]] = defaultdict(list)

    for match in matches:
        # Handle any match with exactly one barcode (regardless of primer count)
        has_b1 = match.has_b1_match()
        has_b2 = match.has_b2_match()

        if has_b1 and not has_b2:
            direction = 'forward'
            barcodes = match.best_b1()
        elif has_b2 and not has_b1:
            direction = 'reverse'
            barcodes = match.best_b2()
        else:
            # Either both barcodes (full match) or neither (unknown) - skip
            continue

        # Add to group for each equally-good barcode
        for barcode in barcodes:
            barcode_groups[(direction, barcode)].append(match)

    # Select best match per barcode group
    results: List[CandidateMatch] = []

    for (direction, barcode), group in barcode_groups.items():
        def sort_key(m: CandidateMatch) -> Tuple[float, int, int, int]:
            # Barcode distance (only the matched one)
            bc_dist = m.b1_distance() if direction == 'forward' else m.b2_distance()

            # Count primer matches (negate so higher count sorts first)
            primer_count = 0
            if m.get_p1():
                primer_count += 1
            if m.get_p2():
                primer_count += 1

            # Combined primer distance
            primer_dist = 0
            file_idx = 0
            if m.get_p1():
                primer_dist += m.primer_distance(Primer.FWD)
                file_idx += m.get_p1().file_index
            if m.get_p2():
                primer_dist += m.primer_distance(Primer.REV)
                file_idx += m.get_p2().file_index

            return (bc_dist, -primer_count, primer_dist, file_idx)

        group.sort(key=sort_key)
        best = group[0]
        results.append(best)

        # Log selection
        if trace_logger:
            trace_logger.log_dereplicate_partial_selected(
                sequence_id, direction, barcode,
                len(group),  # alternatives considered
                sort_key(best)  # winning scores
            )

    return results


def dereplicate_unknown_matches(
        matches: List[CandidateMatch],
        trace_logger: Optional[TraceLogger] = None,
        sequence_id: Optional[str] = None) -> List[CandidateMatch]:
    """Dereplicate unknown matches (no barcodes) to a single best match.

    For matches with no barcodes, select the best match using:
    1. Number of primer matches (2 > 1 > 0) - higher wins
    2. Combined primer distance (p1_distance + p2_distance) - lower wins
    3. File order sum (p1.file_index + p2.file_index) - lower wins

    Args:
        matches: List of unknown CandidateMatches (no barcodes)
        trace_logger: Optional trace logger
        sequence_id: Sequence ID for tracing

    Returns:
        List containing at most one CandidateMatch (the best one)
    """
    if not matches:
        return []

    def sort_key(m: CandidateMatch) -> Tuple[int, int, int]:
        # Count primer matches (negate so higher count sorts first)
        primer_count = 0
        if m.get_p1():
            primer_count += 1
        if m.get_p2():
            primer_count += 1

        # Combined primer distance
        primer_dist = 0
        if m.get_p1():
            primer_dist += m.primer_distance(Primer.FWD)
        if m.get_p2():
            primer_dist += m.primer_distance(Primer.REV)

        # File order sum
        file_idx = 0
        if m.get_p1():
            file_idx += m.get_p1().file_index
        else:
            file_idx += 999
        if m.get_p2():
            file_idx += m.get_p2().file_index
        else:
            file_idx += 999

        return (-primer_count, primer_dist, file_idx)

    matches_sorted = sorted(matches, key=sort_key)
    best = matches_sorted[0]

    if trace_logger and len(matches) > 1:
        scores = sort_key(best)
        trace_logger.log_dereplicate_unknown_selected(
            sequence_id, len(matches), -scores[0], scores[1], scores[2])

    return [best]


def resolve_specimen(match: CandidateMatch, specimens: Specimens,
                     trace_logger: Optional[TraceLogger] = None,
                     sequence_id: Optional[str] = None) -> Tuple[str, ResolutionType]:
    """Determine final sample ID for this match, including specimen resolution and fallback logic.
    
    Returns:
        Tuple of (sample_id, resolution_type)
    """
    sample_id = SampleId.UNKNOWN
    # Start with pool already determined from primers in match_sequence
    pool = match.get_pool()

    is_full_match = match.p1_match and match.p2_match and match.has_b1_match() and match.has_b2_match()

    if is_full_match:
        ids = specimens.specimens_for_barcodes_and_primers(
            match.best_b1(), match.best_b2(), match.get_p1(), match.get_p2())
        if len(ids) > 1:
            # Multiple specimen matches - use first match (this shouldn't happen in new flow)
            sample_id = ids[0]  # Use first specimen ID
            resolution_type = ResolutionType.MULTIPLE_SPECIMENS
            # Use pool from first specimen
            pool = specimens.get_specimen_pool(sample_id)
            match.set_pool(pool)
        elif len(ids) == 1:
            sample_id = ids[0]
            resolution_type = ResolutionType.FULL_MATCH
            # For unique matches, use specimen's pool
            pool = specimens.get_specimen_pool(ids[0])
            match.set_pool(pool)
        else:
            logging.warning(
                f"No Specimens for combo: ({match.best_b1()}, {match.best_b2()}, "
                f"{match.get_p1()}, {match.get_p2()})")
            resolution_type = ResolutionType.UNKNOWN
    else:
        # Partial match - determine type and generate appropriate sample ID
        b1s = match.best_b1()
        b2s = match.best_b2()
        
        if match.has_b1_match() and not match.has_b2_match() and len(b1s) == 1:
            resolution_type = ResolutionType.PARTIAL_FORWARD
            sample_id = SampleId.PREFIX_FWD_MATCH + b1s[0]
        elif match.has_b2_match() and not match.has_b1_match() and len(b2s) == 1:
            resolution_type = ResolutionType.PARTIAL_REVERSE
            sample_id = SampleId.PREFIX_REV_MATCH + b2s[0]
        else:
            resolution_type = ResolutionType.UNKNOWN
            # Keep sample_id as SampleId.UNKNOWN
    
    # Log specimen resolution
    if trace_logger:
        trace_logger.log_specimen_resolved(sequence_id, match, sample_id, 
                                         resolution_type.to_string(), pool or 'none')


    
    return sample_id, resolution_type



def determine_orientation(parameters: MatchParameters, seq: str, rseq: str,
                          fwd_primers: List[PrimerInfo], 
                          rev_primers: List[PrimerInfo]) -> Tuple[Orientation, int, int]:
    """Determine sequence orientation by checking primer matches, returning scores.
    Returns (orientation, forward_score, reverse_score)."""

    forward_matches = 0
    reverse_matches = 0

    for primer in fwd_primers:
        fwd_p1 = align_seq(primer.primer, seq, parameters.max_dist_primers[primer.primer],
                           0, parameters.search_len)
        rev_p1 = align_seq(primer.primer, rseq, parameters.max_dist_primers[primer.primer],
                           0, parameters.search_len)
        if fwd_p1.matched():
            forward_matches += 1
        if rev_p1.matched():
            reverse_matches += 1
    for primer in rev_primers:
        fwd_p2 = align_seq(primer.primer, rseq, parameters.max_dist_primers[primer.primer],
                           0, parameters.search_len)
        rev_p2 = align_seq(primer.primer, seq, parameters.max_dist_primers[primer.primer],
                           0, parameters.search_len)
        if fwd_p2.matched():
            forward_matches += 1
        if rev_p2.matched():
            reverse_matches += 1

    # Determine orientation
    if forward_matches > 0 and reverse_matches == 0:
        orientation = Orientation.FORWARD
    elif reverse_matches > 0 and forward_matches == 0:
        orientation = Orientation.REVERSE
    else:
        orientation = Orientation.UNKNOWN
    
    return orientation, forward_matches, reverse_matches

def get_pool_from_primers(p1: Optional[PrimerInfo], p2: Optional[PrimerInfo]) -> Optional[str]:
    """
    Determine the primer pool based on primers used for match.

    Args:
        p1: Forward primer info (or None)
        p2: Reverse primer info (or None)

    Returns:
        str: Pool name if one can be determined, None otherwise

    Note:
        Uses alphabetical sorting for deterministic selection when multiple pools are possible.
        For full matches, this initial pool assignment may be overridden by the specimen's
        declared pool in resolve_specimen().
    """
    if p1 and p2:
        # Find common pools between primers
        common_pools = set(p1.pools).intersection(p2.pools)
        if common_pools:
            return sorted(common_pools)[0]  # Deterministic alphabetical selection
    elif p1:
        return sorted(p1.pools)[0] if p1.pools else None
    elif p2:
        return sorted(p2.pools)[0] if p2.pools else None
    return None


def find_candidate_matches(prefilter: Optional[BarcodePrefilter], parameters: MatchParameters, seq: SeqRecord,
                           rseq: SeqRecord, specimens: Specimens,
                           trace_logger: Optional[TraceLogger] = None,
                           sequence_id: Optional[str] = None) -> List[CandidateMatch]:
    """Match sequence against primers and barcodes"""
    # extract string versions for performance - roughly 18% improvement
    s = str(seq.seq)
    rs = str(rseq.seq)

    if parameters.preorient:
        orientation, fwd_score, rev_score = determine_orientation(
            parameters, s, rs,
            specimens.get_primers(Primer.FWD),
            specimens.get_primers(Primer.REV))
        
        # Log orientation detection for trace
        if trace_logger:
            confidence = 0.0
            if fwd_score + rev_score > 0:
                confidence = abs(fwd_score - rev_score) / (fwd_score + rev_score)
            trace_logger.log_orientation_detected(sequence_id, orientation.to_string(), 
                                                 fwd_score, rev_score, confidence)
    else:
        orientation = Orientation.UNKNOWN
        # When not pre-orienting, orientation is unknown
        if trace_logger:
            trace_logger.log_orientation_detected(sequence_id, orientation.to_string(), 0, 0, 0.0)

    matches = []
    match_counter = 0

    for fwd_primer in specimens.get_primers(Primer.FWD):
        for rev_primer in specimens.get_paired_primers(fwd_primer.primer):
            if orientation in [Orientation.FORWARD, Orientation.UNKNOWN]:
                candidate_match_id = f"{sequence_id}_match_{match_counter}"
                match = CandidateMatch(seq, specimens.b_length(), candidate_match_id)
                match_one_end(prefilter, match, parameters, rs, True, fwd_primer,
                              Primer.FWD, Barcode.B1, trace_logger, sequence_id)
                match_one_end(prefilter, match, parameters, s, False, rev_primer,
                              Primer.REV, Barcode.B2, trace_logger, sequence_id)
                # Only add matches where at least one primer was found
                if match.p1_match or match.p2_match:
                    # Determine and set pool for this match
                    pool = get_pool_from_primers(fwd_primer, rev_primer)
                    match.set_pool(pool)
                    
                    # Log primer match result
                    if trace_logger:
                        trace_logger.log_primer_matched(sequence_id, match, pool or 'none', 'as_is')
                        trace_logger.log_barcode_matched(sequence_id, match)
                    
                    matches.append(match)
                    match_counter += 1

            if orientation in [Orientation.REVERSE, Orientation.UNKNOWN]:
                candidate_match_id = f"{sequence_id}_match_{match_counter}"
                match = CandidateMatch(rseq, specimens.b_length(), candidate_match_id)
                match_one_end(prefilter, match, parameters, s, True, fwd_primer,
                              Primer.FWD, Barcode.B1, trace_logger, sequence_id)
                match_one_end(prefilter, match, parameters, rs, False, rev_primer,
                              Primer.REV, Barcode.B2, trace_logger, sequence_id)
                # Only add matches where at least one primer was found
                if match.p1_match or match.p2_match:
                    # Determine and set pool for this match
                    pool = get_pool_from_primers(fwd_primer, rev_primer)
                    match.set_pool(pool)
                    
                    # Log primer match result
                    if trace_logger:
                        trace_logger.log_primer_matched(sequence_id, match, pool or 'none', 'reverse_complement')
                        trace_logger.log_barcode_matched(sequence_id, match)
                    
                    matches.append(match)
                    match_counter += 1
    
    # TODO: Consider whether primer pairs that match almost exactly the same extent 
    # should be considered distinct matches or not. This affects multiple match detection
    # for cases where multiple primer pairs cover nearly identical sequence regions.
    return matches

def match_one_end(prefilter: Optional[BarcodePrefilter], match: CandidateMatch, parameters: MatchParameters, sequence: str,
                  reversed_sequence: bool, primer_info: PrimerInfo,
                  which_primer: Primer, which_barcode: Barcode,
                  trace_logger: Optional[TraceLogger] = None,
                  sequence_id: Optional[str] = None) -> None:
    """Match primers and barcodes at one end of the sequence."""

    primer = primer_info.primer
    primer_rc = primer_info.primer_rc
    search_start = len(sequence) - parameters.search_len
    search_end = len(sequence)
    
    # Log primer search attempt
    if trace_logger:
        trace_logger.log_primer_search(sequence_id, primer_info.name, which_primer.to_string(),
                                      search_start, search_end, False, -1, -1)

    primer_match = align_seq(primer_rc, sequence, parameters.max_dist_primers[primer],
                             search_start, search_end)

    # If we found matching primers, look for corresponding barcodes
    if primer_match.matched():
        match.set_primer_match(primer_match, primer_info, reversed_sequence, which_primer)
        
        # Log successful primer match
        if trace_logger:
            match_pos = primer_match.location()[0] if primer_match.location() else -1
            trace_logger.log_primer_search(sequence_id, primer_info.name, which_primer.to_string(),
                                          search_start, search_end, True, primer_match.distance(), match_pos)

        # Get relevant barcodes for this primer pair
        barcodes = primer_info.barcodes

        for b in barcodes:
            b_rc = reverse_complement(b)
            bd = None
            bm = None
            bc = None
            for l in primer_match.locations():
                barcode_search_start = l[1] + 1
                barcode_search_end = len(sequence)
                target_seq = sequence[barcode_search_start:]
                
                # Log barcode search attempt 
                if trace_logger:
                    trace_logger.log_barcode_search(sequence_id, b, which_barcode.to_string(), primer_info.name,
                                                   barcode_search_start, barcode_search_end, False, -1, -1)
                
                if prefilter and not prefilter.match(b_rc, target_seq):
                    continue

                barcode_match = align_seq(b_rc, sequence, parameters.max_dist_index,
                                          barcode_search_start, barcode_search_end, AlignMode.PREFIX)
                if barcode_match.matched():
                    # Log successful barcode search
                    if trace_logger:
                        match_pos = barcode_match.location()[0] if barcode_match.location() else -1
                        trace_logger.log_barcode_search(sequence_id, b, which_barcode.to_string(), primer_info.name,
                                                       barcode_search_start, barcode_search_end, True, 
                                                       barcode_match.distance(), match_pos)
                    
                    if bd is None or barcode_match.distance() < bd:
                        bm = barcode_match
                        bd = barcode_match.distance()
                        bc = b
                        
            if bm:
                match.add_barcode_match(bm, bc, reversed_sequence, which_barcode)
    else:
        # Log failed primer search  
        if trace_logger:
            trace_logger.log_primer_search(sequence_id, primer_info.name, which_primer.to_string(),
                                          search_start, search_end, False, -1, -1)

