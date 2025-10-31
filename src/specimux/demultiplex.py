#!/usr/bin/env python3

"""
Core demultiplexing pipeline logic.

This module contains demultiplex functionality
extracted from the original specimux.py monolithic file.
"""

import argparse
import logging
from enum import Enum
from typing import List, Optional, Tuple

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


def create_write_operation(sample_id, args, seq, match, resolution_type, trace_sequence_id=None):
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
                
                # Process all equivalent matches and create write operations directly
                equivalent_count = len(best_matches)
                has_full_match = False
                for match in best_matches:
                    # Determine final sample ID for this match (includes specimen resolution and fallback logic)
                    final_sample_id, resolution_type = resolve_specimen(match, specimens, trace_logger, sequence_id,
                                                                        args, equivalent_count)
                    
                    # Create and add write operation directly - use the sequence in the matched orientation
                    # This ensures proper orientation normalization
                    op = create_write_operation(final_sample_id, args, match.sequence, match, resolution_type, sequence_id)
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
                op = create_write_operation(SampleId.UNKNOWN, args, match.sequence, match, ResolutionType.UNKNOWN, sequence_id)
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



def resolve_specimen(match: CandidateMatch, specimens: Specimens,
                     trace_logger: Optional[TraceLogger] = None,
                     sequence_id: Optional[str] = None,
                     args: Optional[argparse.Namespace] = None,
                     equivalent_matches_count: int = 1) -> Tuple[str, ResolutionType]:
    """Determine final sample ID for this match, including specimen resolution and fallback logic.
    
    Returns:
        Tuple of (sample_id, resolution_type)
    """
    sample_id = SampleId.UNKNOWN
    # Start with pool already determined from primers in match_sequence
    pool = match.get_pool()

    # Check if we should downgrade full matches to partial
    is_full_match = match.p1_match and match.p2_match and match.has_b1_match() and match.has_b2_match()
    should_downgrade = (is_full_match and equivalent_matches_count > 1 and 
                       args and args.resolve_multiple_matches == MultipleMatchStrategy.DOWNGRADE_FULL)
    
    if should_downgrade:
        # Downgrade full match to partial output
        p1_name = match.get_p1().name if match.get_p1() else 'unknown'
        p2_name = match.get_p2().name if match.get_p2() else 'unknown'
        b1_name = match.best_b1()[0] if match.best_b1() else 'unknown'
        b2_name = match.best_b2()[0] if match.best_b2() else 'unknown'
        sample_id = f"{SampleId.PREFIX_FWD_MATCH}downgraded_{p1_name}_{p2_name}_{b1_name}_{b2_name}"
        resolution_type = ResolutionType.DOWNGRADED_MULTIPLE
        
        # Log downgrade decision
        if trace_logger:
            trace_logger.log_match_discarded(sequence_id, match, 5.0, 'downgraded_multiple_full')
    elif is_full_match:
        ids = specimens.specimens_for_barcodes_and_primers(
            match.best_b1(), match.best_b2(), match.get_p1(), match.get_p2())
        if len(ids) > 1:
            # Multiple specimen matches - use first match (this shouldn't happen in new flow)
            sample_id = ids[0]  # Use first specimen ID
            resolution_type = ResolutionType.MULTIPLE_SPECIMENS
            # Use pool from first specimen
            pool = specimens.get_specimen_pool(sample_id)
        elif len(ids) == 1:
            sample_id = ids[0]
            resolution_type = ResolutionType.FULL_MATCH
            # For unique matches, use specimen's pool
            pool = specimens.get_specimen_pool(ids[0])
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
    """
    if p1 and p2:
        # Find common pools between primers
        common_pools = set(p1.pools).intersection(p2.pools)
        if common_pools:
            return list(common_pools)[0]
    elif p1:
        return p1.pools[0]
    elif p2:
        return p2.pools[0]
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

