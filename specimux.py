#!/usr/bin/env python

"""
Specimux: Demultiplex MinION sequences by dual barcode indexes and primers.

This program is a substantial revision and enhancement of minibar.py, originally developed by:
California Academy of Sciences, Institute for Biodiversity Science & Sustainability
(https://github.com/calacademy-research/minibar)

Original minibar.py license:
BSD 2-Clause License
Copyright (c) 2018, California Academy of Sciences, Institute for Biodiversity Science & Sustainability
All rights reserved.

Specimux extends and modifies the original minibar functionality with algorithmic enhancements
and code restructuring. While it builds upon the core concepts of minibar,
Specimux represents a significant departure from the original codebase.

For full terms of use and distribution, please see the LICENSE file accompanying this program.
"""

import argparse
import csv
import itertools
import logging
import multiprocessing
import os
import sys
import timeit
from collections import Counter
from contextlib import ExitStack
from functools import partial
from operator import itemgetter
from typing import List, Tuple
from typing import NamedTuple
from enum import Enum

from Bio import SeqIO
from Bio.Seq import reverse_complement
from multiprocessing import Pool

MODE_GLOBAL = 'NW'
MODE_INFIX = 'HW'
MODE_PREFIX = 'SHW'

SAMPLE_ID_AMBIGUOUS = "ambiguous"
SAMPLE_ID_UNKNOWN = "unknown"

try:
    import edlib
except:
    logging.error("\n    edlib module not found, to install use:\n\n        pip install edlib\n")
    exit(3)


#wrapper around edlib match result
class MatchResult:
    def __init__(self, edlib_match):
        self._edlib_match = edlib_match

    def matched(self):
        return self._edlib_match['editDistance'] > -1

    def distance(self):
        return self._edlib_match['editDistance']

    def location(self):
        return self._edlib_match['locations'][0]

    def locations(self):
        return self._edlib_match['locations']

    def reverse(self, seq_length):
        if self.distance() == -1: return

        r = self._edlib_match
        r['locations'] = [(seq_length - loc[1] - 1, seq_length - loc[0] - 1) for loc in r['locations']]

    def adjust_start(self, s):
        if self.distance() == -1: return

        self._edlib_match['locations'] = [(loc[0] + s, loc[1] + s) for loc in self._edlib_match['locations']]  # adjust relative match to absolute

    def adjust_distance(self, d):
        self._edlib_match['editDistance'] = d

class Barcode(Enum):
    B1 = 1
    B2 = 2

class Primer(Enum):
    P1 = 3
    P2 = 4

class SequenceMatch:

    def __init__(self, sequence_length, ambiguity_threshold, barcode_length):
        self.sequence_length = sequence_length
        self.p1_match: MatchResult = None
        self._b1_matches: List[Tuple[str, MatchResult, float]] = []
        self.p2_match: MatchResult = None
        self._b2_matches: List[Tuple[str, MatchResult, float]] = []
        self._p1: str = None
        self._p2: str = None
        self.ambiguity_threshold = ambiguity_threshold
        self._barcode_length = barcode_length

    def add_barcode_match(self, match: MatchResult, barcode: str, reverse: bool, which: Barcode):
        if reverse:
            match.reverse(self.sequence_length)

        if which is Barcode.B1:
            self.add_b1_match(match, barcode)
        else:
            self.add_b2_match(match, barcode)

    def add_b1_match(self, match: MatchResult, barcode: str):
        self._b1_matches.append((barcode, match, match.distance()))
        self._b1_matches.sort(key=itemgetter(2))  # Sort by edit distance

    def add_b2_match(self, match: MatchResult, barcode: str):
        self._b2_matches.append((barcode, match, match.distance()))
        self._b2_matches.sort(key=itemgetter(2))  # Sort by edit distance

    def p1_distance(self):
        return self.p1_match.distance() if self.p1_match else -1

    def p2_distance(self):
        return self.p2_match.distance() if self.p2_match else -1

    def b1_distance(self):
        return self._b1_matches[0][2] if len(self._b1_matches) > 0 else -1

    def b2_distance(self):
        return self._b2_matches[0][2] if len(self._b2_matches) > 0 else -1

    def best_b1(self) -> List[str]:
        if not self._b1_matches:
            return []
        best_distance = self.b1_distance()
        return [b for b, _, d in self._b1_matches if abs(d - best_distance) < self.ambiguity_threshold]

    def best_b1_match(self) -> List[MatchResult]:
        if not self._b1_matches:
            return []
        best_distance = self.b1_distance()
        return [m for b, m, d in self._b1_matches if abs(d - best_distance) < self.ambiguity_threshold]

    def best_b2(self) -> List[str]:
        if not self._b2_matches:
            return []
        best_distance = self.b2_distance()
        return [b for b, _, d in self._b2_matches if abs(d - best_distance) < self.ambiguity_threshold]

    def best_b2_match(self) -> List[MatchResult]:
        if not self._b2_matches:
            return []
        best_distance = self.b2_distance()
        return [m for b, m, d in self._b2_matches if abs(d - best_distance) < self.ambiguity_threshold]

    def set_primer_match(self, match: MatchResult, primer: str, reverse: bool, which: Primer):
        if reverse:
            match.reverse(self.sequence_length)
        if which is Primer.P1:
            self.set_p1_match(match, primer)
        else:
            self.set_p2_match(match, primer)

    def set_p1_match(self, match, primer):
        self.p1_match = match
        self._p1 = primer

    def set_p2_match(self, match, primer):
        self.p2_match = match
        self._p2 = primer

    def get_p1_location(self):
        return self.p1_match.location() if self.p1_match else None

    def get_barcode1_location(self):
        return self._b1_matches[0][1].location() if len(self._b1_matches) > 0 else None

    def get_p2_location(self):
        return self.p2_match.location() if self.p2_match else None

    def get_barcode2_location(self):
        return self._b2_matches[0][1].location() if len(self._b2_matches) > 0 else None

    def has_b1_match(self):
        return len(self._b1_matches) > 0

    def has_b2_match(self):
        return len(self._b2_matches) > 0

    def interprimer_extent(self):
        s = 0
        e = self.sequence_length
        if self.p1_match:
            s = self.p1_match.location()[1] + 1
        if self.p2_match:
            e = self.p2_match.location()[0]

        return (s, e)

    def interbarcode_extent(self):
        s = 0
        e = self.sequence_length
        if self.p1_match:
            s = self.p1_match.location()[0]
        if self.p2_match:
            e = self.p2_match.location()[1] + 1

        return (s, e)

    def intertail_extent(self):
        ps, pe = self.interprimer_extent()

        s = -1
        e = -1
        for b in self._b1_matches:
            for l in b[1].locations():
                if s == -1:
                    s = l[0]
                else:
                    s = min(s, l[0])
        for b in self._b2_matches:
            for l in b[1].locations():
                if e == -1:
                    e = l[1] + 1
                else:
                    e = max(e, l[1] + 1)

        if s == -1: s = max(0, ps - self._barcode_length)
        if e == -1: e = min(self.sequence_length, pe + self._barcode_length)

        return (s, e)

    def trim_locations(self, start):
        for b in self._b1_matches: b[1].adjust_start(-1 * start)
        for b in self._b2_matches: b[1].adjust_start(-1 * start)
        if self.p1_match:
            self.p1_match.adjust_start(-1 * start)
        if self.p2_match:
            self.p2_match.adjust_start(-1 * start)

    def distance_code(self):
        p1d = self.p1_match.distance() if self.p1_match else -1
        p2d = self.p2_match.distance() if self.p2_match else -1
        b1d = self.b1_distance()
        b2d = self.b2_distance()
        return format(f"({p1d},{b1d},{p2d},{b2d})")

class Specimens:
    def __init__(self, p1, p2):
        self._specimens = []
        self._b1_specimens_index = {}
        self._b2_specimens_index = {}
        self._p1 = p1
        self._p1_rc = reverse_complement(p1)
        self._p2 = p2
        self._p2_rc = reverse_complement(p2)
        self._barcode_length = 0
        self._b1s = []
        self._b2s = []
        self._b1s_rc = []
        self._b2s_rc = []

    def add_specimen(self, id, b1, b2):
        self._specimens.append((id, b1, b2))

        if b1 not in self._b1_specimens_index:
            self._b1_specimens_index[b1] = set()
        if b2 not in self._b2_specimens_index:
            self._b2_specimens_index[b2] = set()
        self._b1_specimens_index[b1].add(id)
        self._b2_specimens_index[b2].add(id)

        self._barcode_length = max(self._barcode_length, len(b1), len(b2))

    def validate(self):
        self._validate_barcodes_globally_unique()
        self._validate_barcode_lengths()
        self._cache_barcodes()

    def _validate_barcodes_globally_unique(self):
        dups = set(self.barcodes(Barcode.B1)).intersection(set(self.barcodes(Barcode.B2)))
        if len(dups) > 0:
            logging.warning(f"Duplicate Barcodes ({len(dups)}) in Fwd and Rev: {dups}")

    def _validate_barcode_lengths(self):
        if len(set(len(b) for b in self.barcodes(Barcode.B1))) > 1:
            logging.warning("Forward barcodes have inconsistent lengths")
        if len(set(len(b) for b in self.barcodes(Barcode.B2))) > 1:
            logging.warning("Reverse barcodes have inconsistent lengths")

    def _cache_barcodes(self):
        for b in self._b1_specimens_index.keys():
            self._b1s.append(b)
            self._b1s_rc.append(reverse_complement(b))
        for b in self._b2_specimens_index.keys():
            self._b2s.append(b)
            self._b2s_rc.append(reverse_complement(b))

    def specimens_for_barcodes(self, b1_list, b2_list):
        all = set()

        b1s = set()
        for b1 in b1_list: b1s.update(self._b1_specimens_index[b1])
        for b2 in b2_list: all.update(b1s.intersection(self._b2_specimens_index[b2]))

        return list(all)

    def b_length(self):
        return self._barcode_length

    def primer_length(self, which):
        return len(self.primer(which))

    def barcodes(self, which: Barcode, rc = False):
        if not rc:
            if which is Barcode.B1:
                return self._b1s
            elif which is Barcode.B2:
                return self._b2s
        else:
            if which is Barcode.B1:
                return self._b1s_rc
            elif which is Barcode.B2:
                return self._b2s_rc

    def primer(self, which: Primer, rc = False):
        if not rc:
            if which == Primer.P1:
                return self._p1
            elif which == Primer.P2:
                return self._p2
        else:
            if which == Primer.P1:
                return self._p1_rc
            elif which == Primer.P2:
                return self._p2_rc

class MatchParameters:
    def __init__(self, max_dist_primer, max_dist_index, search_len, ambiguity_threshold):
        self.max_dist_primer = max_dist_primer
        self.max_dist_index = max_dist_index
        self.search_len = search_len
        self.ambiguity_threshold = ambiguity_threshold

def homopolymer_compress(sequence):
    return ''.join(ch for ch, _ in itertools.groupby(sequence))

class ClusterDatabase:
    def __init__(self):
        self._next_cluster_id = 0
        self._clusters_for_group = {}
        self._cluster_representatives = {}
        self._cluster_representatives_rc = {}
        self._cluster_abundance = {}
        self._group_abundance = {}
        self._multimatch = 0
        self.match_untrimmed_start = 0
        self.match_untrimmed_end = 0
        self.singletons = 0
        self._group_assigned = 0

    def assign_best_cluster(self, group, sequence):
        self._group_abundance[group] = self._group_abundance.get(group, 0) + 1
        if group not in self._clusters_for_group:
            cluster_id = self._next_cluster_id
            self._next_cluster_id += 1
            self._clusters_for_group[group] = [cluster_id]
            self._cluster_representatives[cluster_id] = sequence
            self._cluster_abundance[cluster_id] = 1
            self.singletons += 1
            return cluster_id
        else:
            best_ed = -1
            best_c = -1
            for c in self._clusters_for_group[group]:
                rep = self._cluster_representatives[c]
                ed = len(rep)*0.05 #TODO 95% similarity
                m = align_seq(sequence, rep, ed, 0, len(sequence), mode=MODE_GLOBAL)
                if m.matched() and (best_ed == -1 or m.distance() < best_ed):
                    best_ed = m.distance()
                    best_c = c
            if best_c != -1:
                if self._cluster_abundance[best_c] == 1: self.singletons -= 1

                self._cluster_abundance[best_c] = self._cluster_abundance.get(best_c, 0) + 1
                return best_c
            else:
                cluster_id = self._next_cluster_id
                self._next_cluster_id += 1
                self._clusters_for_group[group].append(cluster_id)
                self._cluster_representatives[cluster_id] = sequence
                self._cluster_abundance[cluster_id] = 1
                self.singletons += 1
                return cluster_id

    def assign_best_group_and_cluster(self, groups, sequence, trimmed_start, trimmed_end):
        best_ed = -1
        best_c = -1
        best_g = -1
        best_l = (-1, -1)
        for g in groups:
            for c in self._clusters_for_group.get(g, []):
                pct_abd = self._cluster_abundance[c] / float(self._group_abundance[g])
                # if pct_abd < 0.1:  # TODO
                # if self._cluster_abundance[c] < 5 or pct_abd < 0.1:  # TODO
                if self._cluster_abundance[c] < 2: # TODO
                    continue

                rep = self._cluster_representatives[c]
                ed = len(rep) * 0.05 #TODO 95% similarity
                if trimmed_start and trimmed_end:
                    m = align_seq(sequence, rep, ed, 0, len(sequence), mode=MODE_GLOBAL)
                elif trimmed_start:
                    m = align_seq(rep, sequence, ed, 0, len(rep), mode=MODE_PREFIX)
                    if m.matched(): self.match_untrimmed_end += 1
                elif trimmed_end:
                    if c not in self._cluster_representatives_rc:
                        self._cluster_representatives_rc[c] = reverse_complement(rep)
                    rep = self._cluster_representatives_rc[c]
                    m = align_seq(rep, reverse_complement(sequence), ed, 0, len(rep), mode=MODE_PREFIX)
                    if m.matched(): self.match_untrimmed_start += 1
                else:
                    assert False
                if m.matched():
                    if best_ed == -1 or g == best_g:
                        best_ed = m.distance()
                        best_c = c
                        best_g = g
                        best_l = m.location()
                    else: #ambiguous
                        best_ed = -1
                        best_c = -1
                        best_g = -1
                        best_l = (-1, -1)
                        self._multimatch += 1
                        return best_g, best_c, best_l

        if best_c != -1:
            # Don't count non-barcode matches towards cluster abundance, to avoid reinforcing cycles
            # self._cluster_abundance[best_c] = self._cluster_abundance.get(best_c, 0) + 1
            self._group_assigned += 1
        return best_g, best_c, best_l

    def get_cluster_count(self):
        return self._next_cluster_id

    def multimatch(self):
        return self._multimatch

    def group_assigned(self):
        return self._group_assigned

def read_barcode_file(filename: str) -> Specimens:
    """
    Read a tab-separated barcode file and return a Specimens object.

    The file should have the following columns in order:
    sample_id, forward_barcode, forward_primer, reverse_barcode, reverse_primer
    A header line is optional and will be detected and skipped if present.

    Args:
    filename (str): Path to the barcode file

    Returns:
    Specimens: An object containing the parsed barcode information
    """
    specimens = None

    with open(filename, 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')

        # Read the first line
        first_line = next(reader, None)
        if not first_line:
            raise ValueError("The barcode file is empty")

        # Check if it's a header
        if not _is_valid_data_line(first_line):
            # If it's not valid data, assume it's a header and move to the next line
            first_line = next(reader, None)
            if not first_line:
                raise ValueError("The barcode file contains only a header")

        # Process the first line of data
        if len(first_line) != 5:
            raise ValueError(f"First data line does not have 5 columns: {first_line}")

        if not _is_valid_data_line(first_line):
            raise ValueError(f"Invalid data in first line: {first_line}")

        sample_id, forward_barcode, forward_primer, reverse_barcode, reverse_primer = first_line
        specimens = Specimens(forward_primer.upper(), reverse_primer.upper())
        specimens.add_specimen(sample_id, forward_barcode.upper(), reverse_barcode.upper())

        # Process the rest of the lines
        for row_num, row in enumerate(reader, start=2):
            if len(row) != 5:
                raise ValueError(f"Line {row_num} does not have 5 columns: {row}")

            if not _is_valid_data_line(row):
                raise ValueError(f"Invalid data in line {row_num}: {row}")

            sample_id, forward_barcode, forward_primer, reverse_barcode, reverse_primer = row
            specimens.add_specimen(sample_id, forward_barcode.upper(), reverse_barcode.upper())

    if not specimens:
        raise ValueError("No valid data found in the barcode file")

    return specimens


def _is_valid_data_line(line: List[str]) -> bool:
    """
    Check if a line contains valid data.

    This function checks if the barcodes and primers contain only DNA letters (A, T, G, C).

    Args:
    line (List[str]): A list of strings representing a line from the file

    Returns:
    bool: True if the line contains valid data, False otherwise
    """
    if len(line) != 5:
        return False

    _, forward_barcode, forward_primer, reverse_barcode, reverse_primer = line
    dna_sequences = [forward_barcode, forward_primer, reverse_barcode, reverse_primer]

    return all(set(seq.upper()).issubset({'A', 'T', 'G', 'C', 'Y', 'R', 'N', 'W', 'M', 'S', 'K'}) for seq in dna_sequences)

def open_sequence_file(filename, args):
    file_format = "fastq" if filename.endswith((".fastq", ".fq")) else "fasta"
    args.isfastq = file_format == "fastq"
    return SeqIO.parse(filename, file_format)

IUPAC_maps = [("Y", "C"), ("Y", "T"), ("R", "A"), ("R", "G"),
              ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T"),
              ("W", "A"), ("W", "T"), ("M", "A"), ("M", "C"),
              ("S", "C"), ("S", "G"), ("K", "G"), ("K", "T")]

def align_seq(query, target, max_distance, start, end, mode=MODE_INFIX):
    s = 0 if start == -1 else start
    e = len(target) if end == -1 else min(end, len(target))

    t = target[s:e]

    r = edlib.align(query, t, mode, 'locations', max_distance, additionalEqualities=IUPAC_maps)

    #edlib seems to return non-match with high edit distance sometimes.  Needs investigation
    if r['editDistance'] != -1 and r['editDistance'] > max_distance:
        r['editDistance'] = -1

    m = MatchResult(r)
    m.adjust_start(s)
    return m

def get_quality_seq(seq):
    if "phred_quality" in seq.letter_annotations:
        return seq.letter_annotations["phred_quality"]
    else:
        return [40]*len(seq)

ARG_TRIM_PRIMERS = "primers"
ARG_TRIM_BARCODES = "barcodes"
ARG_TRIM_TAILS = "tails"
ARG_TRIM_NONE = "none"

def create_write_operation(sample_id, args, seq, match, specimens):
    formatted_seq = seq.seq
    quality_scores = get_quality_seq(seq)

    s = 0
    e = len(formatted_seq)

    if args.trim == ARG_TRIM_PRIMERS:
        (s, e) = match.interprimer_extent()
    elif args.trim == ARG_TRIM_BARCODES:
        (s, e) = match.interbarcode_extent()
    elif args.trim == ARG_TRIM_TAILS:
        (s, e) = match.intertail_extent()

    if args.trim != ARG_TRIM_NONE:
        formatted_seq = seq.seq[s:e]
        quality_scores = quality_scores[s:e]
        match.trim_locations(s)

    return WriteOperation(
        sample_id=sample_id,
        seq_id=seq.id,
        distance_code=match.distance_code(),
        sequence=str(formatted_seq),
        quality_scores=quality_scores,
        match=match
    )



CLASS_SHORT_SEQ = "Sequence Too Short"
CLASS_LONG_SEQ = "Sequence Too Long"
CLASS_ORIENTATION_FAILED = "Could Not Determine Orientation"
CLASS_NO_BARCODES = "No Barcode Matches"
CLASS_NO_FWD_BARCODE = "No Forward Barcode Matches"
CLASS_FWD_BARCODE_TRUNCATED = "No Forward Barcode Matches (May be truncated)"
CLASS_NO_REV_BARCODE = "No Reverse Barcode Matches"
CLASS_REV_BARCODE_TRUNCATED = "No Reverse Barcode Matches (May be truncated)"
CLASS_NO_PRIMERS = "No Primer Matches"
CLASS_NO_FWD_PRIMER = "No Forward Primer Matches"
CLASS_NO_REV_PRIMER = "No Reverse Primer Matches"
CLASS_AMBIGUOUS_BARCODES = "Multiple Matches for Both Barcodes"
CLASS_AMBIGUOUS_FWD_BARCODE = "Multiple Matches for Forward Barcode"
CLASS_AMBIGUOUS_REV_BARCODE = "Multiple Matches for Reverse Barcode"
CLASS_MATCHED = "Matched"
CLASS_CLUSTER_MATCH = "Matched Cluster"

def process_sequences(seq_records, parameters, specimens, clusters, args):
    classifications = Counter()
    write_ops = []

    for seq in seq_records:
        sample_id = SAMPLE_ID_UNKNOWN

        classification = None
        match = SequenceMatch(len(seq), parameters.ambiguity_threshold, specimens.b_length())

        if args.min_length != -1 and len(seq) < args.min_length:
            classification = CLASS_SHORT_SEQ
        elif args.max_length != -1 and len(seq) > args.max_length:
            classification = CLASS_LONG_SEQ
        else:
            (oriented, seq, rseq) = determine_orientation(parameters, seq, specimens)
            if not oriented:
                classification = CLASS_ORIENTATION_FAILED
            else:
                match, seq = match_sequence(args, parameters, seq, rseq, specimens)
                sample_id = match_sample(match, sample_id, specimens)
                classification = classify_match(match, sample_id, specimens)

                classification, sample_id, seq = cluster_sequence(args, classification, clusters, match, sample_id, seq, specimens)

        classifications[classification] += 1
        if args.debug:
            logging.debug(f"{classification}")
        write_op = create_write_operation(sample_id, args, seq, match, specimens)
        write_ops.append(write_op)

    return write_ops, classifications


def cluster_sequence(args, classification, clusters, match, sample_id, seq, specimens):
    if args.cluster_sequences:
        assert clusters is not None

        (s, e) = match.interprimer_extent()
        ts = match.p1_match is not None
        te = match.p2_match is not None
        g = -1
        c = -1
        l = (-1, -1)
        if sample_id not in [SAMPLE_ID_UNKNOWN, SAMPLE_ID_AMBIGUOUS]:
            c = clusters.assign_best_cluster(sample_id, seq[s:e])
        elif sample_id == SAMPLE_ID_AMBIGUOUS:
            sp = specimens.specimens_for_barcodes(match.best_b1(), match.best_b2())
            g, c, l = clusters.assign_best_group_and_cluster(sp, seq[s:e], ts, te)
        elif match.p1_match and match.has_b1_match() and len(match.best_b1_match()) == 1:
            sp = specimens.specimens_for_barcodes(match.best_b1(), specimens.barcodes(Barcode.B2))
            g, c, l = clusters.assign_best_group_and_cluster(sp, seq[s:e], ts, te)
        elif match.p2_match and match.has_b2_match() and len(match.best_b2_match()) == 1:
            sp = specimens.specimens_for_barcodes(specimens.barcodes(Barcode.B1), match.best_b2())
            g, c, l = clusters.assign_best_group_and_cluster(sp, seq[s:e], ts, te)

        if sample_id in [SAMPLE_ID_UNKNOWN, SAMPLE_ID_AMBIGUOUS] and g != -1:
            sample_id = g
            classification = CLASS_CLUSTER_MATCH
            seq = seq[l[0]:l[1] + 1]  # TODO - cleaner option
    return classification, sample_id, seq


def classify_match(match, sample_id, specimens):
    classification = None
    if sample_id == SAMPLE_ID_UNKNOWN:
        if match.p1_match and match.p2_match:
            if not match.has_b1_match() and not match.has_b2_match():
                classification = CLASS_NO_BARCODES
            elif not match.has_b2_match():
                if (match.sequence_length - match.p2_match.location()[1]) < (specimens.b_length()):
                    classification = CLASS_REV_BARCODE_TRUNCATED
                else:
                    classification = CLASS_NO_REV_BARCODE
            elif not match.has_b1_match():
                if (match.p1_match.location()[0] - 1) < (specimens.b_length()):
                    classification = CLASS_FWD_BARCODE_TRUNCATED
                else:
                    classification = CLASS_NO_FWD_BARCODE
        else:
            if not match.p1_match and not match.p2_match:
                classification = CLASS_NO_PRIMERS
            elif not match.p2_match:
                classification = CLASS_NO_REV_PRIMER
            elif not match.p1_match:
                classification = CLASS_NO_FWD_PRIMER
            else:
                assert(False)
    elif sample_id == SAMPLE_ID_AMBIGUOUS:
        if len(match.best_b1()) > 1 and len(match.best_b2()) > 1:
            classification = CLASS_AMBIGUOUS_BARCODES
        elif len(match.best_b1()) > 1:
            classification = CLASS_AMBIGUOUS_FWD_BARCODE
        elif len(match.best_b2()) > 1:
            classification = CLASS_AMBIGUOUS_REV_BARCODE
        else:
            assert(False)
    else:
        classification = CLASS_MATCHED
    return classification

def match_sample(match, sample_id, specimens):
    if match.p1_match and match.p2_match and match.has_b1_match() and match.has_b2_match():
        ids = specimens.specimens_for_barcodes(match.best_b1(), match.best_b2())
        if len(ids) > 1:
            sample_id = SAMPLE_ID_AMBIGUOUS
        elif len(ids) == 1:
            sample_id = ids[0]
        else:
            logging.warning(f"No Specimens for combo: ({match.best_b1()}, {match.best_b2()})")
    return sample_id


def match_sequence(args, parameters, seq, rseq, specimens):
    match = SequenceMatch(len(seq), parameters.ambiguity_threshold, specimens.b_length())

    match_one_end(match, parameters, specimens, rseq, True, Primer.P1, Barcode.B1)

    match_one_end(match, parameters, specimens, seq, False, Primer.P2, Barcode.B2)

    return match, seq


def match_one_end(match, parameters, specimens, sequence, reversed_sequence, which_primer, which_barcode):

    primer_match = align_seq(specimens.primer(which_primer, True), sequence, parameters.max_dist_primer,
                             len(sequence) - parameters.search_len, len(sequence))
    if primer_match.matched():
        for b, b_rc in zip(specimens.barcodes(which_barcode, False), specimens.barcodes(which_barcode, True)):
            bd = None
            bm = None
            bc = None
            for l in primer_match.locations():
                barcode_match = align_seq(b_rc, sequence, parameters.max_dist_index, l[1] + 1,
                               len(sequence), MODE_PREFIX)
                if barcode_match.matched():
                    if bd is None or barcode_match.distance() < bd:
                        bm = barcode_match
                        bd = barcode_match.distance()
                        bc = b
            if bm:
                match.add_barcode_match(bm, bc, reversed_sequence, which_barcode)
        match.set_primer_match(primer_match, specimens.primer(which_primer, False), reversed_sequence, which_primer)


def determine_orientation(parameters, seq, specimens):
    rseq = seq.reverse_complement()
    rseq.id = seq.id
    rseq.description = seq.description

    p1 = align_seq(specimens.primer(Primer.P1), seq, specimens.primer_length(Primer.P1) * 0.5, 0, parameters.search_len)
    p1r = align_seq(specimens.primer(Primer.P1), rseq, specimens.primer_length(Primer.P1) * 0.5, 0, parameters.search_len)

    p2 = align_seq(specimens.primer(Primer.P2), seq, specimens.primer_length(Primer.P2) * 0.5, 0, parameters.search_len)
    p2r = align_seq(specimens.primer(Primer.P2), rseq, specimens.primer_length(Primer.P2) * 0.5, 0, parameters.search_len)

    vals = []
    if p1.matched(): vals.append(p1.distance())
    if p2.matched(): vals.append(p2.distance())
    if p1r.matched(): vals.append(p1r.distance())
    if p2r.matched(): vals.append(p2r.distance())

    if len(vals) == 0:
        logging.debug(f"Could not orient sequence - no primer matches")
        return False, seq, rseq
    else:
        best = min(vals)
        if (p1.distance() == best or p2r.distance() == best) and (not p2.matched() or p2.distance() > best) and (not p1r.matched() or p1r.distance() > best):
            reverse = False
        elif (p1r.distance() == best or p2.distance() == best) and (not p2r.matched() or p2r.distance() > best) and (not p1.matched() or p1.distance() > best):
            reverse = True
        else:
            logging.debug(f"Could not orient sequence: {p1.distance()} {p2.distance()} {p1r.distance()} {p2r.distance()} {len(seq)}")
            return False, seq, rseq

    if reverse:
        return True, rseq, seq
    return True, seq, rseq

def output_write_operation(write_op, output_manager, args):
    if args.output_to_files:
        fh = output_manager.get_file(write_op.sample_id)
    else:
        fh = sys.stdout

    formatted_seq = write_op.sequence
    if args.color:
        formatted_seq = color_sequence(formatted_seq, write_op.match, write_op.quality_scores)

    # Use '@' for FASTQ and '>' for FASTA
    header_symbol = '@' if args.isfastq else '>'
    fh.write(f"{header_symbol}{write_op.seq_id} {write_op.distance_code} {write_op.sample_id}\n")
    fh.write(f"{formatted_seq}\n")
    if args.isfastq:
        fh.write("+\n")
        fh.write("".join(chr(q+33) for q in write_op.quality_scores) + "\n")

def color_sequence(seq, match, quality_scores):
    blue = "\033[0;34m"
    green = "\033[0;32m"
    red = "\033[0;31m"
    NC = "\033[0m"  # No Color (reset)

    seq_len = len(seq)
    colored_seq = [''] * seq_len  # Initialize a list to hold colored characters
    start = 0
    end = seq_len

    def color_region(location, color):
        if location is not None:
            start, end = location
            if start < 0 or end < 0:
                return

            for i in range(start, end + 1):  # Include the end position
                if i < seq_len:
                    if quality_scores[i] < 10:
                        colored_seq[i] = color + seq[i].lower() + NC
                    else:
                        colored_seq[i] = color + seq[i] + NC

    # Color barcode1 (blue)
    color_region(match.get_barcode1_location(), blue)

    # Color primer1 (green)
    color_region(match.get_p1_location(), green)

    # Color primer2 (green)
    color_region(match.get_p2_location(), green)

    # Color barcode2 (blue)
    color_region(match.get_barcode2_location(), blue)

    # Fill in uncolored regions
    for i in range(seq_len):
        if colored_seq[i] == '':
            if quality_scores[i] < 10:
                colored_seq[i] = seq[i].lower()
            else:
                colored_seq[i] = seq[i]

    return ''.join(colored_seq[start:end])


class OutputManager:
    def __init__(self, output_dir, prefix, is_fastq):
        self.output_dir = output_dir
        self.prefix = prefix
        self.is_fastq = is_fastq
        self.file_stack = ExitStack()
        self.file_handles = {}

    def __enter__(self):
        os.makedirs(self.output_dir, exist_ok=True)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file_stack.close()

    def get_file(self, sample_id):
        if sample_id not in self.file_handles:
            filename = self._make_filename(sample_id)
            file = self.file_stack.enter_context(open(filename, 'w'))
            self.file_handles[sample_id] = file
        return self.file_handles[sample_id]

    def _make_filename(self, sample_id):
        safe_id = "".join(c if c.isalnum() or c in "._-$#" else "_" for c in sample_id)
        extension = '.fastq' if self.is_fastq else '.fasta'
        return os.path.join(self.output_dir, f"{self.prefix}{safe_id}{extension}")

def version():
    # 0.1 September 14, 2024 - rewritten from minibar.py
    return "specimux.py version 0.1"

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Specimux: Demultiplex MinION sequence by dual barcode indexes and primers.")
    
    parser.add_argument("barcode_file", help="File containing barcode information")
    parser.add_argument("sequence_file", help="Sequence file in Fasta or Fastq format, gzipped or plain text")

    parser.add_argument("--min-length", type=int, default=-1, help="Minimum sequence length.  Shorter sequences will be skipped (default: 400)")
    parser.add_argument("--max-length", type=int, default=-1, help="Maximum sequence length.  Longer sequences will be skipped (default: 1400)")
    parser.add_argument("-n", "--num-seqs", type=str, default="-1", help="Number of sequences to read from file (e.g., -n 100 or -n 102,3)")
    parser.add_argument("-p", "--percent-match", type=float, default=0.75, help="Percentage match (default: 0.75)")
    parser.add_argument("-e", "--index-edit-distance", type=int, default=-1, help="Barcode edit distance value, overrides -p")
    parser.add_argument("-E", "--primer-edit-distance", type=int, default=-1, help="Primer edit distance value, overrides -p")
    parser.add_argument("-l", "--search-len", type=int, default=80, help="Length to search for index and primer at start and end of sequence (default: 80)")
    parser.add_argument("-A", "--ambiguity-threshold",type=float, default=1.0, help="Threshold for considering the edit distance between two barcodes to be different")
    parser.add_argument("--cluster-sequences", action="store_true", help="Cluster sequences based on alignment and resolve ambiguities based on cluster identity")
    parser.add_argument("-F", "--output-to-files", action="store_true", help="Create individual sample files for sequences")
    parser.add_argument("-P", "--output-file-prefix", default="sample_", help="Prefix for individual files when using -F (default: sample_)")
    parser.add_argument("-D", "--output-dir", default=".", help="Directory for individual files when using -F (default: .)")
    parser.add_argument("--color", action="store_true", help="Highlight barcode matches in blue, primer matches in green")
    parser.add_argument("--trim", choices=[ARG_TRIM_NONE, ARG_TRIM_TAILS, ARG_TRIM_BARCODES, ARG_TRIM_PRIMERS], default=ARG_TRIM_BARCODES, help="trimming to apply")
    parser.add_argument("--diagnostics", action="store_true", help="Output extra diagnostics")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")

    parser.add_argument("-t", "--threads", type=int, default=-1, help="Number of worker threads to use")
    parser.add_argument("-v", "--version", action="version", version=version())

    args = parser.parse_args(argv[1:])

    if args.num_seqs:
        process_num_seqs(args, parser)

    return args

def process_num_seqs(args, parser):
    if ',' in args.num_seqs:
        start, num = args.num_seqs.split(',')
        try:
            args.start_seq = int(start)
            args.num_seqs = int(num)
        except ValueError:
            parser.error("Invalid format for -n option. Use 'start,num' with integers.")
    else:
        try:
            args.num_seqs = int(args.num_seqs)
            args.start_seq = 1
        except ValueError:
            parser.error("Invalid format for -n option. Use an integer or 'start,num' with integers.")


class WriteOperation(NamedTuple):
    sample_id: str
    seq_id: str
    distance_code: str
    sequence: str
    quality_scores: List[int]
    match: SequenceMatch

class WorkItem(NamedTuple):
    seq_number: int
    seq_records: List  # This now contains actual sequence records, not an iterator
    parameters: MatchParameters
    specimens: Specimens
    args: argparse.Namespace

def worker(work_item, specimens, args):
    write_ops, classifications = process_sequences(work_item.seq_records, work_item.parameters,
                                                   specimens, None, args)
    return write_ops, classifications


def specimux_mp(args):
    specimens = read_barcode_file(args.barcode_file)
    specimens.validate()
    parameters = setup_match_parameters(args, specimens)

    seq_records = open_sequence_file(args.sequence_file, args)

    start_time = timeit.default_timer()

    sequence_block_size = 1000
    last_seq_to_output = args.num_seqs
    all_seqs = last_seq_to_output < 0

    # Skip to the start_seq if necessary
    if args.start_seq > 1:
        for _ in itertools.islice(seq_records, args.start_seq - 1):
            pass

    classifications = Counter()

    if args.cluster_sequences:
        logging.error(f"--cluster-sequences must be used in single-threaded mode (-t 1)")

    with OutputManager(args.output_dir, args.output_file_prefix, args.isfastq) as output_manager:
        if args.threads == -1:
            num_processes = multiprocessing.cpu_count()
        else:
            num_processes = args.threads

        logging.info(f"Will run {num_processes} worker processes")


        with Pool(processes=num_processes) as pool:
            worker_func = partial(worker, specimens=specimens, args=args)

            def producer():
                num_seqs = 0
                while all_seqs or num_seqs < last_seq_to_output:
                    to_read = sequence_block_size if all_seqs else min(sequence_block_size, last_seq_to_output - num_seqs)
                    seq_batch = list(itertools.islice(seq_records, to_read))
                    if not seq_batch:
                        break
                    work_item = WorkItem(num_seqs, seq_batch, parameters, None, args)
                    yield work_item
                    num_seqs += len(seq_batch)

            for write_ops, batch_classifications in pool.imap_unordered(worker_func, producer()):
                classifications += batch_classifications
                for write_op in write_ops:
                    output_write_operation(write_op, output_manager, args)

                total_processed = sum(classifications.values())
                matched = classifications[CLASS_MATCHED]
                pct = matched/total_processed if total_processed > 0 else 0
                print(f"read: {total_processed}\t\tmatched: {matched}\t\t{pct:.2%}", end='\r', file=sys.stderr)
                sys.stderr.flush()

    elapsed = timeit.default_timer() - start_time

    print()
    logging.info(f"Elapsed time: : {elapsed:.2f} seconds")

    output_diagnostics(args, classifications)

def specimux(args):
    specimens = read_barcode_file(args.barcode_file)
    specimens.validate()
    parameters = setup_match_parameters(args, specimens)

    seq_records = open_sequence_file(args.sequence_file, args)

    start_time = timeit.default_timer()

    sequence_block_size = 1000
    last_seq_to_output = args.num_seqs
    all_seqs = last_seq_to_output < 0

    # Skip to the start_seq if necessary
    if args.start_seq > 1:
        for _ in itertools.islice(seq_records, args.start_seq - 1):
            pass

    classifications = Counter()
    clusters = ClusterDatabase()

    with OutputManager(args.output_dir, args.output_file_prefix, args.isfastq) as output_manager:
        num_seqs = 0
        while all_seqs or num_seqs < last_seq_to_output:
            to_read = sequence_block_size if all_seqs else min(sequence_block_size, last_seq_to_output - num_seqs)
            seq_batch = list(itertools.islice(seq_records, to_read))
            if not seq_batch:
                break
            write_ops, batch_classifications = process_sequences(seq_batch, parameters, specimens, clusters, args)
            classifications += batch_classifications
            for write_op in write_ops:
                output_write_operation(write_op, output_manager, args)

            num_seqs += len(seq_batch)
            matched = classifications[CLASS_MATCHED]
            cmatched = classifications[CLASS_CLUSTER_MATCH]
            pct = (matched+cmatched)/num_seqs if num_seqs > 0 else 0
            print(f"read: {num_seqs:>7}\t\tmatched: {matched:>7}/{cmatched:>7}\t\t{pct:.2%}\t-{clusters.multimatch()}\tclusters: {clusters.get_cluster_count()} singletons: {clusters.singletons}", end='\r', file=sys.stderr)
            sys.stderr.flush()

    elapsed = timeit.default_timer() - start_time

    print()
    logging.info(f"Elapsed time: : {elapsed:.2f} seconds")

    output_diagnostics(args, classifications)

def output_diagnostics(args, classifications):
    if args.diagnostics:
        # Sort the classifications by their counts in descending order
        sorted_classifications = sorted(
            classifications.items(),
            key=lambda x: x[1],
            reverse=True
        )

        # Find the length of the longest classification name for padding
        max_name_length = max(len(name) for name, _ in sorted_classifications)
        # Determine the width needed for the count column
        max_count_length = len(str(max(count for _, count in sorted_classifications)))

        total = sum(classifications.values())
        # Print the sorted classifications with aligned columns
        logging.info("Classification Statistics:")
        for classification, count in sorted_classifications:
            logging.info(f"{classification:<{max_name_length}} : {count:>{max_count_length}} ({count / total:2.2%})")


def setup_match_parameters(args, specimens):
    def _calculate_distances(sequences):
        distances = []
        for seq1, seq2 in itertools.combinations(sequences, 2):
            dist = edlib.align(seq1, seq2, task="distance")["editDistance"]
            distances.append(dist)
        return min(distances), Counter(distances)

    def _sanity_check_distance(sequences, max_distance, desc):
        if len(sequences) <= 1:
            return
        m, c = _calculate_distances(sequences)
        logging.info(f"Minimum edit distance is {m} for {desc}")

    max_search_area = args.search_len
    max_dist_index = specimens.b_length() - int(specimens.b_length() * args.percent_match)
    primer_length = min(specimens.primer_length(Primer.P1), specimens.primer_length(Primer.P2))
    max_dist_primer = primer_length - int(primer_length * args.percent_match)
    if args.index_edit_distance != -1: max_dist_index = args.index_edit_distance
    if args.primer_edit_distance != -1: max_dist_primer = args.primer_edit_distance

    if args.diagnostics:
        logging.info(f"Maximum edit distances set to {max_dist_index} (barcode) and {max_dist_primer} (primer)")
        _sanity_check_distance(specimens.barcodes(Barcode.B1), max_dist_index, "Forward Barcodes")
        _sanity_check_distance(specimens.barcodes(Barcode.B2), max_dist_index, "Reverse Barcodes")

        combined_barcodes = list(specimens.barcodes(Barcode.B1)) + [reverse_complement(b) for b in specimens.barcodes(Barcode.B2)]
        _sanity_check_distance(combined_barcodes, max_dist_index, "Forward Barcodes + Reverse Complement of Reverse Barcodes")

        primers = [specimens.primer(Primer.P1, False), specimens.primer(Primer.P1, True),
                   specimens.primer(Primer.P2, False), specimens.primer(Primer.P2, True)]
        combined_barcodes_and_primers = combined_barcodes + primers
        _sanity_check_distance(combined_barcodes_and_primers, min(max_dist_index, max_dist_primer),
                               "Forward Barcodes + Reverse Complement of Reverse Barcodes + Both Primers")

    parameters = MatchParameters(max_dist_primer, max_dist_index, max_search_area, args.ambiguity_threshold)
    return parameters

def main(argv):
    args = parse_args(argv)

    l = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=l, format='%(asctime)s - %(levelname)s - %(message)s')

    if args.threads == 1:
        specimux(args)
    else:
        specimux_mp(args)

if __name__ == "__main__":
    main(sys.argv)

