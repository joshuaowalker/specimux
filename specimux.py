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
import math
import traceback
from _operator import itemgetter

import edlib
from collections import Counter
from contextlib import ExitStack
from functools import partial
from operator import itemgetter
from typing import List, Tuple
from typing import NamedTuple
from typing import Optional, Union
from enum import Enum

from Bio import SeqIO
from Bio.Seq import reverse_complement
from multiprocessing import Pool
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


IUPAC_EQUIV = [("Y", "C"), ("Y", "T"), ("R", "A"), ("R", "G"),
               ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T"),
               ("W", "A"), ("W", "T"), ("M", "A"), ("M", "C"),
               ("S", "C"), ("S", "G"), ("K", "G"), ("K", "T"),
               ("B", "C"), ("B", "G"), ("B", "T"),
               ("D", "A"), ("D", "G"), ("D", "T"),
               ("H", "A"), ("H", "C"), ("H", "T"),
               ("V", "A"), ("V", "C"), ("V", "G"), ]

IUPAC_CODES = set([t[0] for t in IUPAC_EQUIV])
IUPAC_CODES.update(["A", "C", "G", "T"])

class AlignMode:
    GLOBAL = 'NW'
    INFIX = 'HW'
    PREFIX = 'SHW'

class SampleId:
    AMBIGUOUS = "ambiguous"
    UNKNOWN = "unknown"
    PREFIX_FWD_MATCH = "barcode_fwd_"
    PREFIX_REV_MATCH = "barcode_rev_"

class TrimMode:
    PRIMERS = "primers"
    BARCODES = "barcodes"
    TAILS = "tails"
    NONE = "none"

class MatchCode:
    SHORT_SEQ = "1A: Sequence Too Short"
    LONG_SEQ = "1B: Sequence Too Long"
    ORIENTATION_FAILED = " 2: Could Not Determine Orientation"
    NO_PRIMERS = "3C: No Primer Matches"
    NO_FWD_PRIMER = "3A: No Forward Primer Matches"
    NO_REV_PRIMER = "3B: No Reverse Primer Matches"
    NO_BARCODES = "4E: No Barcode Matches"
    NO_FWD_BARCODE = "4A: No Forward Barcode Matches"
    FWD_BARCODE_TRUNCATED = "4B: No Forward Barcode Matches (May be truncated)"
    NO_REV_BARCODE = "4C: No Reverse Barcode Matches"
    REV_BARCODE_TRUNCATED = "4D: No Reverse Barcode Matches (May be truncated)"
    MULTIPLE_PRIMER_PAIRS_MATCHED = "5: Multiple Primer Pair Full Matches"
    AMBIGUOUS_BARCODES = "6C: Multiple Matches for Both Barcodes"
    AMBIGUOUS_FWD_BARCODE = "6A: Multiple Matches for Forward Barcode"
    AMBIGUOUS_REV_BARCODE = "6B: Multiple Matches for Reverse Barcode"
    MATCHED = " 7: Matched"

class Barcode(Enum):
    B1 = 1
    B2 = 2

class Primer(Enum):
    P1 = 3
    P2 = 4

class MatchResult:
    """Wrapper around edlib match result"""
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

    def reversed(self, seq_length):
        """Return a new MatchResult with locations relative to the reversed sequence"""
        m = MatchResult(self._edlib_match.copy())
        m._reverse(seq_length)
        return m

    def _reverse(self, seq_length):
        """Update locations to be relative to the start of the reversed sequence"""
        if self.distance() == -1: return

        r = self._edlib_match
        r['locations'] = [(seq_length - loc[1] - 1, seq_length - loc[0] - 1) for loc in r['locations']]

    def adjust_start(self, s):
        """Update start in locations by s"""
        if self.distance() == -1: return

        self._edlib_match['locations'] = [(loc[0] + s, loc[1] + s) for loc in self._edlib_match['locations']]  # adjust relative match to absolute

class SequenceMatch:

    def __init__(self, sequence_length, barcode_length):
        self.sequence_length = sequence_length
        self.p1_match: Optional[MatchResult] = None
        self._b1_matches: List[Tuple[str, MatchResult, float]] = []
        self.p2_match: Optional[MatchResult] = None
        self._b2_matches: List[Tuple[str, MatchResult, float]] = []
        self._p1: Optional[str] = None
        self._p2: Optional[str] = None
        self.ambiguity_threshold = 1.0
        self._barcode_length = barcode_length

    def add_barcode_match(self, match: MatchResult, barcode: str, reverse: bool, which: Barcode):
        m = match
        if reverse:
            m = m.reversed(self.sequence_length)

        if which is Barcode.B1:
            self._b1_matches.append((barcode, m, m.distance()))
            self._b1_matches.sort(key=itemgetter(2))  # Sort by edit distance

        else:
            self._b2_matches.append((barcode, m, m.distance()))
            self._b2_matches.sort(key=itemgetter(2))  # Sort by edit distance

    def b1_distance(self):
        return self._b1_matches[0][2] if len(self._b1_matches) > 0 else -1

    def b2_distance(self):
        return self._b2_matches[0][2] if len(self._b2_matches) > 0 else -1

    def best_b1(self) -> List[str]:
        if not self._b1_matches:
            return []
        best_distance = self.b1_distance()
        return [b for b, _, d in self._b1_matches if abs(d - best_distance) < self.ambiguity_threshold]

    def best_b2(self) -> List[str]:
        if not self._b2_matches:
            return []
        best_distance = self.b2_distance()
        return [b for b, _, d in self._b2_matches if abs(d - best_distance) < self.ambiguity_threshold]

    def set_primer_match(self, match: MatchResult, primer: str, reverse: bool, which: Primer):
        m = match
        if reverse:
            m = m.reversed(self.sequence_length)
        if which is Primer.P1:
            self.p1_match = m
            self._p1 = primer
        else:
            self.p2_match = m
            self._p2 = primer

    def get_p1(self):
        return self._p1

    def get_p2(self):
        return self._p2

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


class PrimerPairInfo:
    def __init__(self, p1: str, p2: str):
        self.p1 = p1.upper()
        self.p2 = p2.upper()
        self.p1_rc = reverse_complement(p1.upper())
        self.p2_rc = reverse_complement(p2.upper())
        self.specimens = set()  # Set of specimen IDs using this primer pair
        self.b1s = set()  # Set of forward barcodes used with this primer pair
        self.b2s = set()  # Set of reverse barcodes used with this primer pair

    def __eq__(self, other):
        if not isinstance(other, PrimerPairInfo):
            return NotImplemented
        return self.p1 == other.p1 and self.p2 == other.p2

    def __hash__(self):
        return hash((self.p1, self.p2))

class Specimens:
    def __init__(self):
        self._specimens = []  # List of (id, b1, p1, b2, p2) tuples for reference
        self._barcode_length = 0
        self._primer_pairs = {}  # Maps (p1, p2) tuple to PrimerPairInfo
        self._specimen_ids = set()
        self._primer_specimen_index = {}  # Maps (p1,p2) to set of specimen IDs

    def add_specimen(self, specimen_id, b1, p1, b2, p2):
        """Add a specimen with its barcodes and primers"""
        if specimen_id in self._specimen_ids:
            raise ValueError(format(f"Duplicate specimen id in index file: {specimen_id}"))
        self._specimen_ids.add(specimen_id)

        self._specimens.append((specimen_id, b1, p1, b2, p2))
        self._barcode_length = max(self._barcode_length, len(b1), len(b2))

        primer_key = (p1.upper(), p2.upper())
        if primer_key not in self._primer_pairs:
            self._primer_pairs[primer_key] = PrimerPairInfo(p1, p2)

        primer_pair_info = self._primer_pairs[primer_key]
        primer_pair_info.specimens.add(specimen_id)
        primer_pair_info.b1s.add(b1.upper())  # Add forward barcode
        primer_pair_info.b2s.add(b2.upper())  # Add reverse barcode

        if primer_key not in self._primer_specimen_index:
            self._primer_specimen_index[primer_key] = set()
        self._primer_specimen_index[primer_key].add(specimen_id)

    def get_primer_pairs(self):
        """Get list of all unique PrimerPairInfo objects"""
        return list(self._primer_pairs.values())

    def validate(self):
        self._validate_barcodes_globally_unique()
        self._validate_barcode_lengths()
        logging.info(f"Number of unique primer pairs: {len(self._primer_pairs)}")
        for pair in self._primer_pairs.values():
            logging.info(f"Primer pair {pair.p1}/{pair.p2}: {len(pair.specimens)} specimens")

    def _validate_barcodes_globally_unique(self):
        all_b1s = set()
        all_b2s = set()
        for pair in self._primer_pairs.values():
            all_b1s.update(pair.b1s)
            all_b2s.update(pair.b2s)
        dups = all_b1s.intersection(all_b2s)
        if len(dups) > 0:
            logging.warning(f"Duplicate Barcodes ({len(dups)}) in Fwd and Rev: {dups}")

    def _validate_barcode_lengths(self):
        all_b1s = set()
        all_b2s = set()
        for pair in self._primer_pairs.values():
            all_b1s.update(pair.b1s)
            all_b2s.update(pair.b2s)
        if len(set(len(b) for b in all_b1s)) > 1:
            logging.warning("Forward barcodes have inconsistent lengths")
        if len(set(len(b) for b in all_b2s)) > 1:
            logging.warning("Reverse barcodes have inconsistent lengths")

    def specimens_for_barcodes_and_primers(self, b1_list, b2_list, p1_matched, p2_matched):
        matching_specimens = []
        for spec_id, b1, p1, b2, p2 in self._specimens:
            if (p1_matched.upper() == p1 and
                    p2_matched.upper() == p2 and
                    b1.upper() in b1_list and
                    b2.upper() in b2_list):
                matching_specimens.append(spec_id)

        return matching_specimens

    def b_length(self):
        return self._barcode_length

class MatchParameters:
    def __init__(self, max_dist_primer, max_dist_index, search_len):
        self.max_dist_primer = max_dist_primer
        self.max_dist_index = max_dist_index
        self.search_len = search_len

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


import os
from typing import Dict, List, Optional
from cachetools import LRUCache
from collections import defaultdict
import logging
from contextlib import contextmanager


class FileHandleCache(LRUCache):
    """Custom LRU cache that works with CachedFileManager to close file handles on eviction."""

    def __init__(self, maxsize, file_manager):
        super().__init__(maxsize)
        self.file_manager = file_manager

    def popitem(self):
        """Override popitem to ensure file handle cleanup on eviction."""
        key, file_handle = super().popitem()
        # Ensure buffer is flushed before closing
        self.file_manager.flush_buffer(key)
        file_handle.close()
        return key, file_handle

class CachedFileManager:
    """Manages a pool of file handles with LRU caching and write buffering."""

    def __init__(self, max_open_files: int, buffer_size: int):
        """
        Initialize the file manager.

        Args:
            max_open_files: Maximum number of files to keep open at once
            buffer_size: Number of sequences to buffer before writing to disk
        """
        self.max_open_files = max_open_files
        self.buffer_size = buffer_size

        self.file_cache = FileHandleCache(maxsize=max_open_files, file_manager=self)

        self.write_buffers: Dict[str, List[str]] = defaultdict(list)

        # Track if file has been opened before (for append vs truncate)
        self.files_opened: Dict[str, bool] = {}

        # Count of buffer flushes for debugging
        self.flush_count = 0

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Ensure all buffers are flushed and files are closed on exit."""
        try:
            self.flush_all()
        finally:
            self.close_all()

    @contextmanager
    def get_file(self, filename: str) -> Optional[object]:
        """
        Get a file handle for writing, managing the cache as needed.

        Args:
            filename: Path to the file to open

        Returns:
            A context-managed file handle
        """
        try:
            f = self.file_cache[filename]
        except KeyError:
            # If file handle isn't cached, we need to open it.
            # We open in write mode the first time, and then subsequently append
            mode = 'a' if self.files_opened.get(filename, False) else 'w'
            try:
                f = open(filename, mode)
                self.files_opened[filename] = True
                # Add to cache, which may trigger eviction of least recently used file
                self.file_cache[filename] = f
            except OSError as e:
                logging.error(f"Failed to open file {filename}: {e}")
                raise

        yield f

    def write(self, filename: str, data: str):
        """
        Write data to a file through the buffer system.

        Args:
            filename: File to write to
            data: String data to write
        """
        self.write_buffers[filename].append(data)

        if len(self.write_buffers[filename]) >= self.buffer_size:
            self.flush_buffer(filename)

    def flush_buffer(self, filename: str):
        """
        Flush the buffer for a specific file to disk.

        Args:
            filename: File whose buffer to flush
        """
        if not self.write_buffers[filename]:
            return

        buffer_data = ''.join(self.write_buffers[filename])
        with self.get_file(filename) as f:
            f.write(buffer_data)
            f.flush()  # Ensure data is written to OS

        self.write_buffers[filename].clear()
        self.flush_count += 1

    def flush_all(self):
        """Flush all buffers to disk."""
        for filename in list(self.write_buffers.keys()):
            self.flush_buffer(filename)

    def close_all(self):
        """Close all open files."""
        for f in self.file_cache.values():
            try:
                f.close()
            except Exception as e:
                logging.warning(f"Error closing file: {e}")
        self.file_cache.clear()

    def close_file(self, filename: str):
        """
        Close a specific file if it's open.

        Args:
            filename: File to close
        """
        if filename in self.file_cache:
            try:
                self.file_cache[filename].close()
                del self.file_cache[filename]
            except Exception as e:
                logging.warning(f"Error closing file {filename}: {e}")


class OutputManager:
    """Manages output files with caching and buffering."""

    def __init__(self, output_dir: str, prefix: str, is_fastq: bool,
                 max_open_files: int = 200, buffer_size: int = 500):
        self.output_dir = output_dir
        self.prefix = prefix
        self.is_fastq = is_fastq
        self.file_manager = CachedFileManager(max_open_files, buffer_size)

    def __enter__(self):
        os.makedirs(self.output_dir, exist_ok=True)
        self.file_manager.__enter__()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return self.file_manager.__exit__(exc_type, exc_val, exc_tb)

    def _make_filename(self, sample_id: str) -> str:
        """Create a filename for a sample ID, ensuring it's safe for the filesystem."""
        safe_id = "".join(c if c.isalnum() or c in "._-$#" else "_" for c in sample_id)
        extension = '.fastq' if self.is_fastq else '.fasta'
        return os.path.join(self.output_dir, f"{self.prefix}{safe_id}{extension}")

    def write_sequence(self, sample_id: str, header: str, sequence: str,
                       quality_scores: Optional[List[int]] = None):
        """Write a sequence to the appropriate output file."""
        filename = self._make_filename(sample_id)

        # Build the output string
        output = []
        output.append(f"{'@' if self.is_fastq else '>'}{header}\n")
        output.append(f"{sequence}\n")
        if self.is_fastq:
            output.append("+\n")
            output.append("".join(chr(q + 33) for q in quality_scores) + "\n")

        self.file_manager.write(filename, ''.join(output))

class WorkerException(Exception):
    pass

def read_barcode_file(filename: str) -> Specimens:
    """
    Read a tab-separated barcode file and return a Specimens object.
    Each line contains: sample_id, forward_barcode, forward_primer, reverse_barcode, reverse_primer
    """
    specimens = Specimens()

    with open(filename, 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')

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
        specimens.add_specimen(sample_id, forward_barcode.upper(), forward_primer.upper(),
                               reverse_barcode.upper(), reverse_primer.upper())

        # Process the rest of the lines
        for row_num, row in enumerate(reader, start=2):
            if len(row) != 5:
                raise ValueError(f"Line {row_num} does not have 5 columns: {row}")

            if not _is_valid_data_line(row):
                raise ValueError(f"Invalid data in line {row_num}: {row}")

            sample_id, forward_barcode, forward_primer, reverse_barcode, reverse_primer = row
            specimens.add_specimen(sample_id, forward_barcode.upper(), forward_primer.upper(),
                                   reverse_barcode.upper(), reverse_primer.upper())

    if len(specimens._specimens) == 0:
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

    return all(set(seq.upper()).issubset(IUPAC_CODES) for seq in dna_sequences)

def open_sequence_file(filename, args):
    file_format = "fastq" if filename.endswith((".fastq", ".fq")) else "fasta"
    args.isfastq = file_format == "fastq"
    return SeqIO.parse(filename, file_format)

def align_seq(query: Union[str, Seq, SeqRecord],
             target: Union[str, Seq, SeqRecord],
             max_distance: int,
             start: int,
             end: int,
             mode: str = AlignMode.INFIX) -> MatchResult:
    # Convert query to string if it's a Seq or SeqRecord - not functionally necessary, but improves performance
    if isinstance(query, (Seq, SeqRecord)):
        query = str(query.seq if isinstance(query, SeqRecord) else query)

    # Extract target sequence string and handle slicing
    if isinstance(target, (Seq, SeqRecord)):
        target_seq = str(target.seq if isinstance(target, SeqRecord) else target)
    else:
        target_seq = target

    s = 0 if start == -1 else start
    e = len(target_seq) if end == -1 else min(end, len(target_seq))

    t = target_seq[s:e]

    r = edlib.align(query, t, mode, 'locations', max_distance, additionalEqualities=IUPAC_EQUIV)

    # Handle edlib sometimes returning non-match with high edit distance
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

def create_write_operation(sample_id, args, seq, match):
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

    return WriteOperation(
        sample_id=sample_id,
        seq_id=seq.id,
        distance_code=match.distance_code(),
        sequence=str(formatted_seq),
        quality_scores=quality_scores,
        match=match
    )



def process_sequences(seq_records: List[SeqRecord],
                     parameters: MatchParameters,
                     specimens: Specimens,
                     args: argparse.Namespace) -> Tuple[List[WriteOperation], Counter, Counter]:
    classifications = Counter()
    unmatched_barcodes = Counter()
    write_ops = []

    for seq in seq_records:
        sample_id = SampleId.UNKNOWN

        match = SequenceMatch(len(seq), specimens.b_length())

        if args.min_length != -1 and len(seq) < args.min_length:
            classification = MatchCode.SHORT_SEQ
        elif args.max_length != -1 and len(seq) > args.max_length:
            classification = MatchCode.LONG_SEQ
        else:
            (oriented, seq, rseq) = determine_orientation(parameters, seq, specimens)
            if not oriented:
                classification = MatchCode.ORIENTATION_FAILED
            else:
                #extract string versions for performance - roughly 18% improvement
                s = str(seq.seq)
                rs = str(rseq.seq)

                matches = match_sequence(parameters, s, rs, specimens)
                match, multimatch = choose_best_match(matches)
                if multimatch:
                    classification = MatchCode.MULTIPLE_PRIMER_PAIRS_MATCHED
                    sample_id = SampleId.AMBIGUOUS
                else:
                    sample_id = match_sample(match, sample_id, specimens)
                    classification = classify_match(match, sample_id, specimens)

                if args.group_unknowns:
                    sample_id = group_sample(match, sample_id)

                unmatched_barcode = analyze_barcode_region(seq.seq, match, specimens.b_length())
                if unmatched_barcode:
                    unmatched_barcodes[unmatched_barcode] += 1

        classifications[classification] += 1
        if args.debug:
            logging.debug(f"{classification}")
        write_op = create_write_operation(sample_id, args, seq, match)
        write_ops.append(write_op)

    return write_ops, classifications, unmatched_barcodes


def analyze_barcode_region(seq: str, match: SequenceMatch, bc_len: int) -> Optional[str]:
    if match.has_b1_match() and match.p2_match and not match.has_b2_match():
        start = match.p2_match.location()[1]
        end = min(len(seq), start+bc_len)
        bc_region = seq[start:end]
        if len(bc_region) < bc_len:
            bc_region += '-'*(bc_len-len(bc_region))
        return bc_region + "<"
    elif match.has_b2_match() and match.p1_match and not match.has_b1_match():
        end = match.p1_match.location()[0]
        start = max(0, end - bc_len)
        bc_region = seq[start:end]
        if len(bc_region) < bc_len:
            bc_region = '-'*(bc_len-len(bc_region)) + bc_region
        return ">"+bc_region
    else:
        return None


def classify_match(match: SequenceMatch, sample_id: str, specimens: Specimens) -> str:
    if sample_id == SampleId.UNKNOWN:
        if not match.p1_match and not match.p2_match:
            classification = MatchCode.NO_PRIMERS
        elif not match.has_b2_match() and match.has_b1_match():
            if match.p2_match:
                if (match.sequence_length - match.p2_match.location()[1]) < (specimens.b_length()):
                    classification = MatchCode.REV_BARCODE_TRUNCATED
                else:
                    classification = MatchCode.NO_REV_BARCODE
            else:
                classification = MatchCode.NO_REV_PRIMER
        elif not match.has_b1_match() and match.has_b2_match():
            if match.p1_match:
                if (match.p1_match.location()[0] - 1) < (specimens.b_length()):
                    classification = MatchCode.FWD_BARCODE_TRUNCATED
                else:
                    classification = MatchCode.NO_FWD_BARCODE
            else:
                classification = MatchCode.NO_FWD_PRIMER
        elif not match.p2_match:
            classification = MatchCode.NO_REV_PRIMER
        elif not match.p1_match:
            classification = MatchCode.NO_FWD_PRIMER
        elif not match.has_b1_match() and not match.has_b2_match():
            classification = MatchCode.NO_BARCODES
        else:
            raise RuntimeError(f"Unexpected unknown sample state: "
                               f"p1_match={match.p1_match is not None}, "
                               f"p2_match={match.p2_match is not None}, "
                               f"b1_match={match.has_b1_match()}, "
                               f"b2_match={match.has_b2_match()}")
    elif sample_id == SampleId.AMBIGUOUS:
        if len(match.best_b1()) > 1 and len(match.best_b2()) > 1:
            classification = MatchCode.AMBIGUOUS_BARCODES
        elif len(match.best_b1()) > 1:
            classification = MatchCode.AMBIGUOUS_FWD_BARCODE
        elif len(match.best_b2()) > 1:
            classification = MatchCode.AMBIGUOUS_REV_BARCODE
        else:
            raise RuntimeError(f"Unexpected ambiguous sample state: "
                               f"b1_matches={len(match.best_b1())}, "
                               f"b2_matches={len(match.best_b2())}")
    else:
        classification = MatchCode.MATCHED
    return str(classification)

def choose_best_match(matches: List[SequenceMatch]) -> Tuple[SequenceMatch, bool]:
    """In case there were matches under multiple primer pairs, choose the best"""
    best_score = -1
    best_match = None
    multimatch = False

    if not matches:
        raise ValueError("No matches provided to choose_best_match")

    for m in matches:
        score = 0
        if m.p1_match and m.p2_match and m.has_b1_match() and m.has_b2_match():
            score = 4
        elif m.p1_match and m.p2_match and (m.has_b1_match() or m.has_b2_match()):
            score = 3
        elif m.p1_match and m.p2_match:
            score = 2
        elif m.p1_match or m.p2_match:
            score = 1
        if score > best_score:
            best_score = score
            best_match = m
            multimatch = False
        elif score == best_score:
            multimatch = True

    if best_match is None:
        raise RuntimeError("Failed to select best match despite having matches")

    # ambiguity between primer pairs only really matters when there was a full match
    if best_score == 4 and multimatch:
        return best_match, True
    else:
        return best_match, False

def match_sample(match: SequenceMatch, sample_id: str, specimens: Specimens) -> str:
    if match.p1_match and match.p2_match and match.has_b1_match() and match.has_b2_match():
        ids = specimens.specimens_for_barcodes_and_primers(
            match.best_b1(), match.best_b2(), match.get_p1(), match.get_p2())
        if len(ids) > 1:
            sample_id = SampleId.AMBIGUOUS
        elif len(ids) == 1:
            sample_id = ids[0]
        else:
            logging.warning(f"No Specimens for combo: ({match.best_b1()}, {match.best_b2()}, {match.get_p1()}, {match.get_p2()})")
    return sample_id

def group_sample(match, sample_id):
    if sample_id in [SampleId.UNKNOWN, SampleId.AMBIGUOUS]:
        b1s = match.best_b1()
        b2s = match.best_b2()
        if len(b2s) == 1:
            sample_id = SampleId.PREFIX_REV_MATCH + b2s[0]
        elif len(b1s) == 1:
            sample_id = SampleId.PREFIX_FWD_MATCH + b1s[0]
    return sample_id

def match_sequence(parameters: MatchParameters, seq: str, rseq: str, specimens: Specimens) -> List[SequenceMatch]:
    """Match sequence against primers and barcodes"""

    matches = []

    for primer_pair in specimens.get_primer_pairs():
        match = SequenceMatch(len(seq), specimens.b_length())
        match_one_end(match, parameters, rseq, True, primer_pair, Primer.P1, Barcode.B1)
        match_one_end(match, parameters, seq, False, primer_pair, Primer.P2, Barcode.B2)

        matches.append(match)

    return matches


def match_one_end(match: SequenceMatch, parameters: MatchParameters, sequence: str,
                  reversed_sequence: bool, primer_pair: PrimerPairInfo,
                  which_primer: Primer, which_barcode: Barcode) -> None:
    """Match primers and barcodes at one end of the sequence."""

    primer = primer_pair.p1 if which_primer == Primer.P1 else primer_pair.p2
    primer_rc = primer_pair.p1_rc if which_primer == Primer.P1 else primer_pair.p2_rc

    primer_match = align_seq(primer_rc, sequence, parameters.max_dist_primer,
                             len(sequence) - parameters.search_len, len(sequence))

    # If we found matching primers, look for corresponding barcodes
    if primer_match.matched():
        match.set_primer_match(primer_match, primer, reversed_sequence, which_primer)

        # Get relevant barcodes for this primer pair
        barcodes = primer_pair.b1s if which_barcode == Barcode.B1 else primer_pair.b2s

        for b in barcodes:
            b_rc = reverse_complement(b)
            bd = None
            bm = None
            bc = None
            for l in primer_match.locations():
                barcode_match = align_seq(b_rc, sequence, parameters.max_dist_index,
                                          l[1] + 1, len(sequence), AlignMode.PREFIX)
                if barcode_match.matched():
                    if bd is None or barcode_match.distance() < bd:
                        bm = barcode_match
                        bd = barcode_match.distance()
                        bc = b
            if bm:
                match.add_barcode_match(bm, bc, reversed_sequence, which_barcode)


def determine_orientation(parameters: MatchParameters, seq: SeqRecord,
                        specimens: Specimens) -> Tuple[bool, SeqRecord, SeqRecord]:
    """Determine sequence orientation by checking primer pairs"""
    rseq = seq.reverse_complement()
    rseq.id = seq.id
    rseq.description = seq.description

    best_distance = float('inf')
    reverse_sequence = None

    for primer_pair in specimens.get_primer_pairs():
        p1 = primer_pair.p1
        p2 = primer_pair.p2

        # Check sequence in both orientations
        p1_match = align_seq(p1, seq, int(len(p1) * 0.5), 0, parameters.search_len)
        p1r_match = align_seq(p1, rseq, int(len(p1) * 0.5), 0, parameters.search_len)
        p2_match = align_seq(p2, seq, int(len(p2) * 0.5), 0, parameters.search_len)
        p2r_match = align_seq(p2, rseq, int(len(p2) * 0.5), 0, parameters.search_len)

        vals = []
        if p1_match.matched(): vals.append(p1_match.distance())
        if p2_match.matched(): vals.append(p2_match.distance())
        if p1r_match.matched(): vals.append(p1r_match.distance())
        if p2r_match.matched(): vals.append(p2r_match.distance())

        if len(vals) == 0:
            continue  # Try next primer pair

        current_best = min(vals)
        if current_best < best_distance:
            best_distance = current_best

            if (p1_match.distance() == current_best or p2r_match.distance() == current_best) and \
                    (not p2_match.matched() or p2_match.distance() > current_best) and \
                    (not p1r_match.matched() or p1r_match.distance() > current_best):
                reverse_sequence = False
            elif (p1r_match.distance() == current_best or p2_match.distance() == current_best) and \
                    (not p2r_match.matched() or p2r_match.distance() > current_best) and \
                    (not p1_match.matched() or p1_match.distance() > current_best):
                reverse_sequence = True
            else:
                logging.debug(
                    f"Could not determine orientation for primer pair: {p1_match.distance()} {p2_match.distance()} {p1r_match.distance()} {p2r_match.distance()} {len(seq)}")
                continue  # Try next primer pair

    if reverse_sequence is None:
        logging.debug(f"Could not orient sequence - no unambiguous primer matches")
        return False, seq, rseq

    if reverse_sequence:
        return True, rseq, seq
    return True, seq, rseq


def output_write_operation(write_op: WriteOperation,
                           output_manager: OutputManager,
                           args: argparse.Namespace) -> None:
    if not args.output_to_files:
        fh = sys.stdout
        formatted_seq = write_op.sequence
        if args.color:
            formatted_seq = color_sequence(formatted_seq, write_op.match,
                                           write_op.quality_scores)

        header_symbol = '@' if args.isfastq else '>'
        fh.write(f"{header_symbol}{write_op.seq_id} {write_op.distance_code} {write_op.sample_id}\n")
        fh.write(f"{formatted_seq}\n")
        if args.isfastq:
            fh.write("+\n")
            fh.write("".join(chr(q + 33) for q in write_op.quality_scores) + "\n")
    else:
        # Use new buffered file manager
        header = f"{write_op.seq_id} {write_op.distance_code} {write_op.sample_id}"
        output_manager.write_sequence(write_op.sample_id, header,
                                      write_op.sequence, write_op.quality_scores)

def color_sequence(seq, match, quality_scores):
    blue = "\033[0;34m"
    green = "\033[0;32m"
    red = "\033[0;31m"
    color_reset = "\033[0m"  # No Color (reset)

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
                        colored_seq[i] = color + seq[i].lower() + color_reset
                    else:
                        colored_seq[i] = color + seq[i] + color_reset

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


def version():
    # 0.1 September 14, 2024 - rewritten from minibar.py
    # 0.2 November 10, 2024 - support for multiple primer pairs
    # 0.3 December 4, 2024 - code & doc cleanup, write pooling
    return "specimux.py version 0.3"

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Specimux: Demultiplex MinION sequences by dual barcode indexes and primers.")
    
    parser.add_argument("barcode_file", help="File containing barcode information")
    parser.add_argument("sequence_file", help="Sequence file in Fasta or Fastq format, gzipped or plain text")

    parser.add_argument("--min-length", type=int, default=-1, help="Minimum sequence length.  Shorter sequences will be skipped (default: no filtering)")
    parser.add_argument("--max-length", type=int, default=-1, help="Maximum sequence length.  Longer sequences will be skipped (default: no filtering)")
    parser.add_argument("-n", "--num-seqs", type=str, default="-1", help="Number of sequences to read from file (e.g., -n 100 or -n 102,3)")
    parser.add_argument("-e", "--index-edit-distance", type=int, default=-1, help="Barcode edit distance value, default is half of min distance between barcodes")
    parser.add_argument("-E", "--primer-edit-distance", type=int, default=-1, help="Primer edit distance value, default is min distance between primers")
    parser.add_argument("-l", "--search-len", type=int, default=80, help="Length to search for index and primer at start and end of sequence (default: 80)")
    parser.add_argument("--group-unknowns", action="store_true", help="Group unknown sequences based on partial matches and classifications")
    parser.add_argument("-F", "--output-to-files", action="store_true", help="Create individual sample files for sequences")
    parser.add_argument("-P", "--output-file-prefix", default="sample_", help="Prefix for individual files when using -F (default: sample_)")
    parser.add_argument("-O", "--output-dir", default=".", help="Directory for individual files when using -F (default: .)")
    parser.add_argument("--color", action="store_true", help="Highlight barcode matches in blue, primer matches in green")
    parser.add_argument("--trim", choices=[TrimMode.NONE, TrimMode.TAILS, TrimMode.BARCODES, TrimMode.PRIMERS], default=TrimMode.BARCODES, help="trimming to apply")
    parser.add_argument("-d", "--diagnostics", action="store_true", help="Output extra diagnostics")
    parser.add_argument("-D", "--debug", action="store_true", help="Enable debug logging")
    parser.add_argument("--top-unmatched-barcodes", type=int, default=0, help="Display the top N unmatched barcode strings")

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



def worker(work_item: WorkItem, specimens: Specimens, args: argparse.Namespace) -> Tuple[List[WriteOperation], Counter, Counter]:
    try:
        l = logging.DEBUG if args.debug else logging.INFO
        logging.basicConfig(level=l, format='%(asctime)s - %(levelname)s - %(message)s')

        return process_sequences(work_item.seq_records, work_item.parameters, specimens, args)
    except Exception as e:
        logging.error(traceback.format_exc())
        raise WorkerException(e)

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
    unmatched_barcodes = Counter()

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
                    work_item = WorkItem(num_seqs, seq_batch, parameters)
                    yield work_item
                    num_seqs += len(seq_batch)

            try:
                for write_ops, batch_classifications, batch_unmatched in pool.imap_unordered(worker_func, producer()):
                    classifications += batch_classifications
                    unmatched_barcodes += batch_unmatched
                    for write_op in write_ops:
                        output_write_operation(write_op, output_manager, args)

                    total_processed = sum(classifications.values())
                    matched = classifications[MatchCode.MATCHED]
                    pct = matched/total_processed if total_processed > 0 else 0
                    print(f"read: {total_processed}\t\tmatched: {matched}\t\t{pct:.2%}", end='\r', file=sys.stderr)
                    sys.stderr.flush()
            except WorkerException as e:
                logging.error(f"Unexpected error in worker (see details above): {e}")
                sys.exit(1)

    elapsed = timeit.default_timer() - start_time

    print()
    logging.info(f"Elapsed time: : {elapsed:.2f} seconds")

    output_diagnostics(args, classifications, unmatched_barcodes)

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
    unmatched_barcodes = Counter()

    with OutputManager(args.output_dir, args.output_file_prefix, args.isfastq) as output_manager:
        num_seqs = 0
        while all_seqs or num_seqs < last_seq_to_output:
            to_read = sequence_block_size if all_seqs else min(sequence_block_size, last_seq_to_output - num_seqs)
            seq_batch = list(itertools.islice(seq_records, to_read))
            if not seq_batch:
                break
            write_ops, batch_classifications, batch_unmatched = process_sequences(seq_batch, parameters, specimens, args)
            classifications += batch_classifications
            unmatched_barcodes += batch_unmatched
            for write_op in write_ops:
                output_write_operation(write_op, output_manager, args)

            num_seqs += len(seq_batch)

            total_processed = sum(classifications.values())
            matched = classifications[MatchCode.MATCHED]
            pct = matched / total_processed if total_processed > 0 else 0
            print(f"read: {total_processed}\t\tmatched: {matched}\t\t{pct:.2%}", end='\r', file=sys.stderr)
            sys.stderr.flush()

    elapsed = timeit.default_timer() - start_time

    print()
    logging.info(f"Elapsed time: : {elapsed:.2f} seconds")

    output_diagnostics(args, classifications, unmatched_barcodes)

def output_diagnostics(args, classifications, unmatched_barcodes):
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
        logging.info(f"{len(unmatched_barcodes)} distinct barcode strings were unmatched")
        for barcode, count in unmatched_barcodes.most_common(args.top_unmatched_barcodes):
            logging.info(f"Barcode: {barcode}\t\t{count:>{max_count_length}}")

def setup_match_parameters(args, specimens):
    def _calculate_distances(sequences):
        distances = []
        for seq1, seq2 in itertools.combinations(sequences, 2):
            dist = edlib.align(seq1, seq2, task="distance")["editDistance"]
            distances.append(dist)
        return min(distances), Counter(distances)

    def _sanity_check_distance(sequences, desc, args):
        if len(sequences) <= 1:
            return
        m, c = _calculate_distances(sequences)
        if args.diagnostics:
            logging.info(f"Minimum edit distance is {m} for {desc}")
        return m

    # Collect all barcodes across primer pairs
    all_b1s = set()
    all_b2s = set()
    for pair in specimens.get_primer_pairs():
        all_b1s.update(pair.b1s)
        all_b2s.update(pair.b2s)

    _sanity_check_distance(list(all_b1s), "Forward Barcodes", args)
    _sanity_check_distance(list(all_b2s), "Reverse Barcodes", args)

    combined_barcodes = list(all_b1s) + [reverse_complement(b) for b in all_b2s]
    min_bc_dist = _sanity_check_distance(combined_barcodes,
        "Forward Barcodes + Reverse Complement of Reverse Barcodes", args)

    primers = []
    for ppi in specimens.get_primer_pairs():
        primers.append(ppi.p1)
        primers.append(ppi.p2)
        primers.append(ppi.p1_rc)
        primers.append(ppi.p2_rc)
    primers = list(set(primers))
    min_primer_dist = _sanity_check_distance(primers, "All Primers and Reverse Complements", args)

    combined_barcodes_and_primers = combined_barcodes + primers
    _sanity_check_distance(combined_barcodes_and_primers,
                           "Forward Barcodes + Reverse Complement of Reverse Barcodes + All Primers", args)

    max_search_area = args.search_len

    max_dist_index = math.ceil(min_bc_dist / 2.0)
    max_dist_primer = min_primer_dist
    if args.index_edit_distance != -1: max_dist_index = args.index_edit_distance
    if args.primer_edit_distance != -1: max_dist_primer = args.primer_edit_distance

    logging.info(f"Using Edit Distance Thresholds {max_dist_index} (barcode) and {max_dist_primer} (primer)")

    parameters = MatchParameters(max_dist_primer, max_dist_index, max_search_area)
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
