#!/usr/bin/env python3

"""
Sequence alignment and quality utilities.

This module contains alignment functionality
extracted from the original specimux.py monolithic file.
"""

import logging
from typing import List, Optional, Tuple, Union

import edlib
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord

from .constants import AlignMode, IUPAC_CODES, IUPAC_EQUIV
from .models import AlignmentResult


def align_seq(query: Union[str, Seq, SeqRecord],
             target: Union[str, Seq, SeqRecord],
             max_distance: int,
             start: int,
             end: int,
             mode: str = AlignMode.INFIX) -> AlignmentResult:
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

    m = AlignmentResult(r)
    m.adjust_start(s)
    return m

def get_quality_seq(seq):
    if "phred_quality" in seq.letter_annotations:
        return seq.letter_annotations["phred_quality"]
    else:
        return [40]*len(seq)


def color_sequence(seq: str, quality_scores: List[int], p1_location: Tuple[int, int],
                   p2_location: Tuple[int, int], b1_location: Tuple[int, int], b2_location: Tuple[int, int]):
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
            cstart, cend = location
            if cstart < 0 or cend < 0:
                return

            for i in range(cstart, cend + 1):  # Include the end position
                if i < seq_len:
                    if quality_scores[i] < 10:
                        colored_seq[i] = color + seq[i].lower() + color_reset
                    else:
                        colored_seq[i] = color + seq[i] + color_reset

    # Color barcode1 (blue)
    color_region(b1_location, blue)

    # Color primer1 (green)
    color_region(p1_location, green)

    # Color primer2 (green)
    color_region(p2_location, green)

    # Color barcode2 (blue)
    color_region(b2_location, blue)

    # Fill in uncolored regions
    for i in range(seq_len):
        if colored_seq[i] == '':
            if quality_scores[i] < 10:
                colored_seq[i] = seq[i].lower()
            else:
                colored_seq[i] = seq[i]

    return ''.join(colored_seq[start:end])

