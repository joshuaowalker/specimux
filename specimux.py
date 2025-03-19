#!/usr/bin/env python3

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
import fcntl
import hashlib
import itertools
import logging
import math
import multiprocessing
import os
import sys
import timeit
import traceback
from _operator import itemgetter
from collections import Counter
from collections import defaultdict
from contextlib import contextmanager
from enum import Enum
from functools import partial
from multiprocessing import Pool
from operator import itemgetter
from typing import List, Tuple, Dict
from typing import NamedTuple
from typing import Optional, Set
from typing import Union
from typing import Protocol
from tqdm import tqdm

import edlib
import pybloomfilter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import reverse_complement
from Bio.SeqRecord import SeqRecord
from cachetools import LRUCache
import gzip
import io
import mmap
import tempfile
import shutil
import hashlib
import glob

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
    NO_PRIMERS = "2C: No Primer Matches"
    NO_FWD_PRIMER = "2A: No Forward Primer Matches"
    NO_REV_PRIMER = "2B: No Reverse Primer Matches"
    NO_BARCODES = "3E: No Barcode Matches"
    NO_FWD_BARCODE = "3A: No Forward Barcode Matches"
    FWD_BARCODE_TRUNCATED = "3B: No Forward Barcode Matches (May be truncated)"
    NO_REV_BARCODE = "3C: No Reverse Barcode Matches"
    REV_BARCODE_TRUNCATED = "3D: No Reverse Barcode Matches (May be truncated)"
    MULTIPLE_PRIMER_PAIRS_MATCHED = "4: Multiple Primer/Orientation Full Matches"
    AMBIGUOUS_BARCODES = "5C: Multiple Matches for Both Barcodes"
    AMBIGUOUS_FWD_BARCODE = "5A: Multiple Matches for Forward Barcode"
    AMBIGUOUS_REV_BARCODE = "5B: Multiple Matches for Reverse Barcode"
    UNKNOWN_COMBINATION = "5D: No Specimen for Barcodes"
    MATCHED = " 6: Matched"

class Barcode(Enum):
    B1 = 1
    B2 = 2

class Primer(Enum):
    FWD = 3
    REV = 4


class PrimerInfo:
    def __init__(self, name: str, seq: str, direction: Primer, pools: List[str]):
        self.name = name
        self.primer = seq.upper()
        self.direction = direction
        self.primer_rc = reverse_complement(seq.upper())
        self.barcodes = set()
        self.specimens = set()
        self.pools = pools

class PrimerRegistry:
    """Manages primers and their relationships to pools"""

    def __init__(self):
        self._primers: Dict[str, PrimerInfo] = {}  # name -> PrimerInfo
        self._pools: Dict[str, Set[str]] = {}  # pool -> set of primer names
        self._pool_primers: Dict[str, Dict[Primer, List[PrimerInfo]]] = {}  # pool -> {direction -> [PrimerInfo]}

    def add_primer(self, primer: PrimerInfo, pools: List[str]) -> None:
        """
        Add a primer to the registry and its pools

        Args:
            primer: PrimerInfo object to add
            pools: List of pool names this primer belongs to

        Raises:
            ValueError: If primer name already exists
        """
        if primer.name in self._primers:
            raise ValueError(f"Duplicate primer name: {primer.name}")

        # Store primer
        self._primers[primer.name] = primer

        # Add to pools
        for pool in pools:
            if pool not in self._pools:
                self._pools[pool] = set()
                self._pool_primers[pool] = {
                    Primer.FWD: [],
                    Primer.REV: []
                }
            self._pools[pool].add(primer.name)
            self._pool_primers[pool][primer.direction].append(primer)

    def get_primer(self, name: str) -> Optional[PrimerInfo]:
        """Get a primer by name"""
        return self._primers.get(name)

    def get_primers_in_pool(self, pool: str) -> List[PrimerInfo]:
        """Get all primers in a pool"""
        if pool not in self._pools:
            return []
        primers = []
        for direction in [Primer.FWD, Primer.REV]:
            primers.extend(self._pool_primers[pool][direction])
        return primers

    def get_pools(self) -> List[str]:
        """Get list of all pool names"""
        return list(self._pools.keys())

    def get_pool_primers(self, pool: str, direction: Optional[Primer] = None) -> List[PrimerInfo]:
        """
        Get primers in a pool, optionally filtered by direction

        Args:
            pool: Name of the pool
            direction: Optional primer direction to filter by

        Returns:
            List of PrimerInfo objects
        """
        if pool not in self._pools:
            return []
        if direction:
            return self._pool_primers[pool][direction]
        return self.get_primers_in_pool(pool)

    def primer_in_pool(self, primer_name: str, pool: str) -> bool:
        """Check if a primer is in a pool"""
        return pool in self._pools and primer_name in self._pools[pool]

    def validate_pools(self) -> None:
        """
        Validate pool configurations

        Raises:
            ValueError: If validation fails
        """
        for pool in self._pools:
            # Each pool must have at least one forward and one reverse primer
            if not self._pool_primers[pool][Primer.FWD]:
                raise ValueError(f"Pool {pool} has no forward primers")
            if not self._pool_primers[pool][Primer.REV]:
                raise ValueError(f"Pool {pool} has no reverse primers")

    def get_pool_stats(self) -> Dict:
        """Get statistics about pools and primers"""
        stats = {
            'total_primers': len(self._primers),
            'total_pools': len(self._pools),
            'pools': {}
        }

        for pool in self._pools:
            stats['pools'][pool] = {
                'forward_primers': len(self._pool_primers[pool][Primer.FWD]),
                'reverse_primers': len(self._pool_primers[pool][Primer.REV]),
                'total_primers': len(self.get_primers_in_pool(pool))
            }

        return stats


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
    def __init__(self, sequence, barcode_length):
        self.sequence_length = len(sequence)
        self.sequence = sequence
        self.p1_match: Optional[MatchResult] = None
        self._b1_matches: List[Tuple[str, MatchResult, float]] = []
        self.p2_match: Optional[MatchResult] = None
        self._b2_matches: List[Tuple[str, MatchResult, float]] = []
        self._p1: Optional[PrimerInfo] = None
        self._p2: Optional[PrimerInfo] = None
        self._pool: Optional[str] = None  # New: track the pool
        self.ambiguity_threshold = 1.0
        self._barcode_length = barcode_length

    def set_pool(self, pool: str):
        """Set the primer pool for this match"""
        self._pool = pool

    def get_pool(self) -> Optional[str]:
        """Get the primer pool for this match"""
        return self._pool

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

    def set_primer_match(self, match: MatchResult, primer: PrimerInfo, reverse: bool, which: Primer):
        m = match
        if reverse:
            m = m.reversed(self.sequence_length)
        if which is Primer.FWD:
            self.p1_match = m
            self._p1 = primer
        else:
            self.p2_match = m
            self._p2 = primer

    def get_p1(self) -> PrimerInfo:
        return self._p1

    def get_p2(self) -> PrimerInfo:
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


class Specimens:
    def __init__(self, primer_registry: PrimerRegistry):
        self._specimens = []  # List of (id, pool, b1, p1s, b2, p2s) tuples for reference
        self._barcode_length = 0
        self._primers = {}
        self._specimen_ids = set()
        self._primer_pairings = {}
        self._primer_registry = primer_registry
        self._active_pools = set()  # Track pools actually used by specimens

    def add_specimen(self, specimen_id: str, pool: str, b1: str, p1: str, b2: str, p2: str):
        """Add a specimen with its barcodes and primers"""
        if specimen_id in self._specimen_ids:
            raise ValueError(format(f"Duplicate specimen id in index file: {specimen_id}"))
        self._specimen_ids.add(specimen_id)

        # Track active pools
        self._active_pools.add(pool)

        self._barcode_length = max(self._barcode_length, len(b1), len(b2))

        # Handle wildcards and get list of possible primers
        p1_list = self._resolve_primer_name(p1, pool, Primer.FWD)
        p2_list = self._resolve_primer_name(p2, pool, Primer.REV)

        # Register primers and barcodes
        for p1_info in p1_list:
            ps1 = p1_info.primer
            if ps1 not in self._primers:
                self._primers[ps1] = p1_info
            primer_info = self._primers[ps1]
            primer_info.barcodes.add(b1)
            primer_info.specimens.add(specimen_id)

        for p2_info in p2_list:
            ps2 = p2_info.primer
            if ps2 not in self._primers:
                self._primers[ps2] = p2_info
            primer_info = self._primers[ps2]
            primer_info.barcodes.add(b2)
            primer_info.specimens.add(specimen_id)

        self._specimens.append((specimen_id, pool, b1, p1_list, b2, p2_list))

    def prune_unused_pools(self):
        """Remove pools that aren't used by any specimens"""
        # Get list of all unused pools
        all_pools = set(self._primer_registry.get_pools())
        unused_pools = all_pools - self._active_pools

        if unused_pools:
            logging.info(f"Removing unused pools: {unused_pools}")

            # Remove unused pools from all primers
            for primer in self._primers.values():
                primer.pools = [p for p in primer.pools if p in self._active_pools]

            # Update primer registry
            for pool in unused_pools:
                # Remove pool from registry's internal data structures
                if pool in self._primer_registry._pools:
                    del self._primer_registry._pools[pool]
                if pool in self._primer_registry._pool_primers:
                    del self._primer_registry._pool_primers[pool]

            # Update pool stats in logs
            stats = self._primer_registry.get_pool_stats()
            logging.info(f"After pruning: {stats['total_primers']} primers in {stats['total_pools']} pools")
            for pool, pool_stats in stats['pools'].items():
                logging.info(f"Pool {pool}: {pool_stats['forward_primers']} forward, "
                             f"{pool_stats['reverse_primers']} reverse primers")

    def _resolve_primer_name(self, primer_name: str, pool: str, direction: Primer) -> List[PrimerInfo]:
        """Resolve a primer name (including wildcards) to a list of PrimerInfo objects"""
        if primer_name == '-' or primer_name == '*':  # Handle wildcards
            # Get all primers in the specified pool and direction
            primers = []
            for p in self._primer_registry.get_primers_in_pool(pool):
                if p.direction == direction:
                    primers.append(p)
            if not primers:
                raise ValueError(f"No {direction.name} primers found in pool {pool}")
            return primers
        else:
            # Get specific primer
            primer = self._primer_registry.get_primer(primer_name)
            if not primer:
                raise ValueError(f"Primer not found: {primer_name}")
            if primer.direction != direction:
                raise ValueError(f"Primer {primer_name} is not a {direction.name} primer")
            if not self._primer_registry.primer_in_pool(primer_name, pool):
                raise ValueError(f"Primer {primer_name} is not in pool {pool}")
            return [primer]

    def specimens_for_barcodes_and_primers(self, b1_list: List[str], b2_list: List[str],
                                           p1_matched: PrimerInfo, p2_matched: PrimerInfo) -> List[str]:
        matching_specimens = []
        for spec_id, pool, b1, p1s, b2, p2s in self._specimens:
            if (p1_matched in p1s and
                    p2_matched in p2s and
                    b1.upper() in b1_list and
                    b2.upper() in b2_list):
                matching_specimens.append(spec_id)

        return matching_specimens

    def get_primers(self, direction: Primer) -> List[PrimerInfo]:
        return [p for p in self._primers.values() if p.direction == direction]

    def get_paired_primers(self, primer: str) -> List[PrimerInfo]:
        if primer in self._primer_pairings:
            return self._primer_pairings[primer]

        specimens = self._primers[primer].specimens
        direction = self._primers[primer].direction
        rv = []
        for pi in self._primers.values():
            if direction != pi.direction and pi.specimens.intersection(specimens):
                rv.append(pi)

        self._primer_pairings[primer] = rv
        return rv

    def get_specimen_pool(self, specimen_id: str) -> Optional[str]:
        """Get the primer pool for a specimen"""
        for spec_id, pool, _, _, _, _ in self._specimens:
            if spec_id == specimen_id:
                return pool
        return None

    def b_length(self):
        return self._barcode_length

    def validate(self):
        self._validate_barcodes_globally_unique()
        self._validate_barcode_lengths()
        self.prune_unused_pools()

    def _validate_barcodes_globally_unique(self):
        all_b1s = set()
        all_b2s = set()
        for primer in self._primers.values():
            if primer.direction == Primer.FWD:
                all_b1s.update(primer.barcodes)
            else:
                all_b2s.update(primer.barcodes)
        dups = all_b1s.intersection(all_b2s)
        if len(dups) > 0:
            logging.warning(f"Duplicate Barcodes ({len(dups)}) in Fwd and Rev: {dups}")

    def _validate_barcode_lengths(self):
        all_b1s = set()
        all_b2s = set()
        for primer in self._primers.values():
            if primer.direction == Primer.FWD:
                all_b1s.update(primer.barcodes)
            else:
                all_b2s.update(primer.barcodes)
        if len(set(len(b) for b in all_b1s)) > 1:
            logging.warning("Forward barcodes have inconsistent lengths")
        if len(set(len(b) for b in all_b2s)) > 1:
            logging.warning("Reverse barcodes have inconsistent lengths")

class MatchParameters:
    def __init__(self, max_dist_primers: Dict[str, int], max_dist_index: int, search_len: int, preorient: bool):
        self.max_dist_primers = max_dist_primers
        self.max_dist_index = max_dist_index
        self.search_len = search_len
        self.preorient = preorient

class WriteOperation(NamedTuple):
    sample_id: str
    seq_id: str
    distance_code: str
    sequence: str
    quality_sequence: str
    quality_scores: List[int]
    p1_location: Tuple[int, int]
    p2_location: Tuple[int, int]
    b1_location: Tuple[int, int]
    b2_location: Tuple[int, int]
    primer_pool: str
    p1_name: str
    p2_name: str

class WorkItem(NamedTuple):
    seq_number: int
    seq_records: List  # This now contains actual sequence records, not an iterator
    parameters: MatchParameters


class FileHandleCache(LRUCache):
    """Custom LRU cache that works with CachedFileManager to close file handles on eviction."""

    def __init__(self, maxsize, file_manager):
        super().__init__(maxsize)
        self.file_manager = file_manager

    def popitem(self):
        """Override popitem to ensure file handle and lock cleanup on eviction."""
        key, file_handle = super().popitem()
        # Ensure buffer is flushed before closing
        self.file_manager.flush_buffer(key)
        file_handle.close()
        # Also close the corresponding lock if it exists
        self.file_manager.close_lock(key)
        return key, file_handle

class CachedFileManager:
    """Manages a pool of file handles with LRU caching, write buffering, and locking."""

    def __init__(self, max_open_files: int, buffer_size: int, output_dir: str):
        self.max_open_files = max_open_files
        self.buffer_size = buffer_size
        self.file_cache = FileHandleCache(maxsize=max_open_files, file_manager=self)
        self.write_buffers = defaultdict(list)

        # Create lock directory
        self.lock_dir = os.path.join(output_dir, ".specimux_locks")
        os.makedirs(self.lock_dir, exist_ok=True)
        self.locks = {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Ensure all buffers are flushed and files are closed on exit."""
        try:
            self.flush_all()
        finally:
            self.close_all()

    def write(self, filename: str, data: str):
        """Write data to a file through the buffer system."""
        self.write_buffers[filename].append(data)

        if len(self.write_buffers[filename]) >= self.buffer_size:
            self.flush_buffer(filename)

    def flush_buffer(self, filename: str):
        """Flush the buffer for a specific file to disk."""
        if not self.write_buffers[filename]:
            return

        # Get or create lock
        if filename not in self.locks:
            lock_path = self._get_lock_path(filename)
            fd = os.open(lock_path, os.O_CREAT | os.O_RDWR, 0o666)
            self.locks[filename] = fd

        # Hold lock for entire write operation
        fcntl.lockf(self.locks[filename], fcntl.LOCK_EX)
        try:
            buffer_data = ''.join(self.write_buffers[filename])

            try:
                f = self.file_cache[filename]
            except KeyError:
                f = open(filename, 'a')  # Always append mode
                self.file_cache[filename] = f

            f.write(buffer_data)
            f.flush()
            os.fsync(f.fileno())  # Ensure data is written to disk
            self.write_buffers[filename].clear()

        finally:
            fcntl.lockf(self.locks[filename], fcntl.LOCK_UN)

    def flush_all(self):
        """Flush all buffers to disk."""
        for filename in list(self.write_buffers.keys()):
            self.flush_buffer(filename)

    def close_all(self):
        """Close all files and release all locks."""
        try:
            self.flush_all()
            for f in self.file_cache.values():
                f.close()
            self.file_cache.clear()
        finally:
            # Always clean up locks
            for fd in self.locks.values():
                try:
                    os.close(fd)
                except Exception as e:
                    logging.warning(f"Error closing lock file: {e}")
            self.locks.clear()

    def __del__(self):
        """Ensure cleanup on garbage collection"""
        self.close_all()

    def close_lock(self, filename: str):
        """Close a specific lock file if it's open."""
        if filename in self.locks:
            try:
                os.close(self.locks[filename])
                del self.locks[filename]
            except Exception as e:
                logging.warning(f"Error closing lock file for {filename}: {e}")

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

    def _get_lock_path(self, filename: str) -> str:
        """Get path for lock file in shared directory"""
        # Use hash of absolute path to avoid issues with long filenames
        abs_path = os.path.abspath(filename)
        file_hash = hashlib.md5(abs_path.encode()).hexdigest()
        return os.path.join(self.lock_dir, f"{file_hash}.lock")


class OutputManager:
    """Manages output files with pool-based organization."""

    def __init__(self, output_dir: str, prefix: str, is_fastq: bool,
                 max_open_files: int = 200, buffer_size: int = 500):
        self.output_dir = output_dir
        self.prefix = prefix
        self.is_fastq = is_fastq
        self.file_manager = CachedFileManager(max_open_files, buffer_size, output_dir)

    def __enter__(self):
        os.makedirs(self.output_dir, exist_ok=True)
        self.file_manager.__enter__()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return self.file_manager.__exit__(exc_type, exc_val, exc_tb)

    def _make_filename(self, sample_id: str, pool: str, p1: str, p2: str) -> str:
        """Create a filename with pool-based organization."""
        extension = '.fastq' if self.is_fastq else '.fasta'

        # Handle None values
        sample_id = sample_id or SampleId.UNKNOWN
        pool = pool or "unknown"
        p1 = p1 or "unknown"
        p2 = p2 or "unknown"

        safe_id = "".join(c if c.isalnum() or c in "._-$#" else "_" for c in sample_id)
        primer_dir = f"{p1}-{p2}"

        # Determine if this is a full match or a partial/unknown match
        if sample_id == SampleId.UNKNOWN:
            # Complete unknown
            return os.path.join(self.output_dir, pool, primer_dir, "unknown", f"{self.prefix}{safe_id}{extension}")
        elif sample_id == SampleId.AMBIGUOUS:
            # Ambiguous match
            return os.path.join(self.output_dir, pool, primer_dir, "ambiguous", f"{self.prefix}{safe_id}{extension}")
        elif sample_id.startswith(SampleId.PREFIX_FWD_MATCH):
            # Forward barcode matched but reverse didn't
            return os.path.join(self.output_dir, pool, primer_dir, "partial", f"{self.prefix}{safe_id}{extension}")
        elif sample_id.startswith(SampleId.PREFIX_REV_MATCH):
            # Reverse barcode matched but forward didn't
            return os.path.join(self.output_dir, pool, primer_dir, "partial", f"{self.prefix}{safe_id}{extension}")
        else:
            # Full match
            return os.path.join(self.output_dir, pool, primer_dir, "full", f"{self.prefix}{safe_id}{extension}")

    def write_sequence(self, write_op: WriteOperation):
        """Write a sequence to the appropriate output file."""
        filename = self._make_filename(write_op.sample_id, write_op.primer_pool,
                                       write_op.p1_name, write_op.p2_name)

        # Ensure directory exists
        os.makedirs(os.path.dirname(filename), exist_ok=True)

        # Format header to include primer information
        header = (f"{write_op.seq_id} {write_op.distance_code} "
                  f"pool={write_op.primer_pool} "
                  f"primers={write_op.p1_name}+{write_op.p2_name} "
                  f"{write_op.sample_id}")

        output = []
        output.append('@' if self.is_fastq else '>')
        output.append(header + "\n")
        output.append(write_op.sequence + "\n")
        if self.is_fastq:
            output.append("+\n")
            output.append(write_op.quality_sequence + "\n")

        self.file_manager.write(filename, ''.join(output))

class WorkerException(Exception):
    pass

class BarcodePrefilter(Protocol):
    """Protocol defining the interface for barcode prefilters"""
    def match(self, barcode: str, sequence: str) -> bool:
        """Check if barcode potentially matches sequence"""
        ...

class PassthroughPrefilter:

    def match(self, barcode: str, sequence: str) -> bool:
        #always call edlib.align
        return True

class BloomPrefilter:
    """Fast barcode pre-filter using Bloom filters with proper file handling"""

    def __init__(self, barcodes: List[str], max_distance: int, error_rate: float = 0.05,
                 filename: Optional[str] = None):
        """Initialize bloom filter for barcode matching

        Args:
            barcodes: List of barcodes to add to filter
            max_distance: Maximum edit distance for matches
            error_rate: Acceptable false positive rate
            filename: Optional file to store the filter
        """
        if not barcodes:
            raise ValueError("Must provide at least one barcode")

        self.barcodes = list(set(barcodes))  # Deduplicate
        self.barcode_length = len(self.barcodes[0])
        self.max_distance = max_distance
        self.error_rate = error_rate
        self.min_length = self.barcode_length - max_distance

        # Calculate filter size using first barcode as sample
        sample_variants = len(self._generate_variants(self.barcodes[0]))
        total_combinations = sample_variants * len(self.barcodes)
        logging.info(f"Building Bloom filter for estimated {total_combinations} variants")

        # Create and populate filter
        if filename:
            self.bloom_filter = pybloomfilter.BloomFilter(total_combinations, error_rate, filename)
        else:
            self.bloom_filter = pybloomfilter.BloomFilter(total_combinations, error_rate)

        total_variants = 0
        # Add progress bar wrapping barcodes iteration
        for barcode in tqdm(self.barcodes, desc="Building Bloom filter", unit="barcode"):
            variants = self._generate_variants(barcode)
            total_variants += len(variants)
            for variant in variants:
                # Use fixed width format with truncation
                truncated = variant[:self.min_length]
                key = barcode + truncated
                self.bloom_filter.add(key)

        logging.info(f"Added {total_variants} actual variants to Bloom filter")

    def _generate_variants(self, barcode: str) -> Set[str]:
        """Generate all possible variants of a barcode within max_distance."""
        variants = {barcode}
        bases = "ACGT"

        def add_variants_recursive(current: str, distance_left: int, variants: Set[str]):
            if distance_left == 0:
                return

            # Substitutions
            for i in range(len(current)):
                for base in bases:
                    if base != current[i]:
                        new_str = current[:i] + base + current[i + 1:]
                        variants.add(new_str)
                        add_variants_recursive(new_str, distance_left - 1, variants)

            # Deletions
            for i in range(len(current)):
                new_str = current[:i] + current[i + 1:]
                variants.add(new_str)
                add_variants_recursive(new_str, distance_left - 1, variants)

            # Insertions
            for i in range(len(current) + 1):
                for base in bases:
                    new_str = current[:i] + base + current[i:]
                    variants.add(new_str)
                    add_variants_recursive(new_str, distance_left - 1, variants)

        add_variants_recursive(barcode, self.max_distance, variants)
        return variants

    @staticmethod
    def get_cache_path(barcodes: List[str], max_distance: int, error_rate: float = 0.05) -> str:

        """Generate a unique cache filename based on inputs."""
        m = hashlib.sha256()
        for bc in sorted(barcodes):  # Sort for consistency
            m.update(bc.encode())
        m.update(str(max_distance).encode())
        m.update(str(error_rate).encode())

        hash_prefix = m.hexdigest()[:16]
        cache_dir = os.path.join(os.path.expanduser("~"), ".specimux", "cache")
        os.makedirs(cache_dir, exist_ok=True)

        return os.path.join(cache_dir,
                            f"bloom_prefilter_{hash_prefix}_k{max_distance}_e{error_rate}.bf")

    @classmethod
    def create_filter(cls, barcode_rcs: List[str], max_distance: int, error_rate: float = 0.05) -> str:
        """Initialize bloom filter for barcode matching

        Args:
            barcode_rcs: Barcodes in the direction which they will be matched
            max_distance: Maximum edit distance for matches
            error_rate: Acceptable false positive rate

        Returns:
            str: Path to the created/cached bloom filter file
        """

        # Get cache paths
        cache_path = cls.get_cache_path(barcode_rcs, max_distance, error_rate)
        lock_path = cache_path + '.lock'

        # Create lock file if needed
        if not os.path.exists(lock_path):
            open(lock_path, 'w').close()

        # Create filter with proper locking
        with open(lock_path, 'r') as lock_file:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
            try:
                if not os.path.exists(cache_path):
                    logging.info("Creating initial Bloom filter...")
                    prefilter = cls(barcode_rcs, max_distance, error_rate, filename=cache_path)
                    prefilter.save(cache_path)
                    prefilter.close()
                    logging.info(f"Saved Bloom filter to cache: {cache_path}")
            finally:
                fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)

        return cache_path

    @classmethod
    def load_readonly(cls, filename: str, barcodes: List[str], max_distance: int) -> 'BloomPrefilter':
        """Load BloomPrefilter from file in read-only mode"""
        instance = cls.__new__(cls)
        instance.barcodes = list(set(barcodes))
        instance.barcode_length = len(instance.barcodes[0])
        instance.max_distance = max_distance
        instance.min_length = instance.barcode_length - max_distance

        try:
            # Open with read-only mmap
            instance.bloom_filter = pybloomfilter.BloomFilter.open(filename, mode='r')
        except Exception as e:
            raise ValueError(f"Failed to load read-only filter: {e}")

        return instance

    def save(self, filename: str):
        """Save filter to file with proper syncing"""
        self.bloom_filter.sync()

    def match(self, barcode: str, sequence: str) -> bool:
        """Check if barcode matches sequence within max_distance."""
        # Truncate sequence to minimum length
        truncated = sequence[:self.min_length]

        if barcode not in self.barcodes:
            return True

        # Concatenate in same order as when building filter
        key = barcode + truncated
        return key in self.bloom_filter

    def close(self):
        """Close the filter and release resources"""
        if hasattr(self, 'bloom_filter'):
            self.bloom_filter.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
def read_primers_file(filename: str) -> PrimerRegistry:
    """
    Read primers file and build primer registry

    Returns:
        PrimerRegistry object managing all primers and their relationships
    """
    registry = PrimerRegistry()

    for record in SeqIO.parse(filename, "fasta"):
        name = record.id
        sequence = str(record.seq)

        # Parse description for pool and position
        desc = record.description
        pool_names = []
        position = None

        for field in desc.split():
            if field.startswith("pool="):
                # Split pool names on either comma or semicolon
                pool_str = field[5:]
                pool_names = [p.strip() for p in pool_str.replace(';', ',').split(',')]
            elif field.startswith("position="):
                position = field[9:]

        if not pool_names:
            raise ValueError(f"Missing pool specification for primer {name}")
        if not position:
            raise ValueError(f"Missing position specification for primer {name}")

        if position == "forward":
            direction = Primer.FWD
        elif position == "reverse":
            direction = Primer.REV
        else:
            raise ValueError(f"Invalid primer position '{position}' for {name}")

        # Create PrimerInfo and add to registry
        primer = PrimerInfo(name, sequence, direction, pool_names)
        registry.add_primer(primer, pool_names)

    # Validate pool configurations
    registry.validate_pools()

    # Log pool statistics
    stats = registry.get_pool_stats()
    logging.info(f"Loaded {stats['total_primers']} primers in {stats['total_pools']} pools")
    for pool, pool_stats in stats['pools'].items():
        logging.info(f"Pool {pool}: {pool_stats['forward_primers']} forward, "
                    f"{pool_stats['reverse_primers']} reverse primers")

    return registry

def read_specimen_file(filename: str, primer_registry: PrimerRegistry) -> Specimens:
    """
    Read a tab-separated specimen file and return a Specimens object.
    Expected columns: SampleID, PrimerPool, FwIndex, FwPrimer, RvIndex, RvPrimer
    """
    specimens = Specimens(primer_registry)

    expected_columns = {'SampleID', 'PrimerPool', 'FwIndex', 'FwPrimer', 'RvIndex', 'RvPrimer'}

    with open(filename, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')

        # Validate columns
        missing_cols = expected_columns - set(reader.fieldnames)
        if missing_cols:
            raise ValueError(f"Missing required columns in specimen file: {missing_cols}")

        for row_num, row in enumerate(reader, start=1):
            try:
                specimens.add_specimen(
                    specimen_id=row['SampleID'],
                    pool=row['PrimerPool'],
                    b1=row['FwIndex'].upper(),
                    p1=row['FwPrimer'],
                    b2=row['RvIndex'].upper(),
                    p2=row['RvPrimer']
                )
            except (KeyError, ValueError) as e:
                raise ValueError(f"Error processing row {row_num}: {e}")

    if len(specimens._specimens) == 0:
        raise ValueError("No valid data found in the specimen file")

    return specimens

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
        p2_name=p2_name
    )


def process_sequences(seq_records: List[SeqRecord],
                      parameters: MatchParameters,
                      specimens: Specimens,
                      args: argparse.Namespace,
                      prefilter: BarcodePrefilter) -> Tuple[List[WriteOperation], Counter]:
    classifications = Counter()
    write_ops = []

    for seq in seq_records:
        sample_id = SampleId.UNKNOWN
        pool = None

        match = SequenceMatch(seq, specimens.b_length())

        if args.min_length != -1 and len(seq) < args.min_length:
            classification = MatchCode.SHORT_SEQ
        elif args.max_length != -1 and len(seq) > args.max_length:
            classification = MatchCode.LONG_SEQ
        else:
            rseq = seq.reverse_complement()
            rseq.id = seq.id
            rseq.description = seq.description

            matches = match_sequence(prefilter, parameters, seq, rseq, specimens)
            match, multimatch = choose_best_match(matches)
            if multimatch:
                classification = MatchCode.MULTIPLE_PRIMER_PAIRS_MATCHED
                sample_id = SampleId.AMBIGUOUS
            else:
                sample_id, pool = match_sample(match, sample_id, specimens)
                classification = classify_match(match, sample_id, specimens)
                if pool:
                    match.set_pool(pool)  # Set the pool in the match object for output

            sample_id = group_sample(match, sample_id)

        classifications[classification] += 1
        if args.debug:
            logging.debug(f"{classification}")

        write_op = create_write_operation(sample_id, args, match.sequence, match)
        write_ops.append(write_op)

    return write_ops, classifications


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
            classification = MatchCode.UNKNOWN_COMBINATION
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
    best_pool = None

    if not matches:
        raise ValueError("No matches provided to choose_best_match")

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
        if score > best_score:
            best_score = score
            best_match = m
            best_pool = m.get_pool()
            multimatch = False
        elif score == best_score:
            multimatch = True
            if best_pool != m.get_pool():
                best_pool = "unknown"

    if best_match is None:
        raise RuntimeError("Failed to select best match despite having matches")

    best_match.set_pool(best_pool)

    # ambiguity between matches only really matters when there was a full match
    if best_score == 5 and multimatch:
        return best_match, True
    else:
        return best_match, False



def match_sample(match: SequenceMatch, sample_id: str, specimens: Specimens) -> Tuple[str, Optional[str]]:
    """Returns tuple of (sample_id, pool_name)"""
    # Start with previously determined pool from primers
    pool = match.get_pool()

    if match.p1_match and match.p2_match and match.has_b1_match() and match.has_b2_match():
        ids = specimens.specimens_for_barcodes_and_primers(
            match.best_b1(), match.best_b2(), match.get_p1(), match.get_p2())
        if len(ids) > 1:
            sample_id = SampleId.AMBIGUOUS
            # For ambiguous matches, use most common pool from matching specimens
            pools = Counter()
            for id in ids:
                pools[specimens.get_specimen_pool(id)] += 1
            pool = pools.most_common(1)[0][0]
        elif len(ids) == 1:
            sample_id = ids[0]
            # For unique matches, use specimen's pool
            pool = specimens.get_specimen_pool(ids[0])
        else:
            logging.warning(
                f"No Specimens for combo: ({match.best_b1()}, {match.best_b2()}, "
                f"{match.get_p1()}, {match.get_p2()})")

    return sample_id, pool

def group_sample(match, sample_id):
    if sample_id in [SampleId.UNKNOWN, SampleId.AMBIGUOUS]:
        b1s = match.best_b1()
        b2s = match.best_b2()
        if len(b2s) == 1:
            sample_id = SampleId.PREFIX_REV_MATCH + b2s[0]
        elif len(b1s) == 1:
            sample_id = SampleId.PREFIX_FWD_MATCH + b1s[0]
    return sample_id

class Orientation(Enum):
    FORWARD = 1
    REVERSE = 2
    UNKNOWN = 3


def determine_orientation(parameters: MatchParameters, seq: str, rseq: str,
                          fwd_primers: List[PrimerInfo], rev_primers: List[PrimerInfo]) -> Orientation:
    """Determine sequence orientation by checking primer matches.
    Returns FORWARD/REVERSE only when highly confident, otherwise UNKNOWN."""

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
        fwd_p2 = align_seq(primer.primer, seq, parameters.max_dist_primers[primer.primer],
                           len(seq) - parameters.search_len, len(seq))
        rev_p2 = align_seq(primer.primer, rseq, parameters.max_dist_primers[primer.primer],
                           len(seq) - parameters.search_len, len(seq))
        if fwd_p2.matched():
            forward_matches += 1
        if rev_p2.matched():
            reverse_matches += 1

    # Only return an orientation if we have clear evidence in one direction
    # and no evidence in the other
    if forward_matches > 0 and reverse_matches == 0:
        return Orientation.FORWARD
    elif reverse_matches > 0 and forward_matches == 0:
        return Orientation.REVERSE
    else:
        return Orientation.UNKNOWN

def get_pool_from_primers(p1: Optional[PrimerInfo], p2: Optional[PrimerInfo]) -> Optional[str]:
    """
    Determine the primer pool based on matched primers.

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


def match_sequence(prefilter: BarcodePrefilter, parameters: MatchParameters, seq: SeqRecord,
                   rseq: SeqRecord, specimens: Specimens) -> List[SequenceMatch]:
    """Match sequence against primers and barcodes"""
    # extract string versions for performance - roughly 18% improvement
    s = str(seq.seq)
    rs = str(rseq.seq)

    if parameters.preorient:
        orientation = determine_orientation(parameters, s, rs,
                                            specimens.get_primers(Primer.FWD),
                                            specimens.get_primers(Primer.REV))
    else:
        orientation = Orientation.UNKNOWN

    matches = []

    for fwd_primer in specimens.get_primers(Primer.FWD):
        for rev_primer in specimens.get_paired_primers(fwd_primer.primer):
            if orientation in [Orientation.FORWARD, Orientation.UNKNOWN]:
                match = SequenceMatch(seq, specimens.b_length())
                match_one_end(prefilter, match, parameters, rs, True, fwd_primer,
                              Primer.FWD, Barcode.B1)
                match_one_end(prefilter, match, parameters, s, False, rev_primer,
                              Primer.REV, Barcode.B2)
                # Set initial pool based on primers
                pool = get_pool_from_primers(match.get_p1(), match.get_p2())
                if pool:
                    match.set_pool(pool)
                matches.append(match)

            if orientation in [Orientation.REVERSE, Orientation.UNKNOWN]:
                match = SequenceMatch(rseq, specimens.b_length())
                match_one_end(prefilter, match, parameters, s, True, fwd_primer,
                              Primer.FWD, Barcode.B1)
                match_one_end(prefilter, match, parameters, rs, False, rev_primer,
                              Primer.REV, Barcode.B2)
                # Set initial pool based on primers
                pool = get_pool_from_primers(match.get_p1(), match.get_p2())
                if pool:
                    match.set_pool(pool)
                matches.append(match)
    return matches

def match_one_end(prefilter: BarcodePrefilter, match: SequenceMatch, parameters: MatchParameters, sequence: str,
                  reversed_sequence: bool, primer_info: PrimerInfo,
                  which_primer: Primer, which_barcode: Barcode) -> None:
    """Match primers and barcodes at one end of the sequence."""

    primer = primer_info.primer
    primer_rc = primer_info.primer_rc

    primer_match = align_seq(primer_rc, sequence, parameters.max_dist_primers[primer],
                             len(sequence) - parameters.search_len, len(sequence))

    # If we found matching primers, look for corresponding barcodes
    if primer_match.matched():
        match.set_primer_match(primer_match, primer_info, reversed_sequence, which_primer)

        # Get relevant barcodes for this primer pair
        barcodes = primer_info.barcodes

        for b in barcodes:
            b_rc = reverse_complement(b)
            bd = None
            bm = None
            bc = None
            for l in primer_match.locations():
                target_seq = sequence[l[1] + 1:]
                if not prefilter.match(b_rc, target_seq):
                    continue

                barcode_match = align_seq(b_rc, sequence, parameters.max_dist_index,
                                          l[1] + 1, len(sequence), AlignMode.PREFIX)
                if barcode_match.matched():
                    if bd is None or barcode_match.distance() < bd:
                        bm = barcode_match
                        bd = barcode_match.distance()
                        bc = b
            if bm:
                match.add_barcode_match(bm, bc, reversed_sequence, which_barcode)

def output_write_operation(write_op: WriteOperation,
                           output_manager: OutputManager,
                           args: argparse.Namespace) -> None:
    if not args.output_to_files:
        fh = sys.stdout
        formatted_seq = write_op.sequence
        if args.color:
            formatted_seq = color_sequence(formatted_seq, write_op.quality_scores, write_op.p1_location, write_op.p2_location,
                                           write_op.b1_location, write_op.b2_location)

        header_symbol = '@' if args.isfastq else '>'
        fh.write(f"{header_symbol}{write_op.seq_id} {write_op.distance_code} {write_op.sample_id}\n")
        fh.write(f"{formatted_seq}\n")
        if args.isfastq:
            fh.write("+\n")
            fh.write(write_op.quality_sequence + "\n")
    else:
        output_manager.write_sequence(write_op)

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


def version():
    # 0.1 September 14, 2024 - rewritten from minibar.py
    # 0.2 November 10, 2024 - support for multiple primer pairs
    # 0.3 December 4, 2024 - code & doc cleanup, write pooling
    # 0.4 February 1, 2025 - bloom filter acceleration
    return "specimux.py version 0.4"

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Specimux: Demultiplex MinION sequences by dual barcode indexes and primers.")

    parser.add_argument("primer_file", help="Fasta file containing primer information")
    parser.add_argument("specimen_file", help="TSV file containing specimen mapping with barcodes and primers")
    parser.add_argument("sequence_file", help="Sequence file in Fasta or Fastq format, gzipped or plain text")

    parser.add_argument("--min-length", type=int, default=-1, help="Minimum sequence length.  Shorter sequences will be skipped (default: no filtering)")
    parser.add_argument("--max-length", type=int, default=-1, help="Maximum sequence length.  Longer sequences will be skipped (default: no filtering)")
    parser.add_argument("-n", "--num-seqs", type=str, default="-1", help="Number of sequences to read from file (e.g., -n 100 or -n 102,3)")
    parser.add_argument("-e", "--index-edit-distance", type=int, default=-1, help="Barcode edit distance value, default is half of min distance between barcodes")
    parser.add_argument("-E", "--primer-edit-distance", type=int, default=-1, help="Primer edit distance value, default is min distance between primers")
    parser.add_argument("-l", "--search-len", type=int, default=80, help="Length to search for index and primer at start and end of sequence (default: 80)")
    parser.add_argument("-F", "--output-to-files", action="store_true", help="Create individual sample files for sequences")
    parser.add_argument("-P", "--output-file-prefix", default="sample_", help="Prefix for individual files when using -F (default: sample_)")
    parser.add_argument("-O", "--output-dir", default=".", help="Directory for individual files when using -F (default: .)")
    parser.add_argument("--color", action="store_true", help="Highlight barcode matches in blue, primer matches in green")
    parser.add_argument("--trim", choices=[TrimMode.NONE, TrimMode.TAILS, TrimMode.BARCODES, TrimMode.PRIMERS], default=TrimMode.BARCODES, help="trimming to apply")
    parser.add_argument("-d", "--diagnostics", action="store_true", help="Output extra diagnostics")
    parser.add_argument("-D", "--debug", action="store_true", help="Enable debug logging")

    parser.add_argument("--disable-prefilter", action="store_true", help="Disable barcode prefiltering (bloom filter optimization)")
    parser.add_argument("--disable-preorient", action="store_true", help="Disable heuristic pre-orientation")
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


# Global variables for worker processes
_barcode_prefilter = PassthroughPrefilter()
_output_manager = None

def init_worker(specimens: Specimens, max_distance: int, args: argparse.Namespace):
    """Initialize worker process with shared resources"""
    global _output_manager, _barcode_prefilter
    try:
        l = logging.DEBUG if args.debug else logging.INFO
        logging.basicConfig(level=l, format='%(asctime)s - %(levelname)s - %(message)s')

        if not args.disable_prefilter:
            barcodes = barcodes_for_bloom_prefilter(specimens)
            cache_path = BloomPrefilter.get_cache_path(barcodes, max_distance)
            _barcode_prefilter = BloomPrefilter.load_readonly(cache_path, barcodes, max_distance)

        # Create output manager for this worker
        if args.output_to_files:
            _output_manager = OutputManager(args.output_dir, args.output_file_prefix,
                                            args.isfastq, max_open_files=50, buffer_size=100)
            _output_manager.__enter__()

    except Exception as e:
        logging.error(f"Failed to initialize worker: {e}")
        raise

def worker(work_item: WorkItem, specimens: Specimens, args: argparse.Namespace):
    """Process a batch of sequences and write results directly"""
    global _output_manager, _barcode_prefilter
    try:
        write_ops, classifications = process_sequences(
            work_item.seq_records, work_item.parameters, specimens, args, _barcode_prefilter)

        # Write sequences directly from worker
        if args.output_to_files and _output_manager is not None:
            try:
                for write_op in write_ops:
                    output_write_operation(write_op, _output_manager, args)
                # Ensure buffer is flushed after each batch
                _output_manager.file_manager.flush_all()
            except Exception as e:
                logging.error(f"Error writing output: {e}")
                raise

        return classifications

    except Exception as e:
        logging.error(traceback.format_exc())
        # On error, try to clean up immediately rather than waiting for atexit
        if _output_manager is not None:
            try:
                logging.debug("Worker cleaning up output manager")
                _output_manager.__exit__(None, None, None)
            except Exception as e:
                logging.error(f"Error cleaning up worker output manager: {e}")
        raise WorkerException(e)

def barcodes_for_bloom_prefilter(specimens):
    # Collect barcodes
    all_b1s = set()
    all_b2s = set()
    for primer in specimens.get_primers(Primer.FWD):
        all_b1s.update(primer.barcodes)
    for primer in specimens.get_primers(Primer.REV):
        all_b2s.update(primer.barcodes)
    barcode_rcs = [reverse_complement(b) for b in all_b1s] + [reverse_complement(b) for b in all_b2s]
    return barcode_rcs


def estimate_sequence_count(filename: str, args: argparse.Namespace) -> int:
    """Estimate total sequences based on file size and sampling"""
    if args.num_seqs > 0:
        return args.num_seqs

    # Get file size
    file_size = os.path.getsize(filename)

    # Sample first 1000 sequences to get average record size
    sample_size = 1000
    seq_records = open_sequence_file(filename, args)
    total_bytes = 0
    for i, record in enumerate(itertools.islice(seq_records, sample_size)):
        total_bytes += len(str(record.seq))
        if args.isfastq:
            total_bytes += len(record.id) + len(record.description) + 2  # +2 for @ and newlines
            total_bytes += len(record.letter_annotations["phred_quality"]) + 3  # +3 for + and newlines
        else:
            total_bytes += len(record.id) + len(record.description) + 2  # +2 for > and newline

    if i == 0:
        return 0

    avg_record_size = total_bytes / (i + 1)
    raw_estimate = file_size / avg_record_size

    # Round to 2 significant figures, rounding up
    if raw_estimate > 0:
        raw_estimate *= 1.05 # fudge to avoid disappointment
        magnitude = math.floor(math.log10(raw_estimate))
        scaled = raw_estimate / (10 ** (magnitude - 1))
        rounded = math.ceil(scaled) * (10 ** (magnitude - 1))
        estimated_sequences = int(rounded)
    else:
        estimated_sequences = 0

    logging.info(f"Estimated {estimated_sequences:,} sequences based on file size")

    return estimated_sequences


def cleanup_locks(output_dir: str):
    """Remove the shared lock directory after all workers are done"""
    lock_dir = os.path.join(output_dir, ".specimux_locks")
    try:
        shutil.rmtree(lock_dir, ignore_errors=True)
    except Exception as e:
        logging.warning(f"Error cleaning up lock directory: {e}")


def specimux_mp(args):
    primer_registry = read_primers_file(args.primer_file)
    specimens = read_specimen_file(args.specimen_file, primer_registry)
    specimens.validate()
    parameters = setup_match_parameters(args, specimens)

    total_seqs = estimate_sequence_count(args.sequence_file, args)
    seq_records = open_sequence_file(args.sequence_file, args)

    create_output_files(args, specimens)

    start_time = timeit.default_timer()
    sequence_block_size = 1000
    last_seq_to_output = args.num_seqs
    all_seqs = last_seq_to_output < 0

    # Skip to start_seq if necessary
    if args.start_seq > 1:
        for _ in itertools.islice(seq_records, args.start_seq - 1):
            pass

    classifications = Counter()

    num_processes = args.threads if args.threads > 0 else multiprocessing.cpu_count()
    logging.info(f"Will run {num_processes} worker processes")

    with Pool(processes=num_processes,
              initializer=init_worker,
              initargs=(specimens, parameters.max_dist_index, args)) as pool:

        worker_func = partial(worker, specimens=specimens, args=args)

        try:
            pbar = tqdm(total=total_seqs, desc="Processing sequences", unit="seq")
            for batch_class in pool.imap_unordered(
                    worker_func,
                    (WorkItem(i, batch, parameters)
                     for i, batch in enumerate(iter_batches(seq_records, sequence_block_size,
                                                            last_seq_to_output, all_seqs)))):
                classifications += batch_class

                # Update progress
                batch_size = sum(batch_class.values())
                pbar.update(batch_size)

                # Update progress description
                total_processed = sum(classifications.values())
                matched = classifications[MatchCode.MATCHED]
                match_rate = matched / total_processed if total_processed > 0 else 0
                pbar.set_description(f"Processing sequences [Match rate: {match_rate:.1%}]")

            pbar.close()

        except WorkerException as e:
            logging.error(f"Unexpected error in worker (see details above): {e}")
            sys.exit(1)

    elapsed = timeit.default_timer() - start_time
    logging.info(f"Elapsed time: {elapsed:.2f} seconds")
    output_diagnostics(args, classifications, elapsed)

    cleanup_locks(args.output_dir)


def write_primers_fasta(output_dir: str, fwd_primer: PrimerInfo, rev_primer: PrimerInfo):
    """
    Write a primers.fasta file containing the forward and reverse primers to the specified directory.

    Args:
        output_dir: Directory to write the primers.fasta file to
        fwd_primer: Forward primer info
        rev_primer: Reverse primer info
    """
    primer_file_path = os.path.join(output_dir, "primers.fasta")

    with open(primer_file_path, 'w') as f:
        # Write forward primer
        f.write(f">{fwd_primer.name} position=forward pool={','.join(fwd_primer.pools)}\n")
        f.write(f"{fwd_primer.primer}\n")

        # Write reverse primer
        f.write(f">{rev_primer.name} position=reverse pool={','.join(rev_primer.pools)}\n")
        f.write(f"{rev_primer.primer}\n")


def create_output_files(args, specimens):
    """Create directory structure for output files and add primers.fasta files"""
    if args.output_to_files:
        # Create base output directory
        os.makedirs(args.output_dir, exist_ok=True)

        # Create unknown directory for complete unknowns
        os.makedirs(os.path.join(args.output_dir, "unknown", "unknown-unknown"), exist_ok=True)

        # Create pool directories and their substructure
        for pool in specimens._primer_registry.get_pools():
            pool_dir = os.path.join(args.output_dir, pool)

            # Create pool directory
            os.makedirs(pool_dir, exist_ok=True)

            # Create primer pair directories
            for fwd_primer in specimens._primer_registry.get_pool_primers(pool, Primer.FWD):
                for rev_primer in specimens._primer_registry.get_pool_primers(pool, Primer.REV):
                    primer_dir = f"{fwd_primer.name}-{rev_primer.name}"
                    primer_full_dir = os.path.join(pool_dir, primer_dir)

                    # Create directories for all match types
                    os.makedirs(primer_full_dir, exist_ok=True)
                    os.makedirs(os.path.join(primer_full_dir, "full"), exist_ok=True)
                    os.makedirs(os.path.join(primer_full_dir, "partial"), exist_ok=True)
                    os.makedirs(os.path.join(primer_full_dir, "ambiguous"), exist_ok=True)
                    os.makedirs(os.path.join(primer_full_dir, "unknown"), exist_ok=True)

                    # Write primers.fasta file for this primer pair directory and the full subdirectory
                    write_primers_fasta(primer_full_dir, fwd_primer, rev_primer)
                    write_primers_fasta(os.path.join(primer_full_dir, "full"), fwd_primer, rev_primer)

                # Create unknown primer directories
                unknown_primer_dir = f"{fwd_primer.name}-unknown"
                os.makedirs(os.path.join(args.output_dir, pool, unknown_primer_dir), exist_ok=True)

            # Create unknown primer directories for reverse only matches
            for rev_primer in specimens._primer_registry.get_pool_primers(pool, Primer.REV):
                unknown_primer_dir = f"unknown-{rev_primer.name}"
                os.makedirs(os.path.join(args.output_dir, pool, unknown_primer_dir), exist_ok=True)

def iter_batches(seq_records, batch_size: int, max_seqs: int, all_seqs: bool):
    """Helper to iterate over sequence batches"""
    num_seqs = 0
    while all_seqs or num_seqs < max_seqs:
        to_read = batch_size if all_seqs else min(batch_size, max_seqs - num_seqs)
        batch = list(itertools.islice(seq_records, to_read))
        if not batch:
            break
        yield batch
        num_seqs += len(batch)

def specimux(args):
    primer_registry = read_primers_file(args.primer_file)
    specimens = read_specimen_file(args.specimen_file, primer_registry)
    specimens.validate()
    parameters = setup_match_parameters(args, specimens)

    total_seqs = estimate_sequence_count(args.sequence_file, args)
    seq_records = open_sequence_file(args.sequence_file, args)

    create_output_files(args, specimens)

    start_time = timeit.default_timer()

    sequence_block_size = 1000
    last_seq_to_output = args.num_seqs
    all_seqs = last_seq_to_output < 0

    # Skip to the start_seq if necessary
    if args.start_seq > 1:
        for _ in itertools.islice(seq_records, args.start_seq - 1):
            pass

    classifications = Counter()

    with OutputManager(args.output_dir, args.output_file_prefix, args.isfastq) as output_manager:
        num_seqs = 0
        prefilter = PassthroughPrefilter()
        if not args.disable_prefilter:
            barcode_rcs = barcodes_for_bloom_prefilter(specimens)
            cache_path = BloomPrefilter.get_cache_path(barcode_rcs, parameters.max_dist_index)
            prefilter = BloomPrefilter.load_readonly(cache_path, barcode_rcs, parameters.max_dist_index)

        pbar = tqdm(total=total_seqs, desc="Processing sequences", unit="seq")
        while all_seqs or num_seqs < last_seq_to_output:
            to_read = sequence_block_size if all_seqs else min(sequence_block_size, last_seq_to_output - num_seqs)
            seq_batch = list(itertools.islice(seq_records, to_read))
            if not seq_batch:
                break

            write_ops, batch_classifications = process_sequences(seq_batch, parameters, specimens,
                                                                                  args, prefilter)
            classifications += batch_classifications
            for write_op in write_ops:
                output_write_operation(write_op, output_manager, args)

            batch_size = len(seq_batch)
            num_seqs += batch_size
            pbar.update(batch_size)

            # Update progress bar description with match rate
            total_processed = sum(classifications.values())
            matched = classifications[MatchCode.MATCHED]
            match_rate = matched / total_processed if total_processed > 0 else 0
            pbar.set_description(f"Processing sequences [Match rate: {match_rate:.1%}]")

        pbar.close()

    elapsed = timeit.default_timer() - start_time
    logging.info(f"Elapsed time: {elapsed:.2f} seconds")
    output_diagnostics(args, classifications, elapsed)


def output_diagnostics(args, classifications, elapsed_time=None):
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
    
    # Generate the classification statistics
    stats_lines = ["Classification Statistics:"]
    for classification, count in sorted_classifications:
        stats_line = f"{classification:<{max_name_length}} : {count:>{max_count_length}} ({count / total:2.2%})"
        stats_lines.append(stats_line)
        
    # Log to console if diagnostics flag is set
    if args.diagnostics:
        for line in stats_lines:
            logging.info(line)
    
    # Write to log.txt regardless of diagnostics flag
    log_file_path = os.path.join(args.output_dir, "log.txt")
    matched_count = classifications.get(MatchCode.MATCHED, 0)
    match_rate = matched_count / total if total > 0 else 0
    
    with open(log_file_path, 'w') as log_file:
        log_file.write(f"Specimux Run Summary\n")
        log_file.write(f"==================\n\n")
        log_file.write(f"Date: {timeit.time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        log_file.write(f"Command: {' '.join(sys.argv)}\n\n")
        
        log_file.write(f"Input Files:\n")
        log_file.write(f"  Primer File: {args.primer_file}\n")
        log_file.write(f"  Specimen File: {args.specimen_file}\n")
        log_file.write(f"  Sequence File: {args.sequence_file}\n\n")
        
        log_file.write(f"Run Statistics:\n")
        log_file.write(f"  Total Sequences Processed: {total:,}\n")
        log_file.write(f"  Successfully Matched: {matched_count:,} ({match_rate:.2%})\n")
        if elapsed_time is not None:
            seqs_per_sec = total / elapsed_time if elapsed_time > 0 else 0
            log_file.write(f"  Processing Time: {elapsed_time:.2f} seconds\n")
            log_file.write(f"  Processing Rate: {seqs_per_sec:.2f} sequences/second\n\n")
        
        log_file.write("Detailed Classification Statistics:\n")
        for line in stats_lines[1:]:  # Skip the header line
            log_file.write(f"  {line}\n")
            
    logging.info(f"Summary statistics written to {log_file_path}")

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

    def _bp_adjusted_length(primer):
        score = 0
        for b in primer:
            if b in ['A', 'C', 'G', 'T']: score += 3
            elif b in ['K', 'M', 'R', 'S', 'W', 'Y']: score += 2
            elif b in ['B', 'D', 'H', 'V']: score += 1
        return score / 3.0

    # Collect all barcodes across primer pairs
    all_b1s = set()
    all_b2s = set()
    for primer in specimens.get_primers(Primer.FWD):
        all_b1s.update(primer.barcodes)
    for primer in specimens.get_primers(Primer.REV):
        all_b2s.update(primer.barcodes)

    _sanity_check_distance(list(all_b1s), "Forward Barcodes", args)
    _sanity_check_distance(list(all_b2s), "Reverse Barcodes", args)

    combined_barcodes = list(all_b1s) + [reverse_complement(b) for b in all_b2s]
    min_bc_dist = _sanity_check_distance(combined_barcodes,
        "Forward Barcodes + Reverse Complement of Reverse Barcodes", args)

    primers = []
    for pi in specimens.get_primers(Primer.FWD):
        primers.append(pi.primer)
        primers.append(pi.primer_rc)
    for pi in specimens.get_primers(Primer.REV):
        primers.append(pi.primer)
        primers.append(pi.primer_rc)
    primers = list(set(primers))
    min_primer_dist = _sanity_check_distance(primers, "All Primers and Reverse Complements", args)

    max_search_area = args.search_len
    max_dist_index = math.ceil(min_bc_dist / 2.0)
    if args.index_edit_distance != -1:
        max_dist_index = args.index_edit_distance

    logging.info(f"Using Edit Distance Thresholds {max_dist_index} for barcode indexes")

    primer_thresholds = {}
    for pi in specimens.get_primers(Primer.FWD):
        if args.primer_edit_distance != -1:
            primer_thresholds[pi.primer] = args.primer_edit_distance
        else:
            primer_thresholds[pi.primer] = int(_bp_adjusted_length(pi.primer) / 3)
    for pi in specimens.get_primers(Primer.REV):
        if args.primer_edit_distance != -1:
            primer_thresholds[pi.primer] = args.primer_edit_distance
        else:
            primer_thresholds[pi.primer] = int(_bp_adjusted_length(pi.primer) / 3)

    for p, pt in primer_thresholds.items():
        logging.info(f"Using Edit Distance Threshold {pt} for primer {p}")

    if args.disable_preorient:
        logging.info("Sequence pre-orientation disabled, may run slower")
        preorient = False
    else:
        preorient = True

    parameters = MatchParameters(primer_thresholds, max_dist_index, max_search_area, preorient)

    if not args.disable_prefilter:
        if specimens.b_length() > 13:
            logging.warning("Barcode prefilter not tested for barcodes longer than 13 nt.  You may need to use --disable-prefilter")
        if max_dist_index > 3:
            logging.warning("Barcode prefilter not tested for edit distance greater than 3.  You may need to use --disable-prefilter")
        barcode_rcs = barcodes_for_bloom_prefilter(specimens)
        BloomPrefilter.create_filter(barcode_rcs, parameters.max_dist_index)
        logging.info("Using Bloom Filter optimization for barcode matching")
    else:
        logging.info("Barcode prefiltering disabled, may run slower")

    return parameters

def main(argv):
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    if args.output_to_files:
        specimux_mp(args)  # Use multiprocess for file output
    else:
        if args.threads > 1:
            logging.warning(f"Multithreading only supported for file output. Ignoring --threads {args.threads}")
        specimux(args)     # Use single process for console output

if __name__ == "__main__":
    main(sys.argv)


