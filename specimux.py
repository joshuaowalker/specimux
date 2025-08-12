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
from collections import Counter
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from functools import partial
from multiprocessing import Pool
from operator import itemgetter
from pathlib import Path
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
import copy
import subprocess
import re
import os.path

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
    UNKNOWN = "unknown"
    PREFIX_FWD_MATCH = "barcode_fwd_"
    PREFIX_REV_MATCH = "barcode_rev_"

class TrimMode:
    PRIMERS = "primers"
    BARCODES = "barcodes"
    TAILS = "tails"
    NONE = "none"

class MultipleMatchStrategy:
    """Strategy for handling multiple equivalent matches"""
    RETAIN = "retain"  # Default: output all equivalent matches
    DOWNGRADE_FULL = "downgrade-full"  # Downgrade full matches to partial if multiple exist


class ResolutionType(Enum):
    """Types of specimen resolution outcomes."""
    FULL_MATCH = 1  # Both primers + both barcodes matched to a specimen
    PARTIAL_FORWARD = 2  # Forward barcode only (generates FWD_ONLY_ sample ID)
    PARTIAL_REVERSE = 3  # Reverse barcode only (generates REV_ONLY_ sample ID)
    DOWNGRADED_MULTIPLE = 4  # Full match downgraded due to multiple equivalent matches
    MULTIPLE_SPECIMENS = 5  # Multiple specimen matches (shouldn't happen in new flow)
    UNKNOWN = 6  # No resolution possible

    def to_string(self) -> str:
        """Convert resolution type to lowercase string for trace logging."""
        if self == ResolutionType.FULL_MATCH:
            return 'full_match'
        elif self == ResolutionType.PARTIAL_FORWARD:
            return 'partial_forward'
        elif self == ResolutionType.PARTIAL_REVERSE:
            return 'partial_reverse'
        elif self == ResolutionType.DOWNGRADED_MULTIPLE:
            return 'downgraded_multiple_full'
        elif self == ResolutionType.MULTIPLE_SPECIMENS:
            return 'multiple_specimens'
        else:
            return 'unknown'

    def is_full_match(self) -> bool:
        """Check if this represents a successful full match."""
        return self == ResolutionType.FULL_MATCH

    def is_partial_match(self) -> bool:
        """Check if this represents a partial match."""
        return self in [ResolutionType.PARTIAL_FORWARD, ResolutionType.PARTIAL_REVERSE]

    def is_unknown(self) -> bool:
        """Check if this represents an unknown/failed resolution."""
        return self == ResolutionType.UNKNOWN


class Barcode(Enum):
    B1 = 1
    B2 = 2
    
    def to_string(self) -> str:
        """Convert barcode to lowercase string for trace logging."""
        if self == Barcode.B1:
            return 'forward'
        elif self == Barcode.B2:
            return 'reverse'
        else:
            return 'unknown'

class Primer(Enum):
    FWD = 3
    REV = 4
    
    def to_string(self) -> str:
        """Convert primer to lowercase string for trace logging."""
        if self == Primer.FWD:
            return 'forward'
        elif self == Primer.REV:
            return 'reverse'
        else:
            return 'unknown'


class PrimerInfo:
    def __init__(self, name: str, seq: str, direction: Primer, pools: List[str]):
        self.name = name
        self.primer = seq.upper()
        self.direction = direction
        self.primer_rc = reverse_complement(seq.upper())
        self.barcodes = set()
        self.specimens = set()
        self.pools = pools

class PrimerDatabase:
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


class AlignmentResult:
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
        m = AlignmentResult(self._edlib_match.copy())
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

class CandidateMatch:
    def __init__(self, sequence, barcode_length, candidate_match_id: Optional[str] = None):
        self.sequence_length = len(sequence)
        self.sequence = sequence
        self.p1_match: Optional[AlignmentResult] = None
        self._b1_matches: List[Tuple[str, AlignmentResult, float]] = []
        self.p2_match: Optional[AlignmentResult] = None
        self._b2_matches: List[Tuple[str, AlignmentResult, float]] = []
        self._p1: Optional[PrimerInfo] = None
        self._p2: Optional[PrimerInfo] = None
        self._pool: Optional[str] = None  # New: track the pool
        self.match_tolerance = 1.0
        self._barcode_length = barcode_length
        self.candidate_match_id = candidate_match_id  # New: track candidate match ID

    def set_pool(self, pool: str):
        """Set the primer pool for this match"""
        self._pool = pool

    def get_pool(self) -> Optional[str]:
        """Get the primer pool for this match"""
        return self._pool

    def add_barcode_match(self, match: AlignmentResult, barcode: str, reverse: bool, which: Barcode):
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
        return [b for b, _, d in self._b1_matches if abs(d - best_distance) < self.match_tolerance]

    def best_b2(self) -> List[str]:
        if not self._b2_matches:
            return []
        best_distance = self.b2_distance()
        return [b for b, _, d in self._b2_matches if abs(d - best_distance) < self.match_tolerance]

    def set_primer_match(self, match: AlignmentResult, primer: PrimerInfo, reverse: bool, which: Primer):
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
    def __init__(self, primer_registry: PrimerDatabase):
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
    resolution_type: ResolutionType

class SequenceBatch(NamedTuple):
    seq_number: int
    seq_records: List  # This now contains actual sequence records, not an iterator
    parameters: MatchParameters
    start_idx: int     # Starting index for sequence ID generation

@dataclass

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

        # Prepare the data outside of any locks
        buffer_data = ''.join(self.write_buffers[filename])
        self.write_buffers[filename].clear()

        # Get or open the file handle (no locking needed for this step)
        try:
            f = self.file_cache[filename]
        except KeyError:
            # Create directory if needed
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            f = open(filename, 'a')  # Always append mode
            self.file_cache[filename] = f

        # Get or create lock (only for this file)
        if filename not in self.locks:
            lock_path = self._get_lock_path(filename)
            fd = os.open(lock_path, os.O_CREAT | os.O_RDWR, 0o666)
            self.locks[filename] = fd

        # Now acquire lock, write data, and release lock
        fcntl.lockf(self.locks[filename], fcntl.LOCK_EX)
        try:
            f.write(buffer_data)
            f.flush()
            os.fsync(f.fileno())  # Ensure data is written to disk
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

    def _make_filename(self, sample_id: str, pool: str, p1: str, p2: str, resolution_type: ResolutionType) -> str:
        """Create a filename with match-type-first organization."""
        extension = '.fastq' if self.is_fastq else '.fasta'

        # Handle None values
        sample_id = sample_id or SampleId.UNKNOWN
        pool = pool or "unknown"
        p1 = p1 or "unknown"
        p2 = p2 or "unknown"

        safe_id = "".join(c if c.isalnum() or c in "._-$#" else "_" for c in sample_id)
        primer_dir = f"{p1}-{p2}"

        # Determine match type and construct path with match type first
        if resolution_type.is_unknown():
            # Unknown resolution
            return os.path.join(self.output_dir, "unknown", pool, primer_dir, f"{self.prefix}{safe_id}{extension}")
        elif resolution_type.is_partial_match():
            # Partial match (forward or reverse barcode only)
            return os.path.join(self.output_dir, "partial", pool, primer_dir, f"{self.prefix}{safe_id}{extension}")
        else:
            # Full match (including downgraded and multiple specimen matches)
            return os.path.join(self.output_dir, "full", pool, primer_dir, f"{self.prefix}{safe_id}{extension}")

    def write_sequence(self, write_op: WriteOperation, trace_logger: Optional['TraceLogger'] = None):
        """Write a sequence to the appropriate output file."""
        filename = self._make_filename(write_op.sample_id, write_op.primer_pool,
                                       write_op.p1_name, write_op.p2_name, write_op.resolution_type)

        # Log sequence output with actual filename
        if trace_logger:
            # Create relative path from output directory
            relative_path = os.path.relpath(filename, self.output_dir)
            primer_pair = f"{write_op.p1_name}-{write_op.p2_name}"
            
            trace_logger.log_sequence_output(write_op.seq_id, write_op.sample_id,
                                           write_op.primer_pool, primer_pair, relative_path)

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

        output_content = ''.join(output)
        self.file_manager.write(filename, output_content)
        
        # If this is a full match, write it to the pool-level aggregation directory as well
        if write_op.resolution_type.is_full_match():
            
            # Create additional path for pool-level aggregation in full directory
            extension = '.fastq' if self.is_fastq else '.fasta'
            safe_id = "".join(c if c.isalnum() or c in "._-$#" else "_" for c in write_op.sample_id)
            pool_full_path = os.path.join(self.output_dir, "full", write_op.primer_pool, 
                                          f"{self.prefix}{safe_id}{extension}")
            
            # Ensure the pool full directory exists
            os.makedirs(os.path.dirname(pool_full_path), exist_ok=True)
            
            # Write to pool-level aggregation directory
            self.file_manager.write(pool_full_path, output_content)

class WorkerException(Exception):
    pass


class TraceLogger:
    """Manages trace event logging for sequence processing pipeline."""
    
    def __init__(self, enabled: bool, verbosity: int, output_dir: str, worker_id: str, 
                 start_timestamp: str, buffer_size: int = 1000):
        """Initialize trace logger.
        
        Args:
            enabled: Whether trace logging is enabled
            verbosity: Verbosity level (1-3)
            output_dir: Directory for trace files
            worker_id: Unique worker identifier
            start_timestamp: Timestamp for file naming
            buffer_size: Number of events to buffer before writing
        """
        self.enabled = enabled
        self.verbosity = verbosity
        self.worker_id = worker_id
        self.event_counter = 0
        self.buffer = []
        self.buffer_size = buffer_size
        self.file_handle = None
        self.sequence_record_counter = 0
        
        if self.enabled:
            # Create trace directory
            trace_dir = Path(output_dir) / "trace"
            trace_dir.mkdir(exist_ok=True)
            
            # Create trace file
            filename = f"specimux_trace_{start_timestamp}_{worker_id}.tsv"
            self.filepath = trace_dir / filename
            self.file_handle = open(self.filepath, 'w', newline='')
            self.writer = csv.writer(self.file_handle, delimiter='\t')
            
            # Write header
            header = ['timestamp', 'worker_id', 'event_seq', 'sequence_id', 'event_type']
            self.writer.writerow(header)
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    def close(self):
        """Flush buffer and close file."""
        if self.enabled and self.file_handle:
            try:
                # Always flush any remaining events in buffer
                if self.buffer:
                    self._flush_buffer()
                self.file_handle.flush()
                self.file_handle.close()
                self.file_handle = None
                logging.debug(f"Trace logger {self.worker_id} closed successfully")
            except Exception as e:
                logging.error(f"Error closing trace logger {self.worker_id}: {e}")
    
    def _flush_buffer(self):
        """Write buffered events to file."""
        if self.file_handle and self.buffer:
            try:
                buffer_size = len(self.buffer)
                self.writer.writerows(self.buffer)
                self.file_handle.flush()
                self.buffer = []
                if buffer_size > 0:
                    logging.debug(f"Trace logger {self.worker_id} flushed {buffer_size} events")
            except Exception as e:
                logging.error(f"Error flushing trace logger {self.worker_id}: {e}")
    
    def _log_event(self, sequence_id: str, event_type: str, *fields):
        """Log a trace event."""
        if not self.enabled:
            return
        
        self.event_counter += 1
        timestamp = datetime.now().isoformat()
        row = [timestamp, self.worker_id, self.event_counter, sequence_id, event_type] + list(fields)
        self.buffer.append(row)
        
        if len(self.buffer) >= self.buffer_size:
            self._flush_buffer()
    
    def get_sequence_id(self, seq_record, record_num: Optional[int] = None) -> str:
        """Create unique sequence ID."""
        if record_num is None:
            self.sequence_record_counter += 1
            record_num = self.sequence_record_counter
        return f"{seq_record.id}#{record_num:08d}#{self.worker_id}"
    
    # Standard events (verbosity level 1)
    
    def log_sequence_received(self, sequence_id: str, sequence_length: int, sequence_name: str):
        """Log when a sequence enters the pipeline."""
        self._log_event(sequence_id, 'SEQUENCE_RECEIVED', sequence_length, sequence_name)
    
    def log_sequence_filtered(self, sequence_id: str, sequence_length: int, filter_reason: str):
        """Log when a sequence is filtered out."""
        self._log_event(sequence_id, 'SEQUENCE_FILTERED', sequence_length, filter_reason)
    
    def log_orientation_detected(self, sequence_id: str, orientation: str, 
                                 forward_score: int, reverse_score: int, confidence: float):
        """Log orientation detection result."""
        self._log_event(sequence_id, 'ORIENTATION_DETECTED', orientation, 
                       forward_score, reverse_score, f"{confidence:.3f}")
    
    def log_primer_matched(self, sequence_id: str, match: 'CandidateMatch',
                           pool: str, orientation_used: str):
        """Log successful primer match for a specific candidate match."""
        (candidate_match_id, p1_name, p2_name, _, _, 
         _, _, p1_dist, p2_dist, _, _) = self._extract_match_info(match)
        
        # Determine match type from SequenceMatch
        if match.p1_match and match.p2_match:
            match_type = 'both'
        elif match.p1_match:
            match_type = 'forward_only'
        else:
            match_type = 'reverse_only'
        
        self._log_event(sequence_id, 'PRIMER_MATCHED', candidate_match_id, match_type,
                       p1_name, p2_name, p1_dist, p2_dist,
                       pool, orientation_used)
    
    def log_barcode_matched(self, sequence_id: str, match: 'CandidateMatch'):
        """Log barcode match result for a specific candidate match."""
        (candidate_match_id, p1_name, p2_name, b1_name, b2_name, 
         _, _, _, _, b1_dist, b2_dist) = self._extract_match_info(match)
        
        # Determine barcode match type from SequenceMatch
        if match.has_b1_match() and match.has_b2_match():
            barcode_match_type = 'both'
        elif match.has_b1_match():
            barcode_match_type = 'forward_only'
        elif match.has_b2_match():
            barcode_match_type = 'reverse_only'
        else:
            barcode_match_type = 'none'
        
        self._log_event(sequence_id, 'BARCODE_MATCHED', candidate_match_id, barcode_match_type,
                       b1_name, b2_name, b1_dist, b2_dist,
                       p1_name, p2_name)
    
    def _extract_match_info(self, match: 'CandidateMatch') -> tuple:
        """Extract common information from SequenceMatch for trace logging.
        
        Returns: (candidate_match_id, p1_name, p2_name, b1_name, b2_name, 
                 barcode_presence, total_distance, p1_dist, p2_dist, b1_dist, b2_dist)
        """
        candidate_match_id = match.candidate_match_id or 'unknown'
        p1_name = match.get_p1().name if match.get_p1() else 'none'
        p2_name = match.get_p2().name if match.get_p2() else 'none'
        b1_name = match.best_b1()[0] if match.best_b1() else 'none'
        b2_name = match.best_b2()[0] if match.best_b2() else 'none'
        
        # Determine barcode presence
        if match.has_b1_match() and match.has_b2_match():
            barcode_presence = 'both'
        elif match.has_b1_match():
            barcode_presence = 'forward_only'
        elif match.has_b2_match():
            barcode_presence = 'reverse_only'
        else:
            barcode_presence = 'none'
        
        # Calculate distances
        total_distance = 0
        p1_dist = -1
        p2_dist = -1
        b1_dist = -1
        b2_dist = -1
        
        if match.p1_match:
            p1_dist = match.p1_match.distance()
            total_distance += p1_dist
        if match.p2_match:
            p2_dist = match.p2_match.distance()
            total_distance += p2_dist
        if match.has_b1_match():
            b1_dist = match.b1_distance()
            total_distance += b1_dist
        if match.has_b2_match():
            b2_dist = match.b2_distance()
            total_distance += b2_dist
        
        return (candidate_match_id, p1_name, p2_name, b1_name, b2_name, 
                barcode_presence, total_distance, p1_dist, p2_dist, b1_dist, b2_dist)
    
    def log_match_scored(self, sequence_id: str, match: 'CandidateMatch', score: float):
        """Log match scoring for a specific candidate match."""
        (candidate_match_id, p1_name, p2_name, b1_name, b2_name, 
         barcode_presence, total_distance, _, _, _, _) = self._extract_match_info(match)
        self._log_event(sequence_id, 'MATCH_SCORED', candidate_match_id, p1_name, p2_name,
                       b1_name, b2_name, total_distance,
                       barcode_presence, f"{score:.3f}")
    

    
    def log_match_selected(self, sequence_id: str, selection_strategy: str,
                          forward_primer: str, reverse_primer: str,
                          forward_barcode: str, reverse_barcode: str,
                          pool: str, is_unique: bool):
        """Log match selection decision."""
        self._log_event(sequence_id, 'MATCH_SELECTED', selection_strategy,
                       forward_primer, reverse_primer, forward_barcode, reverse_barcode,
                       pool, str(is_unique).lower())
    
    def log_specimen_resolved(self, sequence_id: str, match: 'CandidateMatch',
                              specimen_id: str, resolution_type: str, pool: str):
        """Log specimen resolution."""
        (_, p1_name, p2_name, b1_name, b2_name, 
         _, _, _, _, _, _) = self._extract_match_info(match)
        self._log_event(sequence_id, 'SPECIMEN_RESOLVED', specimen_id, resolution_type,
                       pool, p1_name, p2_name, b1_name, b2_name)
    
    def log_sequence_output(self, sequence_id: str, specimen_id: str,
                           pool: str, primer_pair: str, file_path: str):
        """Log output decision."""
        self._log_event(sequence_id, 'SEQUENCE_OUTPUT', specimen_id,
                       pool, primer_pair, file_path)
    
    def log_no_match_found(self, sequence_id: str, stage_failed: str, reason: str):
        """Log when no matches found."""
        self._log_event(sequence_id, 'NO_MATCH_FOUND', stage_failed, reason)
    
    def log_match_discarded(self, sequence_id: str, match: 'CandidateMatch',
                            score: float, discard_reason: str):
        """Log when a candidate match is discarded."""
        (candidate_match_id, p1_name, p2_name, b1_name, b2_name, 
         _, _, _, _, _, _) = self._extract_match_info(match)
        self._log_event(sequence_id, 'MATCH_DISCARDED', candidate_match_id,
                       p1_name, p2_name, b1_name, b2_name,
                       score, discard_reason)
    
    # Detailed events (verbosity level 2+)
    
    def log_primer_search(self, sequence_id: str, primer_name: str, primer_direction: str,
                         search_start: int, search_end: int, found: bool,
                         edit_distance: int, match_position: int):
        """Log primer search attempt (level 2+)."""
        if self.verbosity >= 2 and (self.verbosity >= 3 or found):
            self._log_event(sequence_id, 'PRIMER_SEARCH', primer_name, primer_direction,
                           search_start, search_end, str(found).lower(),
                           edit_distance, match_position)
    
    def log_barcode_search(self, sequence_id: str, barcode_name: str, barcode_type: str,
                          primer_adjacent: str, search_start: int, search_end: int,
                          found: bool, edit_distance: int, match_position: int):
        """Log barcode search attempt (level 3 only)."""
        if self.verbosity >= 3:
            self._log_event(sequence_id, 'BARCODE_SEARCH', barcode_name, barcode_type,
                           primer_adjacent, search_start, search_end, str(found).lower(),
                           edit_distance, match_position)


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
def read_primers_file(filename: str) -> PrimerDatabase:
    """
    Read primers file and build primer registry

    Returns:
        PrimerRegistry object managing all primers and their relationships
    """
    registry = PrimerDatabase()

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

def read_specimen_file(filename: str, primer_registry: PrimerDatabase) -> Specimens:
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


def detect_file_format(filename: str) -> str:
    """
    Detect file format from filename, handling compressed files.

    Args:
        filename: Path to sequence file

    Returns:
        str: Detected format ('fastq', 'fasta', or other future formats)
    """

    # Strip compression extensions recursively (handles cases like .fastq.gz.gz)
    base_name = os.path.basename(filename)
    compression_exts = ['.gz', '.gzip', '.bz2', '.zip']

    # Keep stripping compression extensions until none are left
    root, ext = os.path.splitext(base_name)
    while ext.lower() in compression_exts:
        base_name = root
        root, ext = os.path.splitext(base_name)

    # Check for known file formats (case-insensitive)
    if base_name.lower().endswith(('.fastq', '.fq')):
        return 'fastq'
    elif base_name.lower().endswith(('.fasta', '.fa', '.fna')):
        return 'fasta'
    else:
        # Check the first few bytes of the file if extension doesn't help
        try:
            # Handle compressed files
            if filename.endswith(('.gz', '.gzip')):
                with gzip.open(filename, 'rt') as f:
                    first_char = f.read(1)
            else:
                with open(filename, 'rt') as f:
                    first_char = f.read(1)

            # Check first character
            if first_char == '@':
                return 'fastq'
            elif first_char == '>':
                return 'fasta'
        except Exception:
            pass

        # Default to FASTA if we can't determine
        return 'fasta'


def open_sequence_file(filename, args):
    """
    Open a sequence file, automatically detecting format and compression.

    Args:
        filename: Path to sequence file
        args: Command line arguments, will update args.isfastq based on format

    Returns:
        SeqIO iterator for the file
    """
    is_gzipped = filename.endswith((".gz", ".gzip"))

    # Detect file format
    file_format = detect_file_format(filename)
    args.isfastq = file_format == 'fastq'

    if is_gzipped:
        handle = gzip.open(filename, "rt")  # Open in text mode
        return SeqIO.parse(handle, file_format)
    else:
        return SeqIO.parse(filename, file_format)

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


def create_write_operation(sample_id, args, seq, match, resolution_type):
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
        resolution_type=resolution_type
    )




def process_sequences(seq_records: List[SeqRecord],
                      parameters: MatchParameters,
                      specimens: Specimens,
                      args: argparse.Namespace,
                      prefilter: BarcodePrefilter,
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
                    
                    # Create and add write operation directly
                    op = create_write_operation(final_sample_id, args, seq, match, resolution_type)
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
                op = create_write_operation(SampleId.UNKNOWN, args, seq, match, ResolutionType.UNKNOWN)
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



class Orientation(Enum):
    FORWARD = 1
    REVERSE = 2
    UNKNOWN = 3
    
    def to_string(self) -> str:
        """Convert orientation to lowercase string for trace logging."""
        return self.name.lower()




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
        fwd_p2 = align_seq(primer.primer, seq, parameters.max_dist_primers[primer.primer],
                           len(seq) - parameters.search_len, len(seq))
        rev_p2 = align_seq(primer.primer, rseq, parameters.max_dist_primers[primer.primer],
                           len(seq) - parameters.search_len, len(seq))
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


def find_candidate_matches(prefilter: BarcodePrefilter, parameters: MatchParameters, seq: SeqRecord,
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
                    pool = get_pool_from_primers(match.get_p1(), match.get_p2())
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
                    pool = get_pool_from_primers(match.get_p1(), match.get_p2())
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

def match_one_end(prefilter: BarcodePrefilter, match: CandidateMatch, parameters: MatchParameters, sequence: str,
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
                
                if not prefilter.match(b_rc, target_seq):
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

def output_write_operation(write_op: WriteOperation,
                           output_manager: OutputManager,
                           args: argparse.Namespace,
                           trace_logger: Optional['TraceLogger'] = None) -> None:
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
        output_manager.write_sequence(write_op, trace_logger)

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
    # 0.5 March 19, 2025 - added Primer Pools, Hierarchical Output with pool-level full match collections, and detailed run log
    # 0.6 August 2025 - multiple match processing, comprehensive trace event system, trace-based statistics framework
    return "specimux.py version 0.6.0-dev"

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
    parser.add_argument("--resolve-multiple-matches", choices=[MultipleMatchStrategy.RETAIN, MultipleMatchStrategy.DOWNGRADE_FULL], default=MultipleMatchStrategy.RETAIN, help="Strategy for handling multiple equivalent matches: 'retain' outputs all matches (default), 'downgrade-full' downgrades full matches to partial when multiple exist")
    parser.add_argument("-d", "--diagnostics", nargs='?', const=1, type=int, choices=[1, 2, 3], 
                        help="Enable diagnostic trace logging: 1=standard (default), 2=detailed, 3=verbose")
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
_trace_logger = None

def init_worker(specimens: Specimens, max_distance: int, args: argparse.Namespace, start_timestamp: str = None):
    """Initialize worker process with shared resources"""
    global _output_manager, _barcode_prefilter, _trace_logger
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

        # Create trace logger for this worker if diagnostics enabled
        if args.diagnostics:
            worker_id = f"worker_{multiprocessing.current_process().name.split('-')[-1]}"
            if start_timestamp is None:
                start_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            _trace_logger = TraceLogger(
                enabled=True, 
                verbosity=args.diagnostics,
                output_dir=args.output_dir,
                worker_id=worker_id,
                start_timestamp=start_timestamp
            )
            _trace_logger.__enter__()
        
        # Note: TraceLogger buffers are flushed after each work batch, similar to OutputManager

    except Exception as e:
        logging.error(f"Failed to initialize worker: {e}")
        raise

def cleanup_worker():
    """Clean up worker resources on error."""
    global _output_manager, _trace_logger
    
    if _output_manager is not None:
        try:
            logging.debug("Worker cleaning up output manager")
            _output_manager.__exit__(None, None, None)
        except Exception as e:
            logging.error(f"Error cleaning up worker output manager: {e}")
    
    if _trace_logger is not None:
        try:
            logging.debug("Worker cleaning up trace logger")
            _trace_logger.close()
        except Exception as e:
            logging.error(f"Error cleaning up worker trace logger: {e}")

def worker(work_item: SequenceBatch, specimens: Specimens, args: argparse.Namespace):
    """Process a batch of sequences and write results directly"""
    global _output_manager, _barcode_prefilter, _trace_logger
    try:
        write_ops, total_count, matched_count = process_sequences(
            work_item.seq_records, work_item.parameters, specimens, args, _barcode_prefilter,
            _trace_logger, work_item.start_idx)

        # Write sequences directly from worker
        if args.output_to_files and _output_manager is not None:
            try:
                for write_op in write_ops:
                    output_write_operation(write_op, _output_manager, args, _trace_logger)
                # Ensure buffer is flushed after each batch
                _output_manager.file_manager.flush_all()
            except Exception as e:
                logging.error(f"Error writing output: {e}")
                raise
        
        # Flush trace logger buffer after each batch (same pattern as OutputManager)
        if _trace_logger is not None:
            _trace_logger._flush_buffer()

        # Return counts for progress tracking
        return total_count, matched_count

    except Exception as e:
        logging.error(traceback.format_exc())
        # On error, try to clean up immediately rather than waiting for atexit
        cleanup_worker()
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
    """
    Estimate total sequences in a file, handling both compressed and uncompressed formats.

    Uses consistent logic for both compressed and uncompressed files, with the only
    difference being how total file size is determined.

    This and related functions are an embarrassing amount of code just to figure out the
    100% mark for the progress bar without reading the entire input file.

    Args:
        filename: Path to the sequence file
        args: Command line arguments

    Returns:
        int: Estimated number of sequences in the file
    """
    if args.num_seqs > 0:
        return args.num_seqs

    is_compressed = filename.endswith((".gz", ".gzip"))

    # Step 1: Determine total uncompressed bytes
    if is_compressed:
        try:
            # Get compressed and uncompressed sizes using gzip -l
            compressed_size, uncompressed_size = get_gzip_info(filename)

            if not compressed_size or not uncompressed_size:
                logging.warning("Could not determine compressed/uncompressed size using gzip -l")
                return 500000

            total_bytes = uncompressed_size
            logging.debug(f"File size: {compressed_size:,} bytes compressed, {uncompressed_size:,} bytes uncompressed")
        except Exception as e:
            logging.warning(f"Error getting gzip info: {e}")
            # Fallback to compressed size with a conservative ratio
            total_bytes = os.path.getsize(filename) * 2.5  # Assume 2.5x compression ratio as fallback
            logging.warning(f"Using fallback uncompressed size estimate: {total_bytes:,} bytes")
    else:
        # For uncompressed files, use the file size directly
        total_bytes = os.path.getsize(filename)

    # Step 2: Sample sequences to calculate average bytes per sequence (same for both types)
    try:
        # Sample first 1000 sequences to get average record size
        sample_size = 1000
        seq_records = open_sequence_file(filename, args)
        records_total_bytes = 0
        record_count = 0

        for record in itertools.islice(seq_records, sample_size):
            record_count += 1
            records_total_bytes += len(str(record.seq))
            if args.isfastq:
                records_total_bytes += len(record.id) + len(record.description) + 2  # +2 for @ and newlines
                records_total_bytes += len(record.letter_annotations["phred_quality"]) + 3  # +3 for + and newlines
            else:
                records_total_bytes += len(record.id) + len(record.description) + 2  # +2 for > and newline

        if record_count == 0:
            logging.warning("Failed to read any sequences from the file")
            return 10000  # Return a larger default

        # Calculate average bytes per sequence
        avg_record_size = records_total_bytes / record_count

        # Log actual sampled info for debugging
        logging.debug(
            f"Sampled {record_count} sequences, total bytes: {records_total_bytes:,}, avg bytes per seq: {avg_record_size:.1f}")

        # Calculate raw estimate based on total bytes and average record size
        raw_estimate = total_bytes / avg_record_size

    except Exception as e:
        logging.warning(f"Error sampling sequences: {e}")
        # typical 2000 bytes per sequence for ITS-length sequences
        raw_estimate = total_bytes / 2000

        logging.warning(f"Using fallback estimation method: estimated {raw_estimate:,} sequences")

    # Step 3: Apply adjustment factor and round to nice number
    if raw_estimate > 0:
        # Apply a slightly conservative adjustment factor
        raw_estimate *= 1.05

        # Round to a pleasing number (2 significant figures for large numbers)
        if raw_estimate > 10000:
            magnitude = math.floor(math.log10(raw_estimate))
            scaled = raw_estimate / (10 ** (magnitude - 1))
            rounded = math.ceil(scaled) * (10 ** (magnitude - 1))
            estimated_sequences = int(rounded)
        else:
            # For smaller numbers, just round up to nearest 100
            estimated_sequences = math.ceil(raw_estimate / 100) * 100
    else:
        estimated_sequences = 500000  # Conservative default

    logging.info(
        f"Estimated {estimated_sequences:,} sequences based on {'compressed' if is_compressed else 'uncompressed'} file")

    return estimated_sequences


def get_gzip_info(filename: str) -> Tuple[int, int]:
    """
    Get compressed and uncompressed file sizes using gzip -l

    Args:
        filename: Path to gzipped file

    Returns:
        Tuple[int, int]: (compressed_size, uncompressed_size)

    Raises:
        subprocess.CalledProcessError: If gzip -l fails
        ValueError: If output parsing fails
    """
    try:
        # Run gzip -l on the file
        result = subprocess.run(['gzip', '-l', filename],
                                capture_output=True,
                                text=True,
                                check=True)

        # Parse the output to extract compressed and uncompressed sizes
        output = result.stdout.strip()

        # Look for the line with sizes using regex
        # Format: compressed uncompressed ratio uncompressed_name
        match = re.search(r'(\d+)\s+(\d+)', output)
        if match:
            compressed_size = int(match.group(1))
            uncompressed_size = int(match.group(2))
            return compressed_size, uncompressed_size
        else:
            raise ValueError(f"Failed to parse gzip output: {output}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running gzip -l: {e}")
        raise
    except Exception as e:
        logging.error(f"Unexpected error analyzing compression: {e}")
        raise

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


    num_processes = args.threads if args.threads > 0 else multiprocessing.cpu_count()
    logging.info(f"Will run {num_processes} worker processes")

    # Create shared timestamp for consistent trace file naming
    start_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    with Pool(processes=num_processes,
              initializer=init_worker,
              initargs=(specimens, parameters.max_dist_index, args, start_timestamp)) as pool:

        worker_func = partial(worker, specimens=specimens, args=args)

        try:
            pbar = tqdm(total=total_seqs, desc="Processing sequences", unit="seq")
            
            # Track totals for match rate
            total_processed = 0
            total_matched = 0
            
            # Create work items with cumulative sequence indexing for trace IDs
            def create_work_items():
                cumulative_idx = 0
                for i, batch in enumerate(iter_batches(seq_records, sequence_block_size,
                                                       last_seq_to_output, all_seqs)):
                    work_item = SequenceBatch(i, batch, parameters, cumulative_idx)
                    cumulative_idx += len(batch)
                    yield work_item
            
            for batch_counts in pool.imap_unordered(worker_func, create_work_items()):
                if batch_counts:
                    batch_total, batch_matched = batch_counts
                    total_processed += batch_total
                    total_matched += batch_matched
                    pbar.update(batch_total)
                    
                    # Update progress with match rate
                    if total_processed > 0:
                        match_rate = total_matched / total_processed
                        pbar.set_description(f"Processing sequences [Match rate: {match_rate:.1%}]")

            pbar.close()
            
            # Log final statistics
            if total_processed > 0:
                final_match_rate = total_matched / total_processed
                logging.info(f"Processed {total_processed:,} sequences, match rate: {final_match_rate:.1%}")

        except WorkerException as e:
            logging.error(f"Unexpected error in worker (see details above): {e}")
            sys.exit(1)

    elapsed = timeit.default_timer() - start_time
    logging.info(f"Elapsed time: {elapsed:.2f} seconds")
    
    # Clean up empty directories after all processing is complete
    if args.output_to_files:
        cleanup_empty_directories(args.output_dir)

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


def write_all_primers_fasta(output_dir: str, fwd_primers: List[PrimerInfo], rev_primers: List[PrimerInfo]):
    """
    Write a primers.fasta file containing all forward and reverse primers to the specified directory.

    Args:
        output_dir: Directory to write the primers.fasta file to
        fwd_primers: List of forward primer info objects
        rev_primers: List of reverse primer info objects
    """
    primer_file_path = os.path.join(output_dir, "primers.fasta")

    with open(primer_file_path, 'w') as f:
        # Write all forward primers
        for fwd_primer in fwd_primers:
            f.write(f">{fwd_primer.name} position=forward pool={','.join(fwd_primer.pools)}\n")
            f.write(f"{fwd_primer.primer}\n")

        # Write all reverse primers
        for rev_primer in rev_primers:
            f.write(f">{rev_primer.name} position=reverse pool={','.join(rev_primer.pools)}\n")
            f.write(f"{rev_primer.primer}\n")


def create_output_files(args, specimens):
    """Create directory structure for output files and add primers.fasta files"""
    if args.output_to_files:
        # Create base output directory
        os.makedirs(args.output_dir, exist_ok=True)

        # Create match-type directories at the top level
        match_types = ["full", "partial", "unknown"]
        for match_type in match_types:
            os.makedirs(os.path.join(args.output_dir, match_type), exist_ok=True)
        
        # Create unknown pool directories for sequences with unrecognized pools or no primers
        for match_type in ["partial", "unknown"]:
            unknown_pool_dir = os.path.join(args.output_dir, match_type, "unknown")
            os.makedirs(unknown_pool_dir, exist_ok=True)
            # Create directories for various primer detection scenarios
            os.makedirs(os.path.join(unknown_pool_dir, "unknown-unknown"), exist_ok=True)
            
            # Also need to handle cases where primers are detected but not in recognized pools
            # These directories will be created dynamically based on detected primers
            # For now, we'll rely on the write_sequence method to create them as needed
        
        # Create pool directories under each match type
        for pool in specimens._primer_registry.get_pools():
            # Create pool directories under each match type
            for match_type in match_types:
                pool_dir = os.path.join(args.output_dir, match_type, pool)
                os.makedirs(pool_dir, exist_ok=True)
            
            # Write primers.fasta file to the pool-level full directory
            pool_full_path = os.path.join(args.output_dir, "full", pool)
            write_all_primers_fasta(pool_full_path, 
                                   specimens._primer_registry.get_pool_primers(pool, Primer.FWD),
                                   specimens._primer_registry.get_pool_primers(pool, Primer.REV))

            # Create primer pair directories under appropriate match types
            for fwd_primer in specimens._primer_registry.get_pool_primers(pool, Primer.FWD):
                for rev_primer in specimens._primer_registry.get_pool_primers(pool, Primer.REV):
                    primer_dir = f"{fwd_primer.name}-{rev_primer.name}"
                    
                    # Create primer-pair directories under each match type
                    for match_type in match_types:
                        primer_full_dir = os.path.join(args.output_dir, match_type, pool, primer_dir)
                        os.makedirs(primer_full_dir, exist_ok=True)
                    
                    # Write primers.fasta file for the full match primer-pair directory
                    full_primer_dir = os.path.join(args.output_dir, "full", pool, primer_dir)
                    write_primers_fasta(full_primer_dir, fwd_primer, rev_primer)

                # Create partial primer directories (forward primer only)
                unknown_primer_dir = f"{fwd_primer.name}-unknown"
                os.makedirs(os.path.join(args.output_dir, "partial", pool, unknown_primer_dir), exist_ok=True)
                os.makedirs(os.path.join(args.output_dir, "unknown", pool, unknown_primer_dir), exist_ok=True)

            # Create partial primer directories (reverse primer only)
            for rev_primer in specimens._primer_registry.get_pool_primers(pool, Primer.REV):
                unknown_primer_dir = f"unknown-{rev_primer.name}"
                os.makedirs(os.path.join(args.output_dir, "partial", pool, unknown_primer_dir), exist_ok=True)
                os.makedirs(os.path.join(args.output_dir, "unknown", pool, unknown_primer_dir), exist_ok=True)

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


    with OutputManager(args.output_dir, args.output_file_prefix, args.isfastq) as output_manager:
        num_seqs = 0
        prefilter = PassthroughPrefilter()
        if not args.disable_prefilter:
            barcode_rcs = barcodes_for_bloom_prefilter(specimens)
            cache_path = BloomPrefilter.get_cache_path(barcode_rcs, parameters.max_dist_index)
            prefilter = BloomPrefilter.load_readonly(cache_path, barcode_rcs, parameters.max_dist_index)

        pbar = tqdm(total=total_seqs, desc="Processing sequences", unit="seq")
        
        # Track totals for match rate
        total_processed = 0
        total_matched = 0
        
        while all_seqs or num_seqs < last_seq_to_output:
            to_read = sequence_block_size if all_seqs else min(sequence_block_size, last_seq_to_output - num_seqs)
            seq_batch = list(itertools.islice(seq_records, to_read))
            if not seq_batch:
                break

            # Create trace logger for single-threaded mode if needed
            trace_logger = None
            if args.diagnostics:
                trace_logger = TraceLogger(
                    enabled=True,
                    verbosity=args.diagnostics,
                    output_dir=args.output_dir,
                    worker_id="main",
                    start_timestamp=datetime.now().strftime("%Y%m%d_%H%M%S")
                )
                trace_logger.__enter__()

            write_ops, batch_total, batch_matched = process_sequences(
                seq_batch, parameters, specimens, args, prefilter, trace_logger, num_seqs)
            
            for write_op in write_ops:
                output_write_operation(write_op, output_manager, args, trace_logger)
            
            # Close trace logger after processing and writing batch
            if trace_logger:
                trace_logger.__exit__(None, None, None)

            num_seqs += batch_total
            total_processed += batch_total
            total_matched += batch_matched
            pbar.update(batch_total)
            
            # Update progress with match rate
            if total_processed > 0:
                match_rate = total_matched / total_processed
                pbar.set_description(f"Processing sequences [Match rate: {match_rate:.1%}]")

        pbar.close()
        
        # Log final statistics
        if total_processed > 0:
            final_match_rate = total_matched / total_processed
            logging.info(f"Processed {total_processed:,} sequences, match rate: {final_match_rate:.1%}")
    
    elapsed = timeit.default_timer() - start_time
    logging.info(f"Elapsed time: {elapsed:.2f} seconds")
    
    # Clean up empty directories after all processing is complete
    if args.output_to_files:
        cleanup_empty_directories(args.output_dir)


def cleanup_empty_directories(output_dir: str):
    """Remove empty directories from the output tree.
    
    Works bottom-up to remove directories that contain no files,
    only removing directories that were created by specimux.
    """
    if not os.path.exists(output_dir):
        return
    
    # List to track directories removed for logging
    removed_dirs = []
    
    # Walk the directory tree bottom-up
    for dirpath, dirnames, filenames in os.walk(output_dir, topdown=False):
        # Skip the base output directory itself
        if dirpath == output_dir:
            continue
            
        # Check if directory is empty (no files and no remaining subdirectories)
        try:
            if not filenames and not os.listdir(dirpath):
                os.rmdir(dirpath)
                removed_dirs.append(dirpath)
        except OSError:
            # Directory not empty or cannot be removed, skip it
            pass
    
    if removed_dirs:
        logging.debug(f"Removed {len(removed_dirs)} empty directories")
        for dir_path in removed_dirs[:10]:  # Show first 10 for debugging
            logging.debug(f"  Removed: {dir_path}")
        if len(removed_dirs) > 10:
            logging.debug(f"  ... and {len(removed_dirs) - 10} more")







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

    # Log multiple match strategy
    logging.info(f"Using multiple match strategy: {args.resolve_multiple_matches}")

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


