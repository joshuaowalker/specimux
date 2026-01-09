#!/usr/bin/env python3

"""
Data models and structures for specimux.

This module contains all data model classes, dataclasses, and named tuples
used throughout the specimux package.
"""

from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional, Set
from typing import NamedTuple

from Bio.Seq import reverse_complement
from operator import itemgetter

from .constants import Primer, Barcode, ResolutionType


class PrimerInfo:
    """Information about a primer sequence and its associations."""

    def __init__(self, name: str, seq: str, direction: Primer, pools: List[str], file_index: int = 0):
        self.name = name
        self.primer = seq.upper()
        self.direction = direction
        self.primer_rc = reverse_complement(seq.upper())
        self.barcodes = set()
        self.specimens = set()
        self.pools = pools
        self.file_index = file_index  # Order in primers.fasta file (0-indexed)


class AlignmentResult:
    """Wrapper around edlib match result."""
    
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
        """Return a new MatchResult with locations relative to the reversed sequence."""
        m = AlignmentResult(self._edlib_match.copy())
        m._reverse(seq_length)
        return m

    def _reverse(self, seq_length):
        """Update locations to be relative to the start of the reversed sequence."""
        if self.distance() == -1: return

        r = self._edlib_match
        r['locations'] = [(seq_length - loc[1] - 1, seq_length - loc[0] - 1) for loc in r['locations']]

    def adjust_start(self, s):
        """Update start in locations by s."""
        if self.distance() == -1: return

        self._edlib_match['locations'] = [(loc[0] + s, loc[1] + s) for loc in self._edlib_match['locations']]  # adjust relative match to absolute


class CandidateMatch:
    """Represents a potential match between a sequence and primers/barcodes."""
    
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
        """Set the primer pool for this match."""
        self._pool = pool

    def get_pool(self) -> Optional[str]:
        """Get the primer pool for this match."""
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

    def has_barcode_match(self, which) -> bool:
        return (len(self._b1_matches) > 0) if which is Barcode.B1 else (len(self._b2_matches) > 0)

    def get_barcode_distances(self, which) -> Dict[str, float]:
        m = self._b1_matches if which is Barcode.B1 else self._b2_matches
        return dict([(bc, d) for bc, _, d in m])

    def get_best_barcodes(self, which) -> List[str]:
        return self.best_b1() if which is Barcode.B1 else self.best_b2()

    def primer_distance(self, which) -> int:
        if which is Primer.FWD:
            return self.p1_match.distance() if self.p1_match else -1
        else:
            return self.p2_match.distance() if self.p2_match else -1

    def get_primer(self, which: Primer) -> Optional[PrimerInfo]:
        if which is Primer.FWD:
            return self._p1
        else:
            return self._p2

    def has_both_primers(self) -> bool:
        return self.p1_match is not None and self.p2_match is not None

    def has_either_primer(self) -> bool:
        return self.p1_match is not None or self.p2_match is not None

    def has_both_barcodes(self) -> bool:
        return self.has_barcode_match(Barcode.B1) and self.has_barcode_match(Barcode.B2)

    def has_full_match(self) -> bool:
        return self.has_both_primers() and self.has_both_barcodes()

    def b1_span(self):
        if not self.has_barcode_match(Barcode.B1):
            return None
        return self._b1_matches[0][1].location()

    def b2_span(self):
        if not self.has_barcode_match(Barcode.B2):
            return None
        return self._b2_matches[0][1].location()

    def p1_span(self):
        if not self.p1_match:
            return None
        return self.p1_match.location()

    def p2_span(self):
        if not self.p2_match:
            return None
        return self.p2_match.location()

    def total_distance(self) -> int:
        """Calculate the total edit distance across all matched elements."""
        distance = 0
        if self.p1_match:
            distance += self.p1_match.distance()
        if self.has_barcode_match(Barcode.B1):
            distance += self.b1_distance()
        if self.has_barcode_match(Barcode.B2):
            distance += self.b2_distance()
        if self.p2_match:
            distance += self.p2_match.distance()
        return distance

    def distance_code(self) -> str:
        """Generate a distance code string summarizing match quality."""
        p1d = self.primer_distance(Primer.FWD)
        b1d = self.b1_distance()
        b2d = self.b2_distance()
        p2d = self.primer_distance(Primer.REV)
        
        p1_str = str(p1d) if p1d >= 0 else 'X'
        b1_str = str(b1d) if b1d >= 0 else 'X'
        b2_str = str(b2d) if b2d >= 0 else 'X'
        p2_str = str(p2d) if p2d >= 0 else 'X'
        
        return f"{p1_str},{b1_str},{b2_str},{p2_str}"

    def num_matches(self) -> int:
        """Count how many components (primers/barcodes) matched."""
        count = 0
        if self.p1_match:
            count += 1
        if self.has_barcode_match(Barcode.B1):
            count += 1
        if self.has_barcode_match(Barcode.B2):
            count += 1
        if self.p2_match:
            count += 1
        return count

    def get_barcode_loc(self, which):
        if which is Barcode.B1:
            return self.b1_span() if self.has_barcode_match(Barcode.B1) else None
        else:
            return self.b2_span() if self.has_barcode_match(Barcode.B2) else None

    def get_primer_loc(self, which):
        if which is Primer.FWD:
            return self.p1_span() if self.p1_match else None
        else:
            return self.p2_span() if self.p2_match else None
    
    # Convenience methods for backward compatibility
    def has_b1_match(self) -> bool:
        """Check if forward barcode matched (convenience method)."""
        return self.has_barcode_match(Barcode.B1)
    
    def has_b2_match(self) -> bool:
        """Check if reverse barcode matched (convenience method)."""
        return self.has_barcode_match(Barcode.B2)
    
    def get_p1(self) -> PrimerInfo:
        """Get forward primer (convenience method)."""
        return self._p1
    
    def get_p2(self) -> PrimerInfo:
        """Get reverse primer (convenience method)."""
        return self._p2
    
    def get_p1_location(self):
        """Get forward primer location (convenience method)."""
        return self.p1_match.location() if self.p1_match else None
    
    def get_p2_location(self):
        """Get reverse primer location (convenience method)."""
        return self.p2_match.location() if self.p2_match else None
    
    def get_barcode1_location(self):
        """Get forward barcode location (convenience method)."""
        return self._b1_matches[0][1].location() if len(self._b1_matches) > 0 else None
    
    def get_barcode2_location(self):
        """Get reverse barcode location (convenience method)."""
        return self._b2_matches[0][1].location() if len(self._b2_matches) > 0 else None
    
    def interprimer_extent(self):
        """Get the extent between primers."""
        s = 0
        e = self.sequence_length
        if self.p1_match:
            s = self.p1_match.location()[1] + 1
        if self.p2_match:
            e = self.p2_match.location()[0]
        
        return (s, e)
    
    def interbarcode_extent(self):
        """Get the extent between barcodes."""
        s = 0
        e = self.sequence_length
        if self.p1_match:
            s = self.p1_match.location()[0]
        if self.p2_match:
            e = self.p2_match.location()[1] + 1
        
        return (s, e)
    
    def intertail_extent(self):
        """Get the extent between tails (outermost matches)."""
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
        """Adjust all location coordinates by start offset."""
        for b in self._b1_matches: b[1].adjust_start(-1 * start)
        for b in self._b2_matches: b[1].adjust_start(-1 * start)
        if self.p1_match:
            self.p1_match.adjust_start(-1 * start)
        if self.p2_match:
            self.p2_match.adjust_start(-1 * start)


class MatchParameters:
    """Parameters for matching sequences against primers and barcodes."""
    
    def __init__(self, max_dist_primers: Dict[str, int], max_dist_index: int, search_len: int, preorient: bool):
        self.max_dist_primers = max_dist_primers
        self.max_dist_index = max_dist_index
        self.search_len = search_len
        self.preorient = preorient


class WriteOperation(NamedTuple):
    """Information needed to write a demultiplexed sequence to output."""
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
    trace_sequence_id: Optional[str] = None


class SequenceBatch(NamedTuple):
    """A batch of sequences to process together."""
    seq_number: int
    seq_records: List  # This now contains actual sequence records, not an iterator
    parameters: MatchParameters
    start_idx: int     # Starting index for sequence ID generation


class WorkerException(Exception):
    """Exception raised in worker processes."""
    pass