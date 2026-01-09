#!/usr/bin/env python3

"""
Constants and enumerations for specimux.

This module contains all constant values and enumeration classes used throughout
the specimux package. No dependencies on other specimux modules.
"""

from enum import Enum

# IUPAC nucleotide codes and their equivalents
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
    """Alignment mode constants for edlib."""
    GLOBAL = 'NW'
    INFIX = 'HW'
    PREFIX = 'SHW'


class SampleId:
    """Sample ID prefixes and constants."""
    UNKNOWN = "unknown"
    PREFIX_FWD_MATCH = "barcode_fwd_"
    PREFIX_REV_MATCH = "barcode_rev_"


class TrimMode:
    """Trimming mode options."""
    PRIMERS = "primers"
    BARCODES = "barcodes"
    TAILS = "tails"
    NONE = "none"


class MultipleMatchStrategy:
    """Strategy for handling multiple equivalent matches."""
    NONE = "none"  # Default: output all equivalent matches
    BEST = "best"  # Group by specimen/barcode, select best match per group


class ResolutionType(Enum):
    """Types of specimen resolution outcomes."""
    FULL_MATCH = 1  # Both primers + both barcodes matched to a specimen
    PARTIAL_FORWARD = 2  # Forward barcode only (generates FWD_ONLY_ sample ID)
    PARTIAL_REVERSE = 3  # Reverse barcode only (generates REV_ONLY_ sample ID)
    MULTIPLE_SPECIMENS = 4  # Multiple specimen matches (shouldn't happen in new flow)
    UNKNOWN = 5  # No resolution possible
    DEREPLICATED_FULL = 6  # Full match selected via dereplication (best match per specimen)

    def to_string(self) -> str:
        """Convert resolution type to lowercase string for trace logging."""
        if self == ResolutionType.FULL_MATCH:
            return 'full_match'
        elif self == ResolutionType.PARTIAL_FORWARD:
            return 'partial_forward'
        elif self == ResolutionType.PARTIAL_REVERSE:
            return 'partial_reverse'
        elif self == ResolutionType.MULTIPLE_SPECIMENS:
            return 'multiple_specimens'
        elif self == ResolutionType.DEREPLICATED_FULL:
            return 'dereplicated_full'
        else:
            return 'unknown'

    def is_full_match(self) -> bool:
        """Check if this represents a successful full match."""
        return self in [ResolutionType.FULL_MATCH, ResolutionType.DEREPLICATED_FULL]

    def is_partial_match(self) -> bool:
        """Check if this represents a partial match."""
        return self in [ResolutionType.PARTIAL_FORWARD, ResolutionType.PARTIAL_REVERSE]

    def is_unknown(self) -> bool:
        """Check if this represents an unknown/failed resolution."""
        return self == ResolutionType.UNKNOWN


class Barcode(Enum):
    """Barcode position identifiers."""
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
    """Primer direction identifiers."""
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


class Orientation(Enum):
    """Sequence orientation states."""
    FORWARD = 1
    REVERSE = 2
    UNKNOWN = 3
    
    def to_string(self) -> str:
        """Convert orientation to lowercase string for trace logging."""
        return self.name.lower()