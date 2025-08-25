#!/usr/bin/env python3

"""
Specimux: Demultiplex MinION sequences by dual barcode indexes and primers.

Copyright (c) 2024-2025 Josh Walker

Specimux is inspired by minibar.py (https://github.com/calacademy-research/minibar),
originally developed by the California Academy of Sciences. While conceptually influenced
by minibar's approach to demultiplexing, Specimux is an independent implementation with
substantial algorithmic enhancements, architectural improvements, and additional features.

This software is released under the BSD 2-Clause License.
For full terms of use and distribution, please see the LICENSE file accompanying this program.
"""

# Import all modules to maintain compatibility with existing code
# This file now serves as a compatibility layer after refactoring

# Constants and enums
from .constants import (
    IUPAC_EQUIV, IUPAC_CODES,
    AlignMode, SampleId, TrimMode, MultipleMatchStrategy,
    ResolutionType, Barcode, Primer, Orientation
)

# Data models
from .models import (
    PrimerInfo, AlignmentResult, CandidateMatch,
    MatchParameters, WriteOperation, SequenceBatch,
    WorkerException
)

# Database and registry classes
from .databases import (
    PrimerDatabase, Specimens,
    BarcodePrefilter, PassthroughPrefilter
)

# File I/O and output management
from .io_utils import (
    FileHandleCache, CachedFileManager, OutputManager,
    read_primers_file, read_specimen_file,
    detect_file_format, open_sequence_file,
    output_write_operation, get_gzip_info,
    cleanup_locks, cleanup_empty_directories
)

# Trace logging
from .trace import TraceLogger

# Bloom filter prefiltering
from .bloom_filter import BloomPrefilter, barcodes_for_bloom_prefilter

# Sequence alignment
from .alignment import (
    align_seq, get_quality_seq, color_sequence
)

# Core demultiplexing pipeline
from .demultiplex import (
    create_write_operation, process_sequences,
    select_best_matches, resolve_specimen,
    determine_orientation, get_pool_from_primers,
    find_candidate_matches, match_one_end
)

# Multiprocessing utilities
from .multiprocessing_utils import (
    init_worker, cleanup_worker, worker
)

# Main orchestration functions
from .orchestration import (
    estimate_sequence_count,
    specimux_mp, specimux,
    write_primers_fasta, write_all_primers_fasta,
    create_output_files,
    subsample_top_quality, iter_batches,
    setup_match_parameters
)

# Export the main functions
__all__ = [
    # Main entry points
    'specimux', 'specimux_mp',
    # Core classes
    'PrimerDatabase', 'Specimens', 'TraceLogger',
    # Key functions
    'read_primers_file', 'read_specimen_file',
    'setup_match_parameters', 'process_sequences',
]