#!/usr/bin/env python3

"""
File I/O and output management utilities for specimux.

This module contains io utils functionality
extracted from the original specimux.py monolithic file.
"""

import argparse
import csv
import fcntl
import gzip
import hashlib
import io
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union, TYPE_CHECKING

if TYPE_CHECKING:
    from .trace import TraceLogger

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from cachetools import LRUCache
from tqdm import tqdm

from .constants import Primer, SampleId, TrimMode, ResolutionType
from .models import PrimerInfo, WriteOperation
from .databases import PrimerDatabase, Specimens


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
            # Full match (including dereplicated and multiple specimen matches)
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
            
            trace_logger.log_sequence_output(write_op.trace_sequence_id, write_op.sample_id,
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

def read_primers_file(filename: str) -> PrimerDatabase:
    """
    Read primers file and build primer registry

    Returns:
        PrimerRegistry object managing all primers and their relationships
    """
    registry = PrimerDatabase()

    for file_index, record in enumerate(SeqIO.parse(filename, "fasta")):
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
        primer = PrimerInfo(name, sequence, direction, pool_names, file_index=file_index)
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

        empty_barcode_errors = []
        for row_num, row in enumerate(reader, start=1):
            try:
                b1 = row['FwIndex'].upper()
                b2 = row['RvIndex'].upper()

                # Check for empty barcodes (single-indexed demultiplexing not supported)
                if not b1.strip() or not b2.strip():
                    empty_barcode_errors.append(
                        f"Row {row_num} ({row['SampleID']}): "
                        f"{'FwIndex is empty' if not b1.strip() else 'RvIndex is empty'}"
                    )
                    continue

                specimens.add_specimen(
                    specimen_id=row['SampleID'],
                    pool=row['PrimerPool'],
                    b1=b1,
                    p1=row['FwPrimer'],
                    b2=b2,
                    p2=row['RvPrimer']
                )
            except (KeyError, ValueError) as e:
                raise ValueError(f"Error processing row {row_num}: {e}")

        if empty_barcode_errors:
            raise ValueError(
                f"Empty barcodes found in {len(empty_barcode_errors)} specimen(s). "
                f"Single-indexed demultiplexing is not supported.\n"
                + "\n".join(empty_barcode_errors[:10])
                + (f"\n... and {len(empty_barcode_errors) - 10} more" if len(empty_barcode_errors) > 10 else "")
            )

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

def output_write_operation(write_op: WriteOperation,
                           output_manager: OutputManager,
                           args: argparse.Namespace,
                           trace_logger: Optional['TraceLogger'] = None) -> None:
    if not args.output_to_files:
        fh = sys.stdout
        formatted_seq = write_op.sequence
        if args.color:
            from .alignment import color_sequence
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

def cleanup_empty_directories(output_dir: str):
    """Remove empty directories from the output tree.
    
    Works bottom-up to remove directories that contain no files,
    only removing directories that were created by specimux.
    Ignores primer metadata files when determining emptiness.
    """
    if not os.path.exists(output_dir):
        return
    
    # Files that don't count as "real" content for emptiness check
    METADATA_FILES = {'primers.fasta', 'primers.txt'}
    
    # List to track directories removed for logging
    removed_dirs = []
    
    # Walk the directory tree bottom-up
    for dirpath, dirnames, filenames in os.walk(output_dir, topdown=False):
        # Skip the base output directory itself
        if dirpath == output_dir:
            continue
            
        # Check if directory has any non-metadata files
        content_files = [f for f in filenames if f not in METADATA_FILES]
        
        # Check if directory is effectively empty (no content files and no subdirectories)
        try:
            # Get list of subdirectories
            subdirs = [d for d in os.listdir(dirpath) 
                      if os.path.isdir(os.path.join(dirpath, d))]
            
            if not content_files and not subdirs:
                # Remove metadata files first, then the directory
                for f in filenames:
                    if f in METADATA_FILES:
                        os.remove(os.path.join(dirpath, f))
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
