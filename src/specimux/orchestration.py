#!/usr/bin/env python3

"""
Main orchestration functions for specimux pipeline.

This module contains orchestration functionality
extracted from the original specimux.py monolithic file.
"""

import argparse
import glob
import itertools
import logging
import math
import multiprocessing
import os
import shutil
import subprocess
import sys
import tempfile
import timeit
from collections import Counter, defaultdict
from datetime import datetime
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from Bio import SeqIO
from Bio.Seq import reverse_complement
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from .constants import Primer, TrimMode, MultipleMatchStrategy
from .databases import BarcodePrefilter, PassthroughPrefilter, PrimerDatabase, Specimens
from .models import MatchParameters, PrimerInfo, SequenceBatch
from .io_utils import (
    CachedFileManager, OutputManager, cleanup_empty_directories, cleanup_locks,
    detect_file_format, open_sequence_file, read_primers_file, read_specimen_file
)
from .trace import TraceLogger
from .bloom_filter import BloomPrefilter, barcodes_for_bloom_prefilter
from .multiprocessing_utils import init_worker, worker, WorkerException
from .demultiplex import process_sequences
from .io_utils import get_gzip_info, output_write_operation
import edlib


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
    
    # Create subsamples if requested
    if args.output_to_files and args.sample_topq > 0:
        subsample_top_quality(args.output_dir, args.sample_topq)

    cleanup_locks(args.output_dir)


def write_primers_fasta(output_dir: str, fwd_primer: PrimerInfo, rev_primer: PrimerInfo):
    """
    Write primers.fasta and primers.txt files containing the forward and reverse primers to the specified directory.
    
    The primers.fasta file includes full headers with position and pool information.
    The primers.txt file includes simplified headers with only the primer ID for compatibility with older tools.

    Args:
        output_dir: Directory to write the primers files to
        fwd_primer: Forward primer info
        rev_primer: Reverse primer info
    """
    # Write primers.fasta with full headers
    primer_fasta_path = os.path.join(output_dir, "primers.fasta")
    with open(primer_fasta_path, 'w') as f:
        # Write forward primer
        f.write(f">{fwd_primer.name} position=forward pool={','.join(fwd_primer.pools)}\n")
        f.write(f"{fwd_primer.primer}\n")

        # Write reverse primer
        f.write(f">{rev_primer.name} position=reverse pool={','.join(rev_primer.pools)}\n")
        f.write(f"{rev_primer.primer}\n")
    
    # Write primers.txt with simplified headers
    primer_txt_path = os.path.join(output_dir, "primers.txt")
    with open(primer_txt_path, 'w') as f:
        # Write forward primer with simplified header
        f.write(f">{fwd_primer.name}\n")
        f.write(f"{fwd_primer.primer}\n")

        # Write reverse primer with simplified header
        f.write(f">{rev_primer.name}\n")
        f.write(f"{rev_primer.primer}\n")


def write_all_primers_fasta(output_dir: str, fwd_primers: List[PrimerInfo], rev_primers: List[PrimerInfo]):
    """
    Write primers.fasta and primers.txt files containing all forward and reverse primers to the specified directory.
    
    The primers.fasta file includes full headers with position and pool information.
    The primers.txt file includes simplified headers with only the primer ID for compatibility with older tools.

    Args:
        output_dir: Directory to write the primers files to
        fwd_primers: List of forward primer info objects
        rev_primers: List of reverse primer info objects
    """
    # Write primers.fasta with full headers
    primer_fasta_path = os.path.join(output_dir, "primers.fasta")
    with open(primer_fasta_path, 'w') as f:
        # Write all forward primers
        for fwd_primer in fwd_primers:
            f.write(f">{fwd_primer.name} position=forward pool={','.join(fwd_primer.pools)}\n")
            f.write(f"{fwd_primer.primer}\n")

        # Write all reverse primers
        for rev_primer in rev_primers:
            f.write(f">{rev_primer.name} position=reverse pool={','.join(rev_primer.pools)}\n")
            f.write(f"{rev_primer.primer}\n")
    
    # Write primers.txt with simplified headers
    primer_txt_path = os.path.join(output_dir, "primers.txt")
    with open(primer_txt_path, 'w') as f:
        # Write all forward primers with simplified headers
        for fwd_primer in fwd_primers:
            f.write(f">{fwd_primer.name}\n")
            f.write(f"{fwd_primer.primer}\n")

        # Write all reverse primers with simplified headers
        for rev_primer in rev_primers:
            f.write(f">{rev_primer.name}\n")
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

def subsample_top_quality(output_dir: str, top_n: int):
    """
    Create subsample directories with top N sequences by average quality score.
    
    Processes all .fastq files in the full/ subdirectory and creates corresponding
    files in subsample/ with only the top N sequences sorted by average quality.
    
    Args:
        output_dir: Base output directory containing full/ subdirectory
        top_n: Number of top-quality sequences to retain in each file
    """
    from Bio import SeqIO
    
    full_dir = os.path.join(output_dir, "full")
    if not os.path.exists(full_dir):
        logging.warning(f"Full directory not found: {full_dir}")
        return
    
    logging.info(f"Creating subsamples with top {top_n} sequences by quality...")
    
    # Walk through all directories in full/
    processed_files = 0
    failed_files = 0
    
    for root, dirs, files in os.walk(full_dir):
        for fname in files:
            if fname.endswith(".fastq"):
                # Create corresponding subsample directory structure
                rel_path = os.path.relpath(root, full_dir)
                subsample_dir = os.path.join(output_dir, "subsample", rel_path)
                os.makedirs(subsample_dir, exist_ok=True)
                
                # Copy primer files to subsample directory
                for primer_file in ["primers.fasta", "primers.txt"]:
                    src_primer = os.path.join(root, primer_file)
                    if os.path.exists(src_primer):
                        import shutil
                        dst_primer = os.path.join(subsample_dir, primer_file)
                        shutil.copy2(src_primer, dst_primer)
                
                # Process the fastq file
                input_path = os.path.join(root, fname)
                output_path = os.path.join(subsample_dir, fname)
                
                try:
                    # Read all sequences and calculate average quality
                    records = list(SeqIO.parse(input_path, "fastq"))
                    
                    if not records:
                        continue
                    
                    # Calculate average quality for each record
                    def avg_quality(record):
                        if "phred_quality" in record.letter_annotations and len(record) > 0:
                            return sum(record.letter_annotations["phred_quality"]) / len(record)
                        return 0
                    
                    # Sort by average quality (highest first) and take top N
                    records.sort(key=avg_quality, reverse=True)
                    top_records = records[:top_n]
                    
                    # Write the subsampled sequences
                    SeqIO.write(top_records, output_path, "fastq")
                    
                    processed_files += 1
                    logging.debug(f"Subsampled {input_path}: {len(records)} → {len(top_records)} sequences")
                    
                except Exception as e:
                    logging.warning(f"Failed to subsample {input_path}: {e}")
                    failed_files += 1
    
    logging.info(f"Subsampling complete: {processed_files} files processed, {failed_files} failed")

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

    # Log dereplication strategy
    logging.info(f"Using dereplication strategy: {args.dereplicate}")

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

