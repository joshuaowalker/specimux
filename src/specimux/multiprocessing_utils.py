#!/usr/bin/env python3

"""
Multiprocessing utilities for parallel sequence processing.

This module contains multiprocessing utils functionality
extracted from the original specimux.py monolithic file.
"""

import argparse
import atexit
import logging
import multiprocessing
import os
import sys
import traceback
from datetime import datetime
from typing import List, Optional, Tuple

from .constants import TrimMode
from .databases import Specimens, PassthroughPrefilter
from .models import MatchParameters, SequenceBatch, WorkerException
from .trace import TraceLogger
from .bloom_filter import BloomPrefilter, barcodes_for_bloom_prefilter
from .io_utils import OutputManager, output_write_operation
from .demultiplex import process_sequences


def init_worker(specimens: Specimens, max_distance: int, args: argparse.Namespace, start_timestamp: str = None):
    """Initialize worker process with shared resources"""
    global _output_manager, _barcode_prefilter, _trace_logger
    try:
        from .cli import setup_logging
        setup_logging(args.debug, args.output_dir if args.output_to_files else None, is_worker=True)

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

# Global variables for worker processes
_barcode_prefilter = None
_output_manager = None
_trace_logger = None

