#!/usr/bin/env python3

"""
Trace logging system for diagnostic output.

This module contains trace functionality
extracted from the original specimux.py monolithic file.
"""

import csv
import hashlib
import logging
import os
from collections import defaultdict
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Set

from .constants import Barcode, Primer, ResolutionType, Orientation
from .models import CandidateMatch


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

    def log_sequence_trim_empty(self, sequence_id: str, trim_mode: str,
                                 trim_start: int, trim_end: int, seq_length: int,
                                 p1_name: str, p2_name: str):
        """Log when trimming would produce an empty sequence."""
        self._log_event(sequence_id, 'SEQUENCE_TRIM_EMPTY', trim_mode,
                       trim_start, trim_end, seq_length, p1_name, p2_name)

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

    def log_dereplicate_expanded(self, sequence_id: str, match_count: int, expanded_count: int):
        """Log dereplication expansion step."""
        self._log_event(sequence_id, 'DEREPLICATE_EXPANDED', match_count, expanded_count)

    def log_dereplicate_selected(self, sequence_id: str, specimen_id: str,
                                 alternatives_count: int,
                                 scores: tuple):
        """Log dereplication selection for a specimen.

        Args:
            scores: Tuple of (barcode_dist, primer_dist, file_idx) for the winning match
        """
        barcode_dist, primer_dist, file_idx = scores
        self._log_event(sequence_id, 'DEREPLICATE_SELECTED', specimen_id,
                       alternatives_count, barcode_dist, primer_dist, file_idx)

    def log_dereplicate_partial_selected(self, sequence_id: str, direction: str,
                                         barcode: str, alternatives_count: int,
                                         scores: tuple):
        """Log partial match dereplication selection.

        Args:
            direction: 'forward' or 'reverse' (which barcode matched)
            barcode: The barcode sequence used for grouping
            scores: Tuple of (barcode_dist, neg_primer_count, primer_dist, file_idx)
        """
        barcode_dist, neg_primer_count, primer_dist, file_idx = scores
        primer_count = -neg_primer_count  # Convert back to positive
        self._log_event(sequence_id, 'DEREPLICATE_PARTIAL_SELECTED', direction,
                       barcode, alternatives_count, barcode_dist, primer_count,
                       primer_dist, file_idx)

    def log_dereplicate_unknown_selected(self, sequence_id: str,
                                          alternatives_count: int,
                                          primer_count: int,
                                          primer_dist: int,
                                          file_idx: int):
        """Log unknown match dereplication selection.

        Args:
            sequence_id: Sequence identifier
            alternatives_count: Number of alternative matches considered
            primer_count: Number of primers matched in winning match (0, 1, or 2)
            primer_dist: Combined primer edit distance
            file_idx: File order sum
        """
        self._log_event(sequence_id, 'DEREPLICATE_UNKNOWN_SELECTED',
                       alternatives_count, primer_count, primer_dist, file_idx)

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

