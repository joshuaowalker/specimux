#!/usr/bin/env python3
"""
Trace-based statistics aggregator for Specimux.

This tool parses trace event files and builds CandidateMatch objects for flexible
statistical analysis. Supports hierarchical text output and JSON export for
downstream visualization tools.

Usage:
    python trace_to_stats.py trace_dir/ --hierarchical pool primer_pair outcome
    python trace_to_stats.py trace_dir/ --sankey-data pool primer_pair outcome --output sankey.json
    python trace_to_stats.py trace_dir/ --count-by sequences --hierarchical orientation outcome
"""

import argparse
import csv
import json
import sys
from collections import defaultdict, Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Tuple
import glob
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class CandidateMatch:
    """Represents a single candidate match attempt or synthesized unmatched sequence.
    
    Design Philosophy:
    - Each CandidateMatch represents one potential primer-pair match location
    - Synthesized matches represent sequences with no primer matches found
    - All dimensions use "none"/"unknown" for missing/unresolved values
    - Computed properties provide derived dimensions for flexible analysis
    - Discard reasons enable detailed contamination and quality analysis
    """
    
    # Primary Keys
    candidate_match_id: Optional[str]  # None for unmatched sequences
    sequence_id: str
    
    # Core Dimensions
    orientation: str = "unknown"               # forward, reverse, unknown
    forward_primer: str = "none"               # primer name or "none"
    reverse_primer: str = "none"               # primer name or "none"
    pool: str = "none"                         # pool name or "none"
    forward_barcode: str = "none"              # barcode name or "none"
    reverse_barcode: str = "none"              # barcode name or "none"
    
    # Classification Dimensions
    match_type: str = "none"                   # both, forward_only, reverse_only, none
    resolution_type: str = "unknown"           # full_match, partial_forward, partial_reverse, unknown
    outcome: str = "unknown"                   # matched, partial, unknown, discarded
    selection_strategy: str = "none"           # unique, first, discarded, none
    discard_reason: str = "none"               # lower_score, none
    
    # Boolean Dimensions (computed from above)
    @property
    def forward_barcode_matched(self) -> bool:
        return self.forward_barcode != "none"
    
    @property
    def reverse_barcode_matched(self) -> bool:
        return self.reverse_barcode != "none"
    
    @property
    def specimen_matched(self) -> bool:
        return self.outcome == "matched"
    
    # Count Dimensions (computed)
    @property
    def barcode_count(self) -> int:
        count = 0
        if self.forward_barcode_matched:
            count += 1
        if self.reverse_barcode_matched:
            count += 1
        return count
    
    # Compound Dimensions (computed)
    @property
    def primer_pair(self) -> str:
        return f"{self.forward_primer}-{self.reverse_primer}"
    
    @property
    def barcode_pair(self) -> str:
        return f"{self.forward_barcode}-{self.reverse_barcode}"
    
    @property
    def outcome_detailed(self) -> str:
        """Detailed outcome including discard reason for analysis."""
        if self.outcome == "discarded" and self.discard_reason != "none":
            return f"discarded_{self.discard_reason}"
        else:
            return self.outcome


class TraceEventParser:
    """Parses trace events and builds CandidateMatch objects.
    
    Design Philosophy:
    - Event-driven reconstruction: Build analysis objects from trace events
    - Robust parsing: Handle incomplete/malformed events gracefully  
    - Synthesized completeness: Create matches for unmatched sequences
    - Single-pass efficiency: Parse all events once, build all matches
    
    Key Behaviors:
    - Groups events by sequence_id for coherent processing
    - Creates one CandidateMatch per PRIMER_MATCHED event (candidate location)
    - Synthesizes CandidateMatch for sequences with no primer matches
    - Enriches matches with barcode, selection, and outcome data from related events
    """
    
    def __init__(self):
        self.candidate_matches: List[CandidateMatch] = []
        self.sequence_events: Dict[str, List[Dict]] = defaultdict(list)
        self.event_stats = Counter()
    
    def parse_trace_files(self, trace_directory: str) -> List[CandidateMatch]:
        """Parse all trace files in directory and return CandidateMatch objects."""
        trace_files = glob.glob(str(Path(trace_directory) / "specimux_trace_*.tsv"))
        
        if not trace_files:
            raise ValueError(f"No trace files found in {trace_directory}")
        
        logger.info(f"Found {len(trace_files)} trace files")
        
        # Read all events
        total_events = 0
        for trace_file in trace_files:
            events_read = self._read_trace_file(trace_file)
            total_events += events_read
            logger.info(f"Read {events_read} events from {Path(trace_file).name}")
        
        logger.info(f"Read {total_events} total events for {len(self.sequence_events)} sequences")
        
        # Process events into CandidateMatch objects
        self._process_events()
        
        logger.info(f"Generated {len(self.candidate_matches)} candidate matches")
        logger.info(f"Event type distribution: {dict(self.event_stats)}")
        
        return self.candidate_matches
    
    def _read_trace_file(self, filepath: str) -> int:
        """Read a single trace file and collect events."""
        events_read = 0
        
        with open(filepath, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            try:
                header = next(reader)  # Skip header
            except StopIteration:
                logger.warning(f"Empty trace file: {filepath}")
                return 0
            
            for row in reader:
                if len(row) < 5:  # Skip malformed rows
                    continue
                
                # Parse base event
                event = {
                    'timestamp': row[0],
                    'worker_id': row[1],
                    'event_seq': int(row[2]) if row[2].isdigit() else 0,
                    'sequence_id': row[3],
                    'event_type': row[4]
                }
                
                # Parse event-specific fields
                if event['event_type'] == 'ORIENTATION_DETECTED' and len(row) >= 9:
                    event.update({
                        'orientation': row[5],
                        'forward_score': row[6],
                        'reverse_score': row[7],
                        'confidence': row[8]
                    })
                elif event['event_type'] == 'PRIMER_MATCHED' and len(row) >= 13:
                    event.update({
                        'candidate_match_id': row[5],
                        'match_type': row[6],
                        'forward_primer': row[7],
                        'reverse_primer': row[8],
                        'forward_distance': row[9],
                        'reverse_distance': row[10],
                        'pool': row[11],
                        'orientation_used': row[12]
                    })
                elif event['event_type'] == 'BARCODE_MATCHED' and len(row) >= 13:
                    event.update({
                        'candidate_match_id': row[5],
                        'match_type': row[6],
                        'forward_barcode': row[7],
                        'reverse_barcode': row[8],
                        'forward_distance': row[9],
                        'reverse_distance': row[10],
                        'forward_primer': row[11],
                        'reverse_primer': row[12]
                    })
                elif event['event_type'] == 'MATCH_SELECTED' and len(row) >= 12:
                    event.update({
                        'selection_strategy': row[5],
                        'forward_primer': row[6],
                        'reverse_primer': row[7],
                        'forward_barcode': row[8],
                        'reverse_barcode': row[9],
                        'pool': row[10],
                        'is_unique': row[11]
                    })
                elif event['event_type'] == 'SPECIMEN_RESOLVED' and len(row) >= 12:
                    event.update({
                        'specimen_id': row[5],
                        'resolution_type': row[6],
                        'pool': row[7],
                        'forward_primer': row[8],
                        'reverse_primer': row[9],
                        'forward_barcode': row[10],
                        'reverse_barcode': row[11]
                    })
                elif event['event_type'] == 'MATCH_DISCARDED' and len(row) >= 12:
                    event.update({
                        'candidate_match_id': row[5],
                        'forward_primer': row[6],
                        'reverse_primer': row[7],
                        'forward_barcode': row[8],
                        'reverse_barcode': row[9],
                        'score': row[10],
                        'discard_reason': row[11]
                    })
                elif event['event_type'] == 'NO_MATCH_FOUND' and len(row) >= 7:
                    event.update({
                        'stage_failed': row[5],
                        'reason': row[6]
                    })
                
                # Add to sequence events
                sequence_id = event['sequence_id']
                self.sequence_events[sequence_id].append(event)
                self.event_stats[event['event_type']] += 1
                events_read += 1
        
        return events_read
    
    def _process_events(self):
        """Process collected events into CandidateMatch objects."""
        for sequence_id, events in self.sequence_events.items():
            matches = self._process_sequence_events(sequence_id, events)
            self.candidate_matches.extend(matches)
    
    def _process_sequence_events(self, sequence_id: str, events: List[Dict]) -> List[CandidateMatch]:
        """Process events for a single sequence into CandidateMatch objects."""
        # Get orientation
        orientation = self._extract_orientation(events)
        
        # Get all candidate matches from PRIMER_MATCHED events
        primer_events = [e for e in events if e['event_type'] == 'PRIMER_MATCHED']
        
        if not primer_events:
            # No candidate matches - synthesize one unmatched sequence
            return [self._create_unmatched_sequence(sequence_id, orientation, events)]
        
        # Create CandidateMatch for each primer match
        candidate_matches = []
        for primer_event in primer_events:
            match = self._create_candidate_match(sequence_id, orientation, primer_event, events)
            candidate_matches.append(match)
        
        return candidate_matches
    
    def _extract_orientation(self, events: List[Dict]) -> str:
        """Extract orientation from events."""
        for event in events:
            if event['event_type'] == 'ORIENTATION_DETECTED':
                return event.get('orientation', 'unknown')
        return 'unknown'
    
    def _create_unmatched_sequence(self, sequence_id: str, orientation: str, events: List[Dict]) -> CandidateMatch:
        """Create CandidateMatch for sequence with no primer matches."""
        return CandidateMatch(
            candidate_match_id=None,
            sequence_id=sequence_id,
            orientation=orientation,
            outcome="unknown",
            selection_strategy="none"
        )
    
    def _create_candidate_match(self, sequence_id: str, orientation: str, 
                               primer_event: Dict, all_events: List[Dict]) -> CandidateMatch:
        """Create CandidateMatch from primer match and related events."""
        candidate_match_id = primer_event['candidate_match_id']
        
        # Start with primer match data
        match = CandidateMatch(
            candidate_match_id=candidate_match_id,
            sequence_id=sequence_id,
            orientation=orientation,
            forward_primer=primer_event.get('forward_primer', 'none'),
            reverse_primer=primer_event.get('reverse_primer', 'none'),
            pool=primer_event.get('pool', 'none'),
            match_type=primer_event.get('match_type', 'none')
        )
        
        # Enrich with barcode data
        self._enrich_with_barcode_data(match, candidate_match_id, all_events)
        
        # Determine final outcome
        self._determine_outcome(match, candidate_match_id, all_events)
        
        return match
    
    def _enrich_with_barcode_data(self, match: CandidateMatch, candidate_match_id: str, events: List[Dict]):
        """Add barcode information to CandidateMatch."""
        for event in events:
            if (event['event_type'] == 'BARCODE_MATCHED' and 
                event.get('candidate_match_id') == candidate_match_id):
                match.forward_barcode = event.get('forward_barcode', 'none')
                match.reverse_barcode = event.get('reverse_barcode', 'none')
                # Update match_type from barcode event if available (more accurate)
                if 'match_type' in event:
                    match.match_type = event['match_type']
                break
    
    def _determine_outcome(self, match: CandidateMatch, candidate_match_id: str, events: List[Dict]):
        """Determine final outcome and selection strategy."""
        # Check if this candidate was discarded
        for event in events:
            if (event['event_type'] == 'MATCH_DISCARDED' and 
                event.get('candidate_match_id') == candidate_match_id):
                match.outcome = "discarded"
                match.selection_strategy = "discarded"
                match.discard_reason = event.get('discard_reason', 'unknown')
                return
        
        # Check for match selection (this candidate was chosen)
        for event in events:
            if event['event_type'] == 'MATCH_SELECTED':
                # Determine if this specific candidate was selected by comparing details
                if (match.forward_primer == event.get('forward_primer', 'none') and
                    match.reverse_primer == event.get('reverse_primer', 'none') and
                    match.forward_barcode == event.get('forward_barcode', 'none') and
                    match.reverse_barcode == event.get('reverse_barcode', 'none')):
                    match.selection_strategy = event.get('selection_strategy', 'unknown')
                    break
        
        # Get final resolution from SPECIMEN_RESOLVED
        for event in events:
            if event['event_type'] == 'SPECIMEN_RESOLVED':
                match.resolution_type = event.get('resolution_type', 'unknown')
                
                # Map resolution type to outcome
                resolution = event.get('resolution_type', 'unknown')
                if resolution == 'full_match':
                    match.outcome = "matched"
                elif resolution in ['partial_forward', 'partial_reverse']:
                    match.outcome = "partial"
                else:
                    match.outcome = "unknown"
                break


class StatsAggregator:
    """Aggregates CandidateMatch objects into hierarchical statistics.
    
    Design Philosophy:
    - Dimension-agnostic: Works with any combination of CandidateMatch dimensions
    - Flexible counting: Count by candidate matches (locations) or unique sequences
    - Generic algorithms: No hardcoded dimension logic, fully general approach
    - Clean separation: Statistics logic independent of visualization concerns
    
    Key Features:
    - Hierarchical aggregation: Build nested trees of any dimension combination
    - Sankey flow generation: Create node/link data for any dimension sequence
    - Validation: Ensure requested dimensions exist and counting mode is valid
    - Performance: Efficient single-pass aggregation algorithms
    """
    
    def __init__(self, matches: List[CandidateMatch]):
        self.matches = matches
        self.available_dimensions = self._get_available_dimensions()
    
    def _get_available_dimensions(self) -> List[str]:
        """Get all available dimensions from CandidateMatch."""
        if not self.matches:
            return []
        
        # Get all attributes that are dimensions (not private/methods)
        sample_match = self.matches[0]
        dimensions = []
        
        for attr in dir(sample_match):
            if not attr.startswith('_') and not callable(getattr(sample_match, attr)):
                dimensions.append(attr)
        
        return sorted(dimensions)
    
    def get_hierarchical_stats(self, dimensions: List[str], count_by: str = "candidate_matches") -> Dict[str, Any]:
        """Build hierarchical statistics tree."""
        self._validate_dimensions(dimensions)
        self._validate_count_by(count_by)
        
        # Get the data to aggregate
        if count_by == "sequences":
            # Group by sequence_id and take first match per sequence
            sequence_matches = {}
            for match in self.matches:
                if match.sequence_id not in sequence_matches:
                    sequence_matches[match.sequence_id] = match
            data = list(sequence_matches.values())
        else:  # candidate_matches
            data = self.matches
        
        # Build hierarchical tree
        result = self._build_tree(data, dimensions, 0)
        
        return {
            "dimensions": dimensions,
            "count_by": count_by,
            "total_count": len(data),
            "data": result
        }
    
    def get_sankey_data(self, dimensions: List[str], count_by: str = "candidate_matches") -> Dict[str, Any]:
        """Generate Sankey flow data."""
        self._validate_dimensions(dimensions)
        self._validate_count_by(count_by)
        
        # Get the data to aggregate  
        if count_by == "sequences":
            sequence_matches = {}
            for match in self.matches:
                if match.sequence_id not in sequence_matches:
                    sequence_matches[match.sequence_id] = match
            data = list(sequence_matches.values())
        else:
            data = self.matches
        
        # Generate nodes and links
        nodes = self._generate_sankey_nodes(data, dimensions)
        links = self._generate_sankey_links(data, dimensions, nodes)
        
        return {
            "dimensions": dimensions,
            "count_by": count_by,
            "total_count": len(data),
            "nodes": nodes,
            "links": links
        }
    
    def _validate_dimensions(self, dimensions: List[str]):
        """Validate that all dimensions exist."""
        invalid = [d for d in dimensions if d not in self.available_dimensions]
        if invalid:
            raise ValueError(f"Invalid dimensions: {invalid}. Available: {self.available_dimensions}")
    
    def _validate_count_by(self, count_by: str):
        """Validate count_by parameter."""
        if count_by not in ["candidate_matches", "sequences"]:
            raise ValueError(f"count_by must be 'candidate_matches' or 'sequences', got: {count_by}")
    
    def _build_tree(self, data: List[CandidateMatch], dimensions: List[str], depth: int) -> Union[int, Dict]:
        """Recursively build hierarchical tree."""
        if depth >= len(dimensions):
            return len(data)
        
        dimension = dimensions[depth]
        grouped = defaultdict(list)
        
        # Group data by current dimension value
        for match in data:
            value = getattr(match, dimension)
            grouped[value].append(match)
        
        # Recursively process each group
        result = {}
        for value, group_data in grouped.items():
            result[value] = self._build_tree(group_data, dimensions, depth + 1)
        
        return result
    
    def _generate_sankey_nodes(self, data: List[CandidateMatch], dimensions: List[str]) -> List[Dict]:
        """Generate nodes for Sankey diagram."""
        nodes = []
        
        for layer_idx, dimension in enumerate(dimensions):
            # Get unique values for this dimension
            values = set(getattr(match, dimension) for match in data)
            
            for value in sorted(values):
                node = {
                    "id": f"{dimension}_{value}",
                    "label": f"{dimension}: {value}",
                    "layer": layer_idx,
                    "dimension": dimension,
                    "value": value
                }
                nodes.append(node)
        
        return nodes
    
    def _generate_sankey_links(self, data: List[CandidateMatch], dimensions: List[str], 
                              nodes: List[Dict]) -> List[Dict]:
        """Generate links for Sankey diagram."""
        if len(dimensions) < 2:
            return []
        
        links = []
        
        # Create node lookup
        node_lookup = {node["id"]: node for node in nodes}
        
        # Generate links between adjacent layers
        for i in range(len(dimensions) - 1):
            source_dim = dimensions[i]
            target_dim = dimensions[i + 1]
            
            # Count flows from source to target
            flow_counts = defaultdict(int)
            
            for match in data:
                source_value = getattr(match, source_dim)
                target_value = getattr(match, target_dim)
                source_id = f"{source_dim}_{source_value}"
                target_id = f"{target_dim}_{target_value}"
                
                flow_counts[(source_id, target_id)] += 1
            
            # Create link objects
            for (source_id, target_id), count in flow_counts.items():
                if source_id in node_lookup and target_id in node_lookup:
                    link = {
                        "source": source_id,
                        "target": target_id,
                        "value": count
                    }
                    links.append(link)
        
        return links


def format_hierarchical_output(stats: Dict[str, Any], indent: str = "     ") -> str:
    """Format hierarchical stats as human-readable text."""
    lines = []
    
    # Header
    dimensions_str = " → ".join(stats["dimensions"])
    lines.append(f"Hierarchical Statistics: {dimensions_str}")
    lines.append(f"Count by: {stats['count_by']}")
    lines.append(f"Total: {stats['total_count']:,}")
    lines.append("")
    
    # Calculate dynamic width for count field based on total
    count_width = len(f"{stats['total_count']:,}")
    
    # Data
    lines.extend(_format_tree_level(stats["data"], stats["dimensions"], 0, indent, stats["total_count"], count_width))
    
    return "\n".join(lines)


def _format_tree_level(data: Union[int, Dict], dimensions: List[str], depth: int, indent: str, parent_total: int, count_width: int) -> List[str]:
    """Recursively format tree levels."""
    if isinstance(data, int):
        return [f"{data:,}"]
    
    lines = []
    dimension_name = dimensions[depth] if depth < len(dimensions) else "value"
    
    # Calculate total for this level
    level_total = _calculate_total(data)
    
    # Sort keys for consistent output
    sorted_keys = sorted(data.keys(), key=str)
    
    for key in sorted_keys:
        value = data[key]
        prefix = indent * depth
        
        if isinstance(value, int):
            # Leaf node - show percentage, count, then label with dimension
            percentage = (value / parent_total * 100) if parent_total > 0 else 0
            lines.append(f"{prefix}{percentage:5.1f}% {value:{count_width},} {key} ({dimension_name})")
        else:
            # Internal node - show percentage, count, then label with dimension
            subtotal = _calculate_total(value)
            percentage = (subtotal / parent_total * 100) if parent_total > 0 else 0
            lines.append(f"{prefix}{percentage:5.1f}% {subtotal:{count_width},} {key} ({dimension_name})")
            lines.extend(_format_tree_level(value, dimensions, depth + 1, indent, subtotal, count_width))
    
    return lines


def _calculate_total(data: Union[int, Dict]) -> int:
    """Calculate total count for a tree node."""
    if isinstance(data, int):
        return data
    
    total = 0
    for value in data.values():
        total += _calculate_total(value)
    return total


def main():
    parser = argparse.ArgumentParser(
        description="Parse Specimux trace events and generate flexible statistics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Hierarchical text output
  %(prog)s trace/ --hierarchical pool primer_pair outcome
  
  # Sequence-level analysis  
  %(prog)s trace/ --hierarchical orientation outcome --count-by sequences
  
  # Sankey flow data (JSON)
  %(prog)s trace/ --sankey-data pool match_type outcome --output flow.json
  
  # List available dimensions
  %(prog)s trace/ --list-dimensions
        """
    )
    
    parser.add_argument('trace_directory', 
                       help='Directory containing trace TSV files')
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--hierarchical', nargs='+', metavar='DIMENSION',
                      help='Generate hierarchical text output with specified dimensions')
    group.add_argument('--sankey-data', nargs='+', metavar='DIMENSION',
                      help='Generate Sankey flow data (JSON) with specified dimensions')
    group.add_argument('--list-dimensions', action='store_true',
                      help='List available dimensions and exit')
    
    parser.add_argument('--count-by', choices=['candidate_matches', 'sequences'], 
                       default='candidate_matches',
                       help='Count by candidate matches or unique sequences (default: candidate_matches)')
    parser.add_argument('--output', '-o', 
                       help='Output file (default: stdout for hierarchical, required for sankey-data)')
    
    args = parser.parse_args()
    
    try:
        # Parse trace files
        parser_obj = TraceEventParser()
        matches = parser_obj.parse_trace_files(args.trace_directory)
        
        if not matches:
            logger.error("No candidate matches found in trace files")
            sys.exit(1)
        
        # Create aggregator
        aggregator = StatsAggregator(matches)
        
        if args.list_dimensions:
            print("Available dimensions:")
            for dim in aggregator.available_dimensions:
                print(f"  {dim}")
            return
        
        if args.hierarchical:
            # Generate hierarchical output
            stats = aggregator.get_hierarchical_stats(args.hierarchical, args.count_by)
            output = format_hierarchical_output(stats)
            
            if args.output:
                with open(args.output, 'w') as f:
                    f.write(output)
                logger.info(f"Hierarchical stats written to {args.output}")
            else:
                print(output)
        
        elif args.sankey_data:
            # Generate Sankey data
            if not args.output:
                logger.error("--output required for --sankey-data")
                sys.exit(1)
            
            data = aggregator.get_sankey_data(args.sankey_data, args.count_by)
            
            with open(args.output, 'w') as f:
                json.dump(data, f, indent=2)
            
            logger.info(f"Sankey data written to {args.output}")
            logger.info(f"Generated {len(data['nodes'])} nodes and {len(data['links'])} links")
    
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()