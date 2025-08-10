#!/usr/bin/env python3
"""
Visualization tool for Specimux processing flow statistics.

This script creates interactive Sankey diagrams from Specimux trace files,
showing the flow of sequences through primer detection, outcome classification,
and pool assignment.

The tool can read from two sources:
- Trace files: Reads trace events from TSV files and reconstructs statistics
- Legacy stats.json: Contains pre-aggregated statistics (deprecated)

Trace files provide complete event-level data for full slice-and-dice analysis,
while the legacy format only contains pre-aggregated statistics.

Requires: pip install plotly

Usage: 
    python visualize_stats.py <trace_directory> [output.html]
    python visualize_stats.py trace/
    python visualize_stats.py trace/ flow_diagram.html
    python visualize_stats.py stats.json  # legacy format
"""

import argparse
import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Any, List, Tuple, Set
import ast
import glob

try:
    import plotly.graph_objects as go
    from plotly.offline import plot
except ImportError:
    print("Error: plotly is required. Install with: pip install plotly")
    sys.exit(1)


def read_trace_files(trace_directory: str) -> Dict[Tuple, int]:
    """Read trace files and convert to unified statistics format.
    
    IMPORTANT LIMITATION: Trace files only capture the final selected match per sequence,
    not all potential matches found during the matching process. This means:
    - Primer pair counts will be lower than original (locations vs sequences)  
    - Match attempts count will be estimated from ambiguity candidate counts
    - Barcode combinations will reflect final selections, not all possibilities
    
    The trace format prioritizes event-level analysis over exact statistical reproduction.
    
    Args:
        trace_directory: Directory containing trace TSV files
        
    Returns:
        Dictionary with tuple keys matching the unified stats format (approximate)
    """
    trace_files = glob.glob(os.path.join(trace_directory, "specimux_trace_*.tsv"))
    if not trace_files:
        print(f"No trace files found in {trace_directory}")
        return {}
    
    print(f"Found {len(trace_files)} trace files")
    
    # Initialize counters
    stats = defaultdict(int)
    sequence_outcomes = {}  # sequence_id -> final outcome info
    sequence_events = defaultdict(list)  # sequence_id -> list of events
    
    # Read all trace events
    for trace_file in trace_files:
        print(f"Reading {trace_file}")
        with open(trace_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            try:
                header = next(reader)  # Read header line
            except StopIteration:
                print(f"  Skipping empty file: {trace_file}")
                continue
            
            for row in reader:
                if len(row) < 5:  # Skip incomplete rows
                    continue
                    
                # Parse base fields
                event_dict = {
                    'timestamp': row[0],
                    'worker_id': row[1], 
                    'event_seq': row[2],
                    'sequence_id': row[3],
                    'event_type': row[4]
                }
                
                # Parse event-specific fields based on event type
                if event_dict['event_type'] == 'SEQUENCE_RECEIVED' and len(row) >= 7:
                    event_dict.update({
                        'sequence_length': row[5],
                        'sequence_name': row[6]
                    })
                elif event_dict['event_type'] == 'ORIENTATION_DETECTED' and len(row) >= 9:
                    event_dict.update({
                        'orientation': row[5],
                        'forward_score': row[6],
                        'reverse_score': row[7],
                        'confidence': row[8]
                    })
                elif event_dict['event_type'] == 'MATCH_SELECTED' and len(row) >= 12:
                    event_dict.update({
                        'selection_strategy': row[5],
                        'forward_primer': row[6],
                        'reverse_primer': row[7],
                        'forward_barcode': row[8],
                        'reverse_barcode': row[9],
                        'pool': row[10],
                        'best_match': row[11]
                    })
                elif event_dict['event_type'] == 'SPECIMEN_RESOLVED' and len(row) >= 12:
                    event_dict.update({
                        'specimen_name': row[5],
                        'match_type': row[6],
                        'pool': row[7],
                        'forward_primer': row[8],
                        'reverse_primer': row[9],
                        'forward_barcode': row[10],
                        'reverse_barcode': row[11]
                    })
                elif event_dict['event_type'] == 'SEQUENCE_OUTPUT' and len(row) >= 11:
                    event_dict.update({
                        'output_category': row[5],
                        'specimen_name': row[6],
                        'pool': row[7],
                        'primer_pair': row[8],
                        'output_file': row[9]
                    })
                elif event_dict['event_type'] == 'AMBIGUITY_DETECTED' and len(row) >= 9:
                    event_dict.update({
                        'ambiguity_type': row[5],
                        'candidate_count': row[6],
                        'pools_involved': row[7],
                        'details': row[8] if len(row) > 8 else ''
                    })
                elif event_dict['event_type'] == 'NO_MATCH_FOUND' and len(row) >= 7:
                    event_dict.update({
                        'reason': row[5],
                        'details': row[6] if len(row) > 6 else ''
                    })
                
                sequence_id = event_dict['sequence_id']
                sequence_events[sequence_id].append(event_dict)
    
    print(f"Processed events for {len(sequence_events)} sequences")
    
    # Process each sequence's events to extract statistics
    for sequence_id, events in sequence_events.items():
        # Count sequence processed
        stats[('sequence_processed',)] += 1
        
        # Find orientation event
        orientation = 'unknown'
        for event in events:
            if event['event_type'] == 'ORIENTATION_DETECTED':
                orientation = event['orientation']
                break
        stats[('orientation', orientation)] += 1
        
        # Find final outcome - prioritize MATCH_SELECTED over SPECIMEN_RESOLVED
        final_outcome = None
        final_pool = 'unknown'
        final_primer_pair = 'unknown-unknown'
        final_match_type = 'unknown'
        forward_barcode = 'none'
        reverse_barcode = 'none'
        forward_primer = 'unknown'
        reverse_primer = 'unknown'
        
        # First get match info from MATCH_SELECTED (has the actual primers)
        for event in events:
            if event['event_type'] == 'MATCH_SELECTED':
                final_pool = event.get('pool', 'unknown')
                forward_primer = event.get('forward_primer', 'none')
                reverse_primer = event.get('reverse_primer', 'none')
                forward_barcode = event.get('forward_barcode', 'none')
                reverse_barcode = event.get('reverse_barcode', 'none')
                final_primer_pair = f"{forward_primer}-{reverse_primer}"
                break
        
        # Then get match type from SPECIMEN_RESOLVED
        for event in events:
            if event['event_type'] == 'SPECIMEN_RESOLVED':
                final_match_type = event.get('match_type', 'unknown')
                # Don't override pool/primers from MATCH_SELECTED
                break
        
        # Handle sequences with no specimen resolution (no primers found)
        if final_match_type == 'unknown':
            no_match_found = False
            for event in events:
                if event['event_type'] == 'NO_MATCH_FOUND':
                    no_match_found = True
                    break
            
            if no_match_found:
                stats[('primer_failure', 'no_primers')] += 1
                final_primer_pair = 'none-none'  # Use 'none' not 'unknown'
                final_match_type = 'ambiguous'  # No primers = ambiguous outcome
        
        # Fix primer naming: convert 'unknown' to 'none' in primer names
        if final_primer_pair != 'none-none':
            # Replace 'unknown' with 'none' in primer pair names
            final_primer_pair = final_primer_pair.replace('unknown', 'none')
        
        # Count match attempts - use ambiguity info to estimate total locations
        ambiguity_candidate_count = 1  # Default to 1 if no ambiguity
        for event in events:
            if event['event_type'] == 'AMBIGUITY_DETECTED':
                try:
                    ambiguity_candidate_count = int(event.get('candidate_count', 1))
                except (ValueError, TypeError):
                    ambiguity_candidate_count = 1
                break
        
        stats[('match_attempts',)] += ambiguity_candidate_count
        
        # Count primer pair matches (this represents match attempts/locations)  
        # Use candidate count to approximate the multiple locations
        if final_primer_pair != 'none-none':
            stats[('primer_pair_match', final_pool, final_primer_pair)] += ambiguity_candidate_count
        
        # Determine barcode combination
        if forward_barcode != 'none' and reverse_barcode != 'none':
            barcode_combo = 'both_barcodes'
        elif forward_barcode != 'none':
            barcode_combo = 'forward_only'
        elif reverse_barcode != 'none':
            barcode_combo = 'reverse_only'
        else:
            barcode_combo = 'no_barcodes'
        
        # Count barcode combinations (using candidate count for locations)
        stats[('barcode_combination', final_pool, final_primer_pair, barcode_combo)] += ambiguity_candidate_count
        
        # Determine ambiguity type from AMBIGUITY_DETECTED events
        ambiguity_type = 'none'
        for event in events:
            if event['event_type'] == 'AMBIGUITY_DETECTED':
                ambiguity_type = event['ambiguity_type']
                break
        stats[('ambiguity_type', ambiguity_type)] += 1
        
        # Determine selection outcome
        selection_outcome = 'unique' if ambiguity_type == 'none' else 'ambiguous'
        stats[('selection_outcome', selection_outcome)] += 1
        
        if selection_outcome == 'ambiguous':
            stats[('ambiguity_strategy_applied', 'first')] += 1
        
        # Count retained matches (all sequences are retained)
        stats[('retained_matches',)] += 1
        
        # Map match_type to outcome categories
        outcome_mapping = {
            'full_match': 'matched',
            'partial_forward': 'partial', 
            'partial_reverse': 'partial',
            'unknown': 'unknown'
        }
        
        if final_match_type in outcome_mapping:
            outcome = outcome_mapping[final_match_type]
        elif ambiguity_type != 'none':
            outcome = 'ambiguous'
        else:
            outcome = 'unknown'
        
        # Count final outcomes
        stats[('outcome', final_pool, final_primer_pair, outcome)] += 1
    
    return dict(stats)


def convert_raw_stats_to_flow(raw_stats: Dict[str, Any]) -> Dict[str, Any]:
    """Convert raw unified statistics to flow diagram format.
    
    Takes raw unified statistics with string keys (converted from tuples) and
    converts them into a flow diagram with nodes and links.
    
    Args:
        raw_stats: Dictionary with 'metadata' and 'unified_stats' keys
        
    Returns:
        Dictionary with 'metadata', 'nodes', and 'links' keys for flow diagram
    """
    # Extract unified_stats and convert string keys back to tuples
    unified_stats = {}
    for str_key, value in raw_stats['unified_stats'].items():
        # Convert string representation back to tuple
        try:
            tuple_key = ast.literal_eval(str_key)
            unified_stats[tuple_key] = value
        except (ValueError, SyntaxError):
            print(f"Warning: Could not parse key: {str_key}")
            continue
    
    if not unified_stats:
        return {
            "metadata": raw_stats['metadata'].copy(),
            "nodes": [],
            "links": []
        }
    
    def calculate_proportional_positions(node_values, start_y=0.15, available_height=0.8, padding_fraction=0.02):
        """Calculate proportional y-positions for nodes based on their values.
        
        Args:
            node_values: List of (name, value) tuples in desired order
            start_y: Y coordinate to start from (reserves space for title)
            available_height: Total height available for nodes
            padding_fraction: Fraction of available height to use as padding between nodes
        
        Returns:
            Dictionary mapping node names to y positions (centers of nodes)
        """
        if not node_values:
            return {}
        
        if len(node_values) == 1:
            return {node_values[0][0]: start_y + available_height / 2}
        
        total_value = sum(value for _, value in node_values)
        if total_value == 0:
            # Fallback to even spacing if all values are 0
            positions = {}
            for i, (name, _) in enumerate(node_values):
                y_pos = start_y + (i * available_height) / (len(node_values) - 1)
                positions[name] = y_pos
            return positions
        
        total_padding = padding_fraction * available_height * (len(node_values) - 1)
        available_for_nodes = available_height - total_padding
        
        positions = {}
        current_y = start_y
        
        for name, value in node_values:
            # Calculate proportional height for this node
            node_height = (value / total_value) * available_for_nodes
            
            # Position at center of the node's allocated space
            node_center_y = current_y + node_height / 2
            positions[name] = node_center_y
            
            # Move to start of next node (including padding)
            current_y += node_height + (padding_fraction * available_height)
        
        return positions
    
    nodes = []
    links = []
    
    # Define color scheme for visualization
    colors = {
        # Layer 1 - Input
        'input': '#2E86AB',
        'filtered': '#B71C1C',
        
        # Layer 2 - Orientation (updated: forward and reverse are both green - successful orientation)
        'orient_forward': '#4CAF50',
        'orient_reverse': '#4CAF50', 
        'orient_unknown': '#FF9800',
        
        # Layer 3 - Primer Pairs (updated colors)
        'primer_full': '#2E7D32',    # Both primers found - dark green
        'primer_partial': '#FF9800', # One primer found - orange
        'primer_none': '#D32F2F',    # No primers found - red
        
        # Layer 4 - Barcode Combinations
        'barcode_both': '#1B5E20',     # Both barcodes
        'barcode_forward': '#8BC34A',  # Forward only
        'barcode_reverse': '#8BC34A',  # Reverse only
        'barcode_none': '#B71C1C',     # No barcodes
        
        # Layer 5 - Outcomes
        'outcome_matched': '#1B5E20',
        'outcome_partial': '#FF9800',
        'outcome_ambiguous': '#9C27B0',
        'outcome_unknown': '#B71C1C',
        
        # Layer 6 - Pools
        'pool_ITS': '#1976D2',
        'pool_ITS2': '#00796B',
        'pool_unknown': '#D32F2F'
    }
    
    # Define layer x-coordinates
    layer_x = {
        'input': 0.0,      # Layer 1 - Input
        'orientation': 0.2, # Layer 2 - Orientation
        'primer': 0.4,     # Layer 3 - Primer Pairs
        'barcode': 0.6,    # Layer 4 - Barcode Combinations
        'outcome': 0.8,    # Layer 5 - Outcomes
        'pool': 1.0        # Layer 6 - Pools
    }
    
    # Layer 1: Input
    total_input = unified_stats.get(('sequence_processed',), 0)
    nodes.append({
        "id": "input", 
        "name": f"Input ({total_input:,} seqs)",
        "color": colors['input'],
        "x": layer_x['input'],
        "y": 0.5  # Center vertically
    })
    
    # Handle filtered sequences (they don't go through the pipeline)
    filtered_short = unified_stats.get(('filtered', 'too_short'), 0)
    filtered_long = unified_stats.get(('filtered', 'too_long'), 0)
    total_filtered = filtered_short + filtered_long
    if total_filtered > 0:
        nodes.append({
            "id": "filtered", 
            "name": f"Filtered ({total_filtered:,} seqs)",
            "color": colors['filtered'],
            "x": layer_x['input'],
            "y": 0.8  # Below input
        })
        links.append({"source": "input", "target": "filtered", "value": total_filtered})
    
    # Calculate sequences that go through processing
    processed_sequences = total_input - total_filtered
    
    # Layer 2: Orientation Decisions
    orientation_totals = {}
    for key, count in unified_stats.items():
        if key[0] == 'orientation':
            orientation_totals[key[1]] = count
    
    # Create orientation nodes and links from input (ordered: forward, reverse, unknown)
    orientation_order = ['forward', 'reverse', 'unknown']
    
    # Calculate proportional positions for orientation nodes
    orientation_values = []
    for orientation in orientation_order:
        if orientation in orientation_totals:
            orientation_values.append((orientation, orientation_totals[orientation]))
    
    orientation_y_positions = calculate_proportional_positions(orientation_values)
    
    for orientation in orientation_order:
        if orientation in orientation_totals:
            count = orientation_totals[orientation]
            node_id = f"orient_{orientation}"
            color_key = f"orient_{orientation}"
            nodes.append({
                "id": node_id, 
                "name": f"Orientation: {orientation.title()} ({count:,} seqs)",
                "color": colors.get(color_key, colors['orient_unknown']),
                "x": layer_x['orientation'],
                "y": orientation_y_positions[orientation]
            })
            links.append({"source": "input", "target": node_id, "value": count})
    
    # Layer 3: Primer Pairs Detected
    # Collect all primer pairs from barcode combinations and outcomes
    primer_pair_totals = {}
    
    # From barcode combinations (this includes both full and partial matches)
    for key, count in unified_stats.items():
        if key[0] == 'barcode_combination':  # ('barcode_combination', pool, primer_pair, combo)
            primer_pair = key[2]
            if primer_pair not in primer_pair_totals:
                primer_pair_totals[primer_pair] = 0
            primer_pair_totals[primer_pair] += count
    
    # Add sequences with no primers detected
    no_primers_count = unified_stats.get(('primer_failure', 'no_primers'), 0)
    if no_primers_count > 0:
        primer_pair_totals['no-primers'] = no_primers_count
    
    # Helper function to determine primer pair color based on unknown primers
    def get_primer_pair_color(primer_pair):
        if primer_pair == 'no-primers':
            return colors['primer_none']
        elif 'unknown-unknown' in primer_pair:
            return colors['primer_none']  # Red - no valid primers
        elif 'unknown' in primer_pair:
            return colors['primer_partial']  # Orange - one primer unknown
        else:
            return colors['primer_full']  # Green - both primers valid
    
    # Create primer pair nodes (ordered: both-primers, single-primer, no-primers)
    # Sort primer pairs to get proper ordering
    primer_pairs_sorted = []
    both_primers = []
    single_primer = []
    no_primers = []
    
    for primer_pair in primer_pair_totals.keys():
        if primer_pair == 'no-primers' or 'unknown-unknown' in primer_pair:
            no_primers.append(primer_pair)
        elif 'unknown' in primer_pair:
            single_primer.append(primer_pair)
        else:
            both_primers.append(primer_pair)
    
    primer_pairs_sorted = sorted(both_primers) + sorted(single_primer) + sorted(no_primers)
    
    # Calculate proportional positions for primer pairs
    primer_values = [(pp, primer_pair_totals[pp]) for pp in primer_pairs_sorted]
    primer_y_positions = calculate_proportional_positions(primer_values)
    
    for primer_pair in primer_pairs_sorted:
        count = primer_pair_totals[primer_pair]
        node_id = f"primer_{primer_pair}"
        color = get_primer_pair_color(primer_pair)
        
        if primer_pair == 'no-primers':
            display_name = f"No Primers ({count:,} locs)"
        else:
            display_name = f"{primer_pair} ({count:,} locs)"
        
        nodes.append({
            "id": node_id, 
            "name": display_name,
            "color": color,
            "x": layer_x['primer'],
            "y": primer_y_positions[primer_pair]
        })
    
    # Create links from orientation to primer pairs
    # For now, distribute proportionally based on total counts
    # This is a simplification - ideally we'd track orientation->primer_pair flows directly
    total_orientation_count = sum(orientation_totals.values())
    for orientation, orientation_count in orientation_totals.items():
        orientation_id = f"orient_{orientation}"
        
        # Distribute each primer pair proportionally across orientations
        for primer_pair, pp_count in primer_pair_totals.items():
            # Calculate proportional flow for this orientation
            if total_orientation_count > 0:
                flow_value = int(pp_count * (orientation_count / total_orientation_count))
                if flow_value > 0:
                    primer_id = f"primer_{primer_pair}"
                    links.append({"source": orientation_id, "target": primer_id, "value": flow_value})
    
    # Layer 4: Barcode Match Combinations (simplified - no pool/primer info)
    # Aggregate counts across all primer pairs for each barcode combination type
    barcode_totals = {
        'both_barcodes': 0,
        'forward_only': 0, 
        'reverse_only': 0,
        'no_barcodes': 0
    }
    
    for key, count in unified_stats.items():
        if key[0] == 'barcode_combination':  # ('barcode_combination', pool, primer_pair, combo)
            combo = key[3]
            if combo in barcode_totals:
                barcode_totals[combo] += count
    
    # Add no-primers to no_barcodes category
    if no_primers_count > 0:
        barcode_totals['no_barcodes'] += no_primers_count
    
    # Create simplified barcode combination nodes (ordered: both, single (forward/reverse), neither)
    barcode_display_names = {
        'both_barcodes': 'Both Barcodes',
        'forward_only': 'Forward Only',
        'reverse_only': 'Reverse Only', 
        'no_barcodes': 'No Barcodes'
    }
    
    barcode_order = ['both_barcodes', 'forward_only', 'reverse_only', 'no_barcodes']
    
    # Calculate proportional positions for barcode combinations
    barcode_values = []
    for combo in barcode_order:
        count = barcode_totals.get(combo, 0)
        if count > 0:  # Only include nodes that have actual data
            barcode_values.append((combo, count))
    
    barcode_y_positions = calculate_proportional_positions(barcode_values)
    
    for combo in barcode_order:
        count = barcode_totals.get(combo, 0)
        if count > 0:  # Only create nodes that have actual data
            node_id = f"barcode_{combo}"
            display_name = barcode_display_names[combo]
            
            # Determine color based on barcode combination
            if combo == 'both_barcodes':
                color = colors['barcode_both']
            elif combo in ['forward_only', 'reverse_only']:
                color = colors['barcode_forward'] if combo == 'forward_only' else colors['barcode_reverse']
            else:  # no_barcodes
                color = colors['barcode_none']
                
            nodes.append({
                "id": node_id,
                "name": f"{display_name} ({count:,} locs)",
                "color": color,
                "x": layer_x['barcode'],
                "y": barcode_y_positions[combo]
            })
    
    # Create links from primer pairs to barcode combinations
    # Since we simplified barcode combinations, we need to aggregate the flows
    primer_to_barcode_flows = {}
    
    for key, count in unified_stats.items():
        if key[0] == 'barcode_combination':  # ('barcode_combination', pool, primer_pair, combo)
            primer_pair = key[2]
            combo = key[3]
            
            primer_id = f"primer_{primer_pair}"
            barcode_id = f"barcode_{combo}"
            
            flow_key = (primer_id, barcode_id)
            if flow_key not in primer_to_barcode_flows:
                primer_to_barcode_flows[flow_key] = 0
            primer_to_barcode_flows[flow_key] += count
    
    # Add no-primers flow to no_barcodes
    if no_primers_count > 0:
        flow_key = ("primer_no-primers", "barcode_no_barcodes")
        primer_to_barcode_flows[flow_key] = no_primers_count
    
    # Add the aggregated flows as links
    for (primer_id, barcode_id), count in primer_to_barcode_flows.items():
        primer_exists = any(n['id'] == primer_id for n in nodes)
        barcode_exists = any(n['id'] == barcode_id for n in nodes)
        if primer_exists and barcode_exists:
            links.append({"source": primer_id, "target": barcode_id, "value": count})
    
    # Layer 5: Match Outcomes
    outcome_totals = {}
    for key, count in unified_stats.items():
        if key[0] == 'outcome':  # ('outcome', pool, primer_pair, match_type)
            match_type = key[3]
            if match_type not in outcome_totals:
                outcome_totals[match_type] = 0
            outcome_totals[match_type] += count
    
    # Create outcome nodes (ordered: matched, partial, ambiguous, unknown)
    outcome_order = ['matched', 'partial', 'ambiguous', 'unknown']
    
    # Calculate proportional positions for outcome nodes
    outcome_values = []
    for match_type in outcome_order:
        count = outcome_totals.get(match_type, 0)
        if count > 0:  # Only include nodes that have actual data
            outcome_values.append((match_type, count))
    
    outcome_y_positions = calculate_proportional_positions(outcome_values)
    
    for match_type in outcome_order:
        count = outcome_totals.get(match_type, 0)
        if count > 0:  # Only create nodes that have actual data
            node_id = f"outcome_{match_type}"
            color_key = f"outcome_{match_type}"
            nodes.append({
                "id": node_id, 
                "name": f"Outcome: {match_type.title()} ({count:,} seqs)",
                "color": colors.get(color_key, colors['outcome_unknown']),
                "x": layer_x['outcome'],
                "y": outcome_y_positions[match_type]
            })
    
    # Create links from barcode combinations to outcomes
    # With simplified barcode combinations, we need to map outcomes to the right barcode types
    barcode_to_outcome_flows = {}
    
    for key, count in unified_stats.items():
        if key[0] == 'outcome':  # ('outcome', pool, primer_pair, match_type)
            pool = key[1]
            primer_pair = key[2]
            match_type = key[3]
            
            # Map outcomes to their expected barcode combinations
            if match_type == 'matched':
                barcode_type = 'both_barcodes'
            elif match_type == 'partial':
                # Partial matches come from forward_only or reverse_only barcodes
                # We need to check which barcode combinations contributed to this outcome
                for bc_type in ['forward_only', 'reverse_only']:
                    bc_count = unified_stats.get(('barcode_combination', pool, primer_pair, bc_type), 0)
                    if bc_count > 0:
                        flow_key = (f"barcode_{bc_type}", f"outcome_{match_type}")
                        if flow_key not in barcode_to_outcome_flows:
                            barcode_to_outcome_flows[flow_key] = 0
                        barcode_to_outcome_flows[flow_key] += bc_count
                continue
            else:  # unknown, ambiguous
                barcode_type = 'no_barcodes'
            
            flow_key = (f"barcode_{barcode_type}", f"outcome_{match_type}")
            if flow_key not in barcode_to_outcome_flows:
                barcode_to_outcome_flows[flow_key] = 0
            barcode_to_outcome_flows[flow_key] += count
    
    # Add the flows as links
    for (barcode_id, outcome_id), count in barcode_to_outcome_flows.items():
        barcode_exists = any(n['id'] == barcode_id for n in nodes)
        outcome_exists = any(n['id'] == outcome_id for n in nodes)
        if barcode_exists and outcome_exists:
            links.append({"source": barcode_id, "target": outcome_id, "value": count})
    
    # Layer 6: Final Pool Assignment
    pool_totals = {}
    for key, count in unified_stats.items():
        if key[0] == 'outcome':  # ('outcome', pool, primer_pair, match_type)
            pool = key[1]
            if pool not in pool_totals:
                pool_totals[pool] = 0
            pool_totals[pool] += count
    
    # Create pool nodes (ordered with unknown last)
    known_pools = []
    unknown_pools = []
    
    for pool in pool_totals.keys():
        if pool == 'None' or pool.lower() == 'unknown':
            unknown_pools.append(pool)
        else:
            known_pools.append(pool)
    
    pool_order = sorted(known_pools) + sorted(unknown_pools)
    
    # Calculate proportional positions for pool nodes
    pool_values = [(pool, pool_totals[pool]) for pool in pool_order]
    pool_y_positions = calculate_proportional_positions(pool_values)
    
    for pool in pool_order:
        count = pool_totals[pool]
        node_id = f"pool_{pool}"
        pool_display = pool if pool != 'None' else 'unknown'
        
        # Determine color based on pool type
        pool_key = f'pool_{pool_display}'
        if pool_key in colors:
            color = colors[pool_key]
        else:
            color = colors['pool_unknown']
            
        nodes.append({
            "id": node_id, 
            "name": f"Pool: {pool_display.upper() if pool_display != 'unknown' else 'Unknown'} ({count:,} seqs)",
            "color": color,
            "x": layer_x['pool'],
            "y": pool_y_positions[pool]
        })
    
    # Create links from outcomes to pools
    for key, count in unified_stats.items():
        if key[0] == 'outcome':  # ('outcome', pool, primer_pair, match_type)
            pool = key[1]
            match_type = key[3]
            
            outcome_id = f"outcome_{match_type}"
            pool_id = f"pool_{pool}"
            
            # Only add link if both nodes exist
            outcome_exists = any(n['id'] == outcome_id for n in nodes)
            pool_exists = any(n['id'] == pool_id for n in nodes)
            if outcome_exists and pool_exists:
                links.append({"source": outcome_id, "target": pool_id, "value": count})
    
    # Create flow stats with updated metadata
    flow_stats = {
        "metadata": {
            "title": "Specimux Detailed Processing Flow",
            "description": "6-layer flow: Input → Orientation → Primer Pairs → Barcode Combinations → Outcomes → Pools",
            "total_sequences": total_input,
            "generated_by": "visualize_stats (from raw specimux data)",
            "timestamp": raw_stats['metadata'].get('timestamp', 'unknown'),
            "layers": [
                "Input sequences",
                "Orientation decisions", 
                "Primer pair detection",
                "Barcode match combinations",
                "Classification outcomes",
                "Pool assignments"
            ],
            "unit_explanation": {
                "sequences": "Individual DNA sequences being processed (seqs)",
                "locations": "Potential primer-pair matches - amplification occurs here as one sequence can match multiple primer pairs (locs)",
                "note": "Numbers may increase from Layer 2 to Layer 3 due to amplification (one sequence matching multiple primer pairs), then decrease back to sequence counts in subsequent layers."
            }
        },
        "nodes": nodes,
        "links": links
    }
    
    return flow_stats


def load_stats_data(input_path: str) -> Dict[str, Any]:
    """Load and validate stats from trace files or legacy JSON."""
    if os.path.isdir(input_path):
        # This is a directory - look for trace files
        print("Detected trace directory - reading trace files...")
        unified_stats = read_trace_files(input_path)
        
        if not unified_stats:
            print("Error: No trace data found")
            sys.exit(1)
        
        # Create metadata
        total_sequences = unified_stats.get(('sequence_processed',), 0)
        metadata = {
            "title": "Specimux Raw Statistics",
            "description": "Raw unified statistics reconstructed from trace events",
            "total_sequences": total_sequences,
            "generated_by": "visualize_stats (from trace files)",
            "timestamp": "reconstructed from trace"
        }
        
        # Convert to string keys for compatibility with existing flow conversion
        raw_stats = {
            "metadata": metadata,
            "unified_stats": {str(k): v for k, v in unified_stats.items()}
        }
        
        return convert_raw_stats_to_flow(raw_stats)
        
    elif os.path.isfile(input_path):
        # This is a file - assume it's a JSON file (legacy format)
        print("Detected file - assuming legacy JSON format...")
        try:
            with open(input_path, 'r') as f:
                data = json.load(f)
        except FileNotFoundError:
            print(f"Error: File '{input_path}' not found.")
            sys.exit(1)
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in '{input_path}': {e}")
            sys.exit(1)
        
        # Check if this is raw stats or flow data format
        if 'unified_stats' in data:
            # This is raw stats format - convert to flow format
            print("Detected raw stats format - converting to flow diagram...")
            return convert_raw_stats_to_flow(data)
        elif 'nodes' in data and 'links' in data:
            # This is already flow format - validate and return
            print("Detected flow diagram format - using directly...")
            required_keys = ['metadata', 'nodes', 'links']
            for key in required_keys:
                if key not in data:
                    print(f"Error: Missing required key '{key}' in JSON file")
                    sys.exit(1)
            
            if 'title' not in data['metadata'] or 'total_sequences' not in data['metadata']:
                print("Error: Missing required metadata fields (title, total_sequences)")
                sys.exit(1)
                
            return data
        else:
            print("Error: JSON file must contain either 'unified_stats' (raw format) or 'nodes'+'links' (flow format)")
            sys.exit(1)
    else:
        print(f"Error: '{input_path}' is neither a file nor a directory")
        sys.exit(1)


def create_node_colors(nodes: List[Dict[str, str]]) -> List[str]:
    """Extract colors from JSON nodes, with fallback for missing color fields."""
    colors = []
    
    for node in nodes:
        # Use color from JSON node if available, otherwise fall back to default
        node_color = node.get('color', '#CCCCCC')
        colors.append(node_color)
            
    return colors


def create_flow_diagram(data: Dict[str, Any], output_file: str, width: int = 1200, height: int = 600) -> None:
    """Create a Sankey diagram from specimux flow statistics."""
    
    # Extract nodes and create index mapping
    nodes = data['nodes']
    if not nodes:
        print("Error: No nodes found in the data")
        sys.exit(1)
    
    node_ids = [node['id'] for node in nodes]
    node_names = [node['name'] for node in nodes]
    
    # Extract x and y coordinates if available
    node_x = []
    node_y = []
    for node in nodes:
        node_x.append(node.get('x', None))  # Use None if no x coordinate
        node_y.append(node.get('y', None))  # Use None if no y coordinate
    
    # Create hover text with additional details
    hover_text = []
    for node in nodes:
        # Extract count from node name if available
        name = node['name']
        hover_info = f"ID: {node['id']}<br>Name: {name}"
        if 'color' in node:
            hover_info += f"<br>Color: {node['color']}"
        if 'x' in node and 'y' in node:
            hover_info += f"<br>Position: ({node['x']:.1f}, {node['y']:.1f})"
        hover_text.append(hover_info)
    
    # Create index mapping
    id_to_index = {node_id: i for i, node_id in enumerate(node_ids)}
    
    # Extract links with indices
    links = data['links']
    if not links:
        print("Warning: No links found in the data")
    
    try:
        source_indices = [id_to_index[link['source']] for link in links]
        target_indices = [id_to_index[link['target']] for link in links]
        values = [link['value'] for link in links]
    except KeyError as e:
        print(f"Error: Invalid node reference in links: {e}")
        available_nodes = list(id_to_index.keys())
        print(f"Available node IDs: {available_nodes}")
        sys.exit(1)
    
    # Generate colors from JSON nodes
    node_colors = create_node_colors(nodes)
    
    # Validate that we have colors for all nodes
    if len(node_colors) != len(nodes):
        print(f"Warning: Color count ({len(node_colors)}) doesn't match node count ({len(nodes)})")
    
    # Create node configuration with x and y coordinates if available
    node_config = dict(
        pad=20,  # Increased spacing between nodes
        thickness=25,  # Slightly thicker nodes
        line=dict(color="black", width=0.5),
        label=node_names,
        color=node_colors,
        hovertemplate='%{label}<br>%{customdata}<extra></extra>',
        customdata=hover_text
    )
    
    # Add x and y coordinates if all nodes have them
    if all(x is not None for x in node_x) and all(y is not None for y in node_y):
        node_config['x'] = node_x
        node_config['y'] = node_y
        print("Using explicit node positioning")
    else:
        print("Using automatic node positioning (some nodes missing coordinates)")
    
    # Create the Sankey diagram with enhanced styling
    fig = go.Figure(data=[go.Sankey(
        node=node_config,
        link=dict(
            source=source_indices,
            target=target_indices,
            value=values,
            hovertemplate='%{source.label} → %{target.label}<br>Count: %{value:,}<extra></extra>'
        )
    )])
    
    # Set title and layout with enhanced formatting
    title = data['metadata']['title']
    total_seqs = data['metadata']['total_sequences']
    description = data['metadata'].get('description', '')
    
    # Create subtitle with description if available
    subtitle = f"({total_seqs:,} sequences processed)"
    if description:
        subtitle = f"{description}<br>{subtitle}"
    
    fig.update_layout(
        title_text=f"{title}<br><sub>{subtitle}</sub>",
        title_font_size=16,
        font_size=11,
        width=width,
        height=height,
        margin=dict(t=80, b=40, l=40, r=40),  # Better margins
        hoverlabel=dict(
            bgcolor="white",
            bordercolor="black",
            font_size=12,
            font_family="monospace"
        )
    )
    
    # Save the diagram
    plot(fig, filename=output_file, auto_open=False)
    print(f"Flow diagram saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Create interactive Sankey diagrams from Specimux flow statistics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s trace/
  %(prog)s trace/ custom_output.html
  %(prog)s trace/ --width 1600 --height 800
  %(prog)s stats.json  # legacy format
        """
    )
    
    parser.add_argument('input_path', 
                       help='Path to trace directory or legacy stats.json file')
    parser.add_argument('output_file', nargs='?', 
                       help='Output HTML file (default: flow_diagram.html)')
    parser.add_argument('--width', type=int, default=1200,
                       help='Diagram width in pixels (default: 1200)')
    parser.add_argument('--height', type=int, default=600,
                       help='Diagram height in pixels (default: 600)')
    
    args = parser.parse_args()
    
    # Validate input path
    if not os.path.exists(args.input_path):
        print(f"Error: Input path '{args.input_path}' does not exist")
        sys.exit(1)
    
    # Set default output file if not provided
    output_file = args.output_file or 'flow_diagram.html'
    
    # Validate output directory is writable
    output_dir = os.path.dirname(os.path.abspath(output_file))
    if not os.access(output_dir, os.W_OK):
        print(f"Error: Cannot write to output directory '{output_dir}'")
        sys.exit(1)
    
    # Load data and create diagram
    data = load_stats_data(args.input_path)
    create_flow_diagram(data, output_file, args.width, args.height)


if __name__ == "__main__":
    main()