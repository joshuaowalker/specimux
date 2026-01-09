#!/usr/bin/env python3
"""
Sankey diagram generator for trace-based statistics.

This tool takes JSON output from trace_to_stats.py and generates interactive
plotly Sankey diagrams. Keeps visualization logic separate from stats aggregation.

Usage:
    python stats_to_sankey.py sankey_data.json output.html
    python stats_to_sankey.py sankey_data.json --width 1600 --height 800
    python stats_to_sankey.py sankey_data.json --theme dark
"""

import argparse
import json
import sys
from typing import Dict, List, Any
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

try:
    import plotly.graph_objects as go
    from plotly.offline import plot
except ImportError:
    logger.error("plotly is required. Install with: pip install plotly")
    sys.exit(1)


class SankeyVisualizer:
    """Generates plotly Sankey diagrams from trace statistics JSON.
    
    Design Philosophy:
    - Modular separation: Complete independence from statistics computation
    - Semantic coloring: Use processing pipeline knowledge for intuitive colors
    - Universal applicability: Work with any dimensions without hardcoded logic
    - Robust layout: Use plotly's adaptive algorithms instead of brittle positioning
    
    Coloring Heuristics:
    - Success gradient: Red (failed) → Orange (partial) → Yellow → Green (complete)
    - Semantic mapping: forward=green, reverse=blue, discarded=gray
    - Pool diversity: Hash-based color wheel for arbitrary user pool names
    - Pattern recognition: none/partial/full primers get red/orange/green respectively
    - Consistent assignment: Same values get same colors across runs via hashing
    """
    
    # Base color palettes for heuristic coloring
    PALETTES = {
        'light': {
            # Success gradient: red -> orange -> yellow -> green
            'success_colors': ['#D32F2F', '#FF9800', '#FFC107', '#4CAF50', '#2E7D32'],
            # Distinct color wheel for arbitrary categorical values (pools, primers, etc.)
            'categorical_colors': [
                '#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#6A994E', 
                '#A663CC', '#4ECDC4', '#FF6B6B', '#4ECDC4', '#45B7D1'
            ],
            # Orientation colors
            'orientation_forward': '#4CAF50',   # Green - forward processing
            'orientation_reverse': '#2196F3',   # Blue - reverse processing  
            'orientation_unknown': '#FF9800',   # Orange - uncertain
            # Special semantic colors
            'none_missing': '#D32F2F',          # Red for none/missing values
            'discarded': '#424242',             # Gray for discarded items
            'default': '#9E9E9E'                # Light gray fallback
        },
        'dark': {
            'success_colors': ['#F44336', '#FFC107', '#FFEB3B', '#8BC34A', '#4CAF50'],
            'categorical_colors': [
                '#4FC3F7', '#E91E63', '#FF9800', '#E53935', '#8BC34A',
                '#BA68C8', '#4DB6AC', '#FF7043', '#81C784', '#64B5F6'
            ],
            'orientation_forward': '#8BC34A',
            'orientation_reverse': '#64B5F6', 
            'orientation_unknown': '#FFC107',
            'none_missing': '#F44336',
            'discarded': '#616161',
            'default': '#BDBDBD'
        }
    }
    
    def __init__(self, theme: str = 'light'):
        self.theme = theme
        self.palette = self.PALETTES.get(theme, self.PALETTES['light'])
        # Cache for consistent pool colors
        self._pool_color_cache = {}
    
    def generate_diagram(self, data: Dict[str, Any], output_file: str, 
                        width: int = 1200, height: int = 600) -> None:
        """Generate Sankey diagram from JSON data."""
        
        # Validate input data
        self._validate_data(data)
        
        # Extract data
        nodes = data['nodes']
        links = data['links']
        dimensions = data['dimensions']
        count_by = data['count_by']
        total_count = data['total_count']
        
        logger.info(f"Generating Sankey with {len(nodes)} nodes and {len(links)} links")
        
        # Prepare node data for plotly
        node_labels = [node['label'] for node in nodes]
        node_colors = [self._get_node_color(node, nodes) for node in nodes]
        
        # Use plotly's automatic layout instead of explicit positioning
        # This is more robust and adapts better to different data
        
        # Create node index mapping
        node_index = {node['id']: i for i, node in enumerate(nodes)}
        
        # Prepare link data for plotly
        source_indices = []
        target_indices = []
        values = []
        link_colors = []
        
        for link in links:
            if link['source'] in node_index and link['target'] in node_index:
                source_indices.append(node_index[link['source']])
                target_indices.append(node_index[link['target']])
                values.append(link['value'])
                link_colors.append(self._get_link_color(link, nodes))
        
        # Create hover text for nodes
        hover_text = []
        for i, node in enumerate(nodes):
            # Calculate total flow through this node
            total_in = sum(link['value'] for link in links if link['target'] == node['id'])
            total_out = sum(link['value'] for link in links if link['source'] == node['id'])
            
            # For source nodes, use outgoing flow; for sink nodes, use incoming flow
            flow_count = total_out if total_out > 0 else total_in
            
            hover_info = f"<b>{node['label']}</b><br>"
            hover_info += f"Layer: {node['layer'] + 1}<br>"
            hover_info += f"Flow: {flow_count:,} {count_by.replace('_', ' ')}"
            hover_text.append(hover_info)
        
        # Create the Sankey diagram with automatic layout
        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=20,
                thickness=25,
                line=dict(color="rgba(0,0,0,0.3)", width=0.5),
                label=node_labels,
                color=node_colors,
                hovertemplate='%{customdata}<extra></extra>',
                customdata=hover_text
            ),
            link=dict(
                source=source_indices,
                target=target_indices,
                value=values,
                color=link_colors,
                hovertemplate='<b>%{source.label}</b> → <b>%{target.label}</b><br>' +
                             'Flow: %{value:,} ' + count_by.replace('_', ' ') + '<extra></extra>'
            )
        )])
        
        # Set title and layout
        dimensions_str = " → ".join(dimensions)
        title = f"Specimux Flow Diagram: {dimensions_str}"
        subtitle = f"({total_count:,} {count_by.replace('_', ' ')} processed)"
        
        fig.update_layout(
            title=dict(
                text=f"{title}<br><sub>{subtitle}</sub>",
                font_size=16
            ),
            font_size=11,
            width=width,
            height=height,
            margin=dict(t=80, b=40, l=40, r=40),
            hoverlabel=dict(
                bgcolor="white",
                bordercolor="black", 
                font_size=12,
                font_family="monospace"
            )
        )
        
        # Save the diagram
        plot(fig, filename=output_file, auto_open=False)
        logger.info(f"Sankey diagram saved to {output_file}")
    
    def _validate_data(self, data: Dict[str, Any]) -> None:
        """Validate input JSON data structure."""
        required_keys = ['dimensions', 'count_by', 'total_count', 'nodes', 'links']
        for key in required_keys:
            if key not in data:
                raise ValueError(f"Missing required key in JSON data: {key}")
        
        # Validate nodes
        if not data['nodes']:
            raise ValueError("No nodes found in data")
        
        node_required = ['id', 'label', 'layer', 'dimension', 'value']
        for i, node in enumerate(data['nodes']):
            for key in node_required:
                if key not in node:
                    raise ValueError(f"Node {i} missing required key: {key}")
        
        # Validate links
        if data['links']:  # Links might be empty for single-dimension data
            link_required = ['source', 'target', 'value']
            for i, link in enumerate(data['links']):
                for key in link_required:
                    if key not in link:
                        raise ValueError(f"Link {i} missing required key: {key}")
    
    def _get_node_color(self, node: Dict[str, Any], all_nodes: List[Dict]) -> str:
        """Determine color using semantic heuristics based on dimension and value."""
        dimension = node['dimension']
        value = str(node['value'])
        
        # Pool colors: Use categorical color wheel for arbitrary pool names
        if dimension == 'pool':
            return self._get_pool_color(value)
        
        # Primer-related dimensions: Pattern-based coloring
        elif dimension in ['primer_pair', 'forward_primer', 'reverse_primer']:
            return self._get_primer_color(value)
        
        # Barcode dimensions: Count or presence-based
        elif dimension == 'barcode_count':
            return self._get_barcode_count_color(value)
        elif dimension in ['forward_barcode', 'reverse_barcode', 'barcode_pair']:
            return self._get_barcode_color(value)
        elif dimension in ['forward_barcode_matched', 'reverse_barcode_matched']:
            return self._get_boolean_color(value)
        
        # Orientation: Semantic direction colors
        elif dimension == 'orientation':
            return self._get_orientation_color(value)
        
        # Match type: Success gradient based on completeness
        elif dimension == 'match_type':
            return self._get_match_type_color(value)
        
        # Outcome/resolution: Success gradient with detailed discard reasons
        elif dimension in ['outcome', 'resolution_type', 'outcome_detailed']:
            return self._get_outcome_color(value)
        
        # Selection strategy: Process-based coloring
        elif dimension == 'selection_strategy':
            return self._get_selection_color(value)
        
        # Boolean dimensions
        elif dimension == 'specimen_matched':
            return self._get_boolean_color(value)
        
        # Categorical dimensions: Use color wheel based on unique values in layer
        else:
            return self._get_categorical_color(dimension, value, all_nodes)
    
    def _get_pool_color(self, pool_name: str) -> str:
        """Get color for pool using consistent categorical assignment."""
        if pool_name.lower() in ['none', 'unknown']:
            return self.palette['none_missing']
        
        # Use cached color for consistency
        if pool_name not in self._pool_color_cache:
            # Assign color based on hash for consistency across runs
            color_idx = hash(pool_name) % len(self.palette['categorical_colors'])
            self._pool_color_cache[pool_name] = self.palette['categorical_colors'][color_idx]
        
        return self._pool_color_cache[pool_name]
    
    def _get_primer_color(self, primer_value: str) -> str:
        """Color primers based on completeness pattern."""
        if primer_value.lower() in ['none', 'none-none']:
            return self.palette['none_missing']
        elif 'none' in primer_value.lower() or 'unknown' in primer_value.lower():
            return self.palette['success_colors'][1]  # Orange - partial
        else:
            return self.palette['success_colors'][4]  # Green - complete
    
    def _get_barcode_count_color(self, count_str: str) -> str:
        """Color barcode counts on success gradient."""
        try:
            count = int(count_str)
            if count == 0:
                return self.palette['success_colors'][0]  # Red
            elif count == 1:
                return self.palette['success_colors'][2]  # Yellow
            else:  # 2 or more
                return self.palette['success_colors'][4]  # Green
        except ValueError:
            return self.palette['default']
    
    def _get_barcode_color(self, barcode_value: str) -> str:
        """Color individual barcodes."""
        if barcode_value.lower() in ['none', 'unknown']:
            return self.palette['none_missing']
        else:
            return self.palette['success_colors'][3]  # Light green - present
    
    def _get_orientation_color(self, orientation: str) -> str:
        """Color orientations with semantic meaning."""
        orientation_lower = orientation.lower()
        if orientation_lower == 'forward':
            return self.palette['orientation_forward']
        elif orientation_lower == 'reverse':
            return self.palette['orientation_reverse']
        else:
            return self.palette['orientation_unknown']
    
    def _get_match_type_color(self, match_type: str) -> str:
        """Color match types based on completeness."""
        match_lower = match_type.lower()
        if match_lower == 'both':
            return self.palette['success_colors'][4]  # Green - both found
        elif match_lower in ['forward_only', 'reverse_only']:
            return self.palette['success_colors'][2]  # Yellow - partial
        else:
            return self.palette['success_colors'][0]  # Red - none
    
    def _get_outcome_color(self, outcome: str) -> str:
        """Color outcomes on success gradient with detailed discard reason support."""
        outcome_lower = outcome.lower()
        if outcome_lower == 'matched':
            return self.palette['success_colors'][4]  # Green - success
        elif outcome_lower == 'partial':
            return self.palette['success_colors'][2]  # Yellow - partial success
        elif outcome_lower.startswith('discarded'):
            # Handle detailed discard reasons with color variations
            if 'lower_score' in outcome_lower:
                return '#757575'  # Medium gray - competitive discard
            else:
                return self.palette['discarded']  # Default gray - general discard
        else:  # unknown
            return self.palette['success_colors'][0]  # Red - failed
    
    def _get_selection_color(self, strategy: str) -> str:
        """Color selection strategies."""
        strategy_lower = strategy.lower()
        if strategy_lower == 'unique':
            return self.palette['success_colors'][4]  # Green - ideal case
        elif strategy_lower == 'first':
            return self.palette['success_colors'][2]  # Yellow - ambiguity resolved
        elif strategy_lower == 'discarded':
            return self.palette['discarded']          # Gray - not selected
        else:
            return self.palette['default']
    
    def _get_boolean_color(self, bool_str: str) -> str:
        """Color boolean values."""
        if bool_str.lower() in ['true', '1', 'yes']:
            return self.palette['success_colors'][4]  # Green - true
        else:
            return self.palette['success_colors'][0]  # Red - false
    
    def _get_categorical_color(self, dimension: str, value: str, all_nodes: List[Dict]) -> str:
        """Assign colors to arbitrary categorical values using color wheel."""
        # Get all unique values for this dimension to ensure consistent assignment
        dimension_values = sorted(set(
            str(node['value']) for node in all_nodes if node['dimension'] == dimension
        ))
        
        try:
            value_index = dimension_values.index(value)
            color_index = value_index % len(self.palette['categorical_colors'])
            return self.palette['categorical_colors'][color_index]
        except ValueError:
            return self.palette['default']
    
    def _get_link_color(self, link: Dict[str, Any], nodes: List[Dict]) -> str:
        """Determine color for a link (semi-transparent version of source node)."""
        # Find source node
        source_node = next((n for n in nodes if n['id'] == link['source']), None)
        if not source_node:
            return 'rgba(200,200,200,0.3)'
        
        # Get source color and make it semi-transparent
        source_color = self._get_node_color(source_node, nodes)
        
        # Convert hex to rgba with alpha
        if source_color.startswith('#'):
            # Convert hex to rgba
            r = int(source_color[1:3], 16)
            g = int(source_color[3:5], 16)  
            b = int(source_color[5:7], 16)
            return f'rgba({r},{g},{b},0.4)'
        else:
            return 'rgba(200,200,200,0.3)'  # Fallback
    


def main():
    parser = argparse.ArgumentParser(
        description="Generate interactive Sankey diagrams from trace statistics JSON",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s sankey_data.json output.html
  %(prog)s sankey_data.json output.html --width 1600 --height 800
  %(prog)s sankey_data.json output.html --theme dark
        """
    )
    
    parser.add_argument('json_file', 
                       help='JSON file from trace_to_stats.py --sankey-data')
    parser.add_argument('output_file', nargs='?', default='sankey_diagram.html',
                       help='Output HTML file (default: sankey_diagram.html)')
    parser.add_argument('--width', type=int, default=1200,
                       help='Diagram width in pixels (default: 1200)')
    parser.add_argument('--height', type=int, default=600,  
                       help='Diagram height in pixels (default: 600)')
    parser.add_argument('--theme', choices=['light', 'dark'], default='light',
                       help='Color theme (default: light)')
    
    args = parser.parse_args()
    
    try:
        # Load JSON data
        with open(args.json_file, 'r') as f:
            data = json.load(f)
        
        logger.info(f"Loaded data with {len(data['nodes'])} nodes and {len(data['links'])} links")
        
        # Generate diagram
        visualizer = SankeyVisualizer(theme=args.theme)
        visualizer.generate_diagram(data, args.output_file, args.width, args.height)
        
    except FileNotFoundError:
        logger.error(f"JSON file not found: {args.json_file}")
        sys.exit(1)
    except json.JSONDecodeError as e:
        logger.error(f"Invalid JSON file: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()