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
    """Generates plotly Sankey diagrams from trace statistics JSON."""
    
    # Color schemes
    COLOR_SCHEMES = {
        'light': {
            # Pools
            'pool_ITS': '#2E86AB',
            'pool_ITS2': '#A23B72', 
            'pool_none': '#D32F2F',
            
            # Primer pairs (green gradient for full pairs, orange for partial, red for none)
            'primer_full': '#2E7D32',      # Both primers - dark green
            'primer_partial': '#FF9800',   # One primer - orange  
            'primer_none': '#D32F2F',      # No primers - red
            
            # Barcode counts
            'barcode_2': '#1B5E20',        # Both barcodes - darkest green
            'barcode_1': '#8BC34A',        # One barcode - light green
            'barcode_0': '#B71C1C',        # No barcodes - dark red
            
            # Outcomes
            'outcome_matched': '#1B5E20',       # Matched - dark green
            'outcome_partial': '#FF9800',       # Partial - orange
            'outcome_unknown': '#B71C1C',       # Unknown - red
            'outcome_discarded': '#424242',     # Discarded - dark gray
            
            # Default fallbacks
            'default': '#CCCCCC'
        },
        'dark': {
            # Darker theme with more saturated colors
            'pool_ITS': '#4FC3F7',
            'pool_ITS2': '#E91E63', 
            'pool_none': '#F44336',
            
            'primer_full': '#4CAF50',
            'primer_partial': '#FFC107',
            'primer_none': '#F44336',
            
            'barcode_2': '#2E7D32',
            'barcode_1': '#8BC34A',
            'barcode_0': '#E53935',
            
            'outcome_matched': '#2E7D32',
            'outcome_partial': '#FFC107',
            'outcome_unknown': '#E53935',
            'outcome_discarded': '#616161',
            
            'default': '#9E9E9E'
        }
    }
    
    def __init__(self, theme: str = 'light'):
        self.colors = self.COLOR_SCHEMES.get(theme, self.COLOR_SCHEMES['light'])
    
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
        node_colors = [self._get_node_color(node) for node in nodes]
        node_x, node_y = self._calculate_node_positions(nodes, len(dimensions))
        
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
        
        # Create the Sankey diagram
        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=20,
                thickness=25,
                line=dict(color="rgba(0,0,0,0.3)", width=0.5),
                label=node_labels,
                color=node_colors,
                x=node_x,
                y=node_y,
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
    
    def _get_node_color(self, node: Dict[str, Any]) -> str:
        """Determine color for a node based on its dimension and value."""
        dimension = node['dimension']
        value = node['value']
        
        # Pool colors
        if dimension == 'pool':
            return self.colors.get(f'pool_{value}', self.colors['default'])
        
        # Primer pair colors
        elif dimension == 'primer_pair':
            if value == 'none-none':
                return self.colors['primer_none']
            elif 'none' in value:
                return self.colors['primer_partial']
            else:
                return self.colors['primer_full']
        
        # Barcode count colors  
        elif dimension == 'barcode_count':
            return self.colors.get(f'barcode_{value}', self.colors['default'])
        
        # Outcome colors
        elif dimension == 'outcome':
            return self.colors.get(f'outcome_{value}', self.colors['default'])
        
        # Default color for other dimensions
        else:
            return self.colors['default']
    
    def _get_link_color(self, link: Dict[str, Any], nodes: List[Dict]) -> str:
        """Determine color for a link (semi-transparent version of source node)."""
        # Find source node
        source_node = next((n for n in nodes if n['id'] == link['source']), None)
        if not source_node:
            return 'rgba(200,200,200,0.3)'
        
        # Get source color and make it semi-transparent
        source_color = self._get_node_color(source_node)
        
        # Convert hex to rgba with alpha
        if source_color.startswith('#'):
            # Convert hex to rgba
            r = int(source_color[1:3], 16)
            g = int(source_color[3:5], 16)  
            b = int(source_color[5:7], 16)
            return f'rgba({r},{g},{b},0.4)'
        else:
            return 'rgba(200,200,200,0.3)'  # Fallback
    
    def _calculate_node_positions(self, nodes: List[Dict], num_layers: int) -> tuple:
        """Calculate x and y positions for nodes."""
        # Group nodes by layer
        layers = {}
        for node in nodes:
            layer = node['layer']
            if layer not in layers:
                layers[layer] = []
            layers[layer].append(node)
        
        node_x = []
        node_y = []
        
        # Calculate positions
        for node in nodes:
            layer = node['layer']
            layer_nodes = layers[layer]
            
            # X position: evenly distribute layers
            if num_layers > 1:
                x = layer / (num_layers - 1)
            else:
                x = 0.5
            
            # Y position: evenly distribute nodes within layer
            if len(layer_nodes) > 1:
                node_index = layer_nodes.index(node)
                y = node_index / (len(layer_nodes) - 1)
            else:
                y = 0.5
            
            node_x.append(x)
            node_y.append(y)
        
        return node_x, node_y


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