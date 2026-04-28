#!/usr/bin/env python3
"""
create_visual_gnn_data.py <dataset> <index> [latent] [--output OUTPUT]
create_visual_gnn_data.py <graph_file.pyg> [--output OUTPUT]

Creates an interactive 3D HTML visualization of a graph from the GNN pipeline.
Displays all nodes (hits) and edges, with edges colored by their truth labels:
  - Red: true edges (hits from same particle)
  - Black: false edges (candidate connections)

Usage examples:
  python create_visual_gnn_data.py trainset 1
    # Visualizes the 1st graph in trainset (learned latent space graphs - DEFAULT)
  
  python create_visual_gnn_data.py trainset 1 simple
    # Visualizes the 1st graph in trainset (simple radius+KNN graphs)
  
  python create_visual_gnn_data.py valset 3
    # Visualizes the 3rd graph in valset (learned latent space graphs - DEFAULT)
  
  python create_visual_gnn_data.py data/graph_constructed_latent/trainset/event000000001-graph.pyg
    # Can still use full path if needed
  
Dependencies: torch, torch_geometric, plotly, pandas, numpy
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import torch
import yaml

try:
    import plotly.graph_objects as go
    import plotly.express as px
    import pandas as pd
except ImportError:
    print("Error: plotly and pandas are required. Install with: pip install plotly pandas")
    sys.exit(1)

from visual_utils import (
    cylindrical_to_cartesian,
    get_standard_scene_layout,
    create_node_hover_text,
    build_hit_particle_type_map,
    get_particle_label
)


def load_graph(graph_path):
    """Load a PyTorch Geometric graph from file."""
    try:
        graph = torch.load(graph_path, weights_only=False)
    except Exception as e:
        raise RuntimeError(f"Error loading graph file: {e}")


    return graph


def create_visualization(graph, r_max=None, k_max=None):
    """
    Create an interactive Plotly figure of the graph.
    Visualizes the graph exactly as constructed - showing all nodes (hits) and edges
    that will be used as input to the GNN model.
    
    Args:
        graph: PyTorch Geometric Data object with nodes and edges
        r_max: Radius parameter used in graph construction (for title)
        k_max: Max neighbors parameter used in graph construction (for title)
    
    Returns:
        plotly.graph_objects.Figure
    """
    # Extract node positions (convert to numpy)
    r = graph.r.numpy()
    phi = graph.phi.numpy()
    z_cyl = graph.z.numpy()
    
    # Convert cylindrical to Cartesian coordinates for plotting
    x, y, z = cylindrical_to_cartesian(r, phi, z_cyl)
    
    num_nodes = len(x)
    
    # Extract edge information
    edge_index = graph.edge_index.numpy()
    edge_y = graph.edge_y.numpy()
    
    # Separate true and false edges
    true_edges_mask = (edge_y == 1)
    false_edges_mask = (edge_y == 0)

    true_edges_raw = edge_index[:, true_edges_mask]
    false_edges = edge_index[:, false_edges_mask]

    # Deduplicate true edges: (i,j) and (j,i) are the same line in 3D.
    # edge_y labels both directions true when one_in_one_out=True, so canonicalize.
    if true_edges_raw.shape[1] > 0:
        canonical = np.stack([true_edges_raw.min(axis=0), true_edges_raw.max(axis=0)])
        _, unique_idx = np.unique(canonical.T, axis=0, return_index=True)
        true_edges = true_edges_raw[:, unique_idx]
    else:
        true_edges = true_edges_raw

    num_true = true_edges.shape[1]
    num_false = false_edges.shape[1]
    
    # Create figure
    fig = go.Figure()
    
    # Add false edges (plot first so they're in background)
    edge_x_false = []
    edge_y_false = []
    edge_z_false = []
    
    for i in range(false_edges.shape[1]):
        src = false_edges[0, i]
        dst = false_edges[1, i]
        edge_x_false.extend([x[src], x[dst], None])
        edge_y_false.extend([y[src], y[dst], None])
        edge_z_false.extend([z[src], z[dst], None])
    
    fig.add_trace(go.Scatter3d(
        x=edge_x_false,
        y=edge_y_false,
        z=edge_z_false,
        mode='lines',
        line=dict(color='black', width=1),
        opacity=0.1,
        name=f'False edges ({num_false})',
        hoverinfo='skip',
        showlegend=True
    ))
    
    # Add true edges grouped by particle ID (if available)
    # Try to get particle_id if available (but only if it matches num_nodes)
    particle_ids = None
    if hasattr(graph, 'particle_id') and graph.particle_id is not None:
        pid_tensor = graph.particle_id
        if len(pid_tensor) == num_nodes:
            particle_ids = pid_tensor.numpy()
        else:
            print(f"Warning: particle_id length ({len(pid_tensor)}) doesn't match num_nodes ({num_nodes}), skipping particle grouping")
            particle_ids = None

    # Build particle type mapping (particle_id → PDG code) for labels
    particle_type_map = build_hit_particle_type_map(graph)

    if particle_ids is not None:

        # Group true edges by particle ID
        # Get particle IDs for source nodes of true edges
        true_edge_src_pids = particle_ids[true_edges[0]]
        unique_particles = np.unique(true_edge_src_pids)
        unique_particles = unique_particles[unique_particles > 0]  # Remove noise (pid == 0)

        # Use qualitative colors for different particles
        colors = px.colors.qualitative.Set3

        # Create a trace for each particle
        for i, pid in enumerate(unique_particles):
            # Find edges belonging to this particle
            pid_mask = true_edge_src_pids == pid
            pid_edges = true_edges[:, pid_mask]

            if pid_edges.shape[1] == 0:
                continue

            # Build edge coordinates for this particle
            pid_edge_x = []
            pid_edge_y = []
            pid_edge_z = []

            for j in range(pid_edges.shape[1]):
                src = pid_edges[0, j]
                dst = pid_edges[1, j]
                pid_edge_x.extend([x[src], x[dst], None])
                pid_edge_y.extend([y[src], y[dst], None])
                pid_edge_z.extend([z[src], z[dst], None])

            color = colors[i % len(colors)]

            # Generate label with particle species if available
            particle_label = get_particle_label(int(pid), particle_type_map)
            legend_label = f'{particle_label} ({pid_edges.shape[1]} edges)'

            fig.add_trace(go.Scatter3d(
                x=pid_edge_x,
                y=pid_edge_y,
                z=pid_edge_z,
                mode='lines',
                line=dict(color=color, width=2),
                opacity=0.7,
                name=legend_label,
                hoverinfo='skip',
                showlegend=True
            ))
    else:
        # Fallback: single trace for all true edges if particle_id not available
        edge_x_true = []
        edge_y_true = []
        edge_z_true = []
        
        for i in range(true_edges.shape[1]):
            src = true_edges[0, i]
            dst = true_edges[1, i]
            edge_x_true.extend([x[src], x[dst], None])
            edge_y_true.extend([y[src], y[dst], None])
            edge_z_true.extend([z[src], z[dst], None])
        
        fig.add_trace(go.Scatter3d(
            x=edge_x_true,
            y=edge_y_true,
            z=edge_z_true,
            mode='lines',
            line=dict(color='red', width=2),
            opacity=0.6,
            name=f'True edges ({num_true})',
            hoverinfo='skip',
            showlegend=True
        ))
    
    # Add nodes with hover text
    # Extract time attributes if available
    hit_t = graph.hit_t.numpy() if hasattr(graph, 'hit_t') else None
    hit_t_smeared = graph.hit_t_smeared.numpy() if hasattr(graph, 'hit_t_smeared') else None
    hit_has_time = graph.hit_has_time.numpy() if hasattr(graph, 'hit_has_time') else None

    hover_text = []
    for i in range(num_nodes):
        additional_fields = {}
        if hit_t is not None:
            additional_fields['t'] = float(hit_t[i])
        if hit_t_smeared is not None:
            additional_fields['t_smeared'] = float(hit_t_smeared[i])
        if hit_has_time is not None:
            additional_fields['has_time'] = int(hit_has_time[i])
        if particle_ids is not None:
            pid = int(particle_ids[i])
            additional_fields['particle_id'] = pid
            # Add particle species if available
            if pid in particle_type_map:
                additional_fields['particle_species'] = get_particle_label(pid, particle_type_map)

        hover_text.append(
            create_node_hover_text(i, {
                'x': x[i], 'y': y[i], 'z_cart': z[i],
                'r': r[i], 'phi': phi[i]
            }, additional_fields=additional_fields)
        )
    
    fig.add_trace(go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode='markers',
        marker=dict(
            size=4,
            color='steelblue',
            opacity=0.8,
            line=dict(width=0.5, color='white')
        ),
        name=f'Nodes ({num_nodes})',
        text=hover_text,
        hoverinfo='text',
        showlegend=True
    ))
    
    # Update layout with standard scene configuration
    title_text = f"Graph: {num_nodes} nodes, {num_true} true edges, {num_false} false edges"
    if r_max is not None and k_max is not None:
        title_text += f"<br>Graph Construction: r_max={r_max}, k_max={k_max}"

    fig.update_layout(**get_standard_scene_layout(title_text, margin_t=40, legend_opacity=0.8))
    
    return fig


def main():
    parser = argparse.ArgumentParser(
        description='Create interactive 3D HTML visualization of a graph from the GNN pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python create_visual_gnn_data.py trainset 1
    # Visualizes the 1st graph in data/graph_constructed_latent/trainset/
  
  python create_visual_gnn_data.py trainset 1 latent
    # Visualizes the 1st graph in data/graph_constructed_latent/trainset/
  
  python create_visual_gnn_data.py valset 3
    # Visualizes the 3rd graph in data/graph_constructed_latent/valset/
  
  python create_visual_gnn_data.py data/graph_constructed_latent/trainset/event000000003-graph.pyg
    # Can still use full path if needed
        """
    )
    parser.add_argument(
        'input',
        type=str,
        help='Either: (1) "<dataset> <index>" like "trainset 1" or "valset 3", or (2) full path to .pyg file'
    )
    parser.add_argument(
        'index',
        type=str,
        nargs='?',
        default=None,
        help='Index number (if first argument is dataset name). Ignored if first argument is a file path.'
    )
    parser.add_argument(
        'mode',
        type=str,
        nargs='?',
        default=None,
        help='Optional: "simple" to use graph_constructed directory instead of graph_constructed_latent (default: latent)'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        default=None,
        help='Output HTML file path (default: data/visuals/gnn_data/<dataset>/<filename>.html)'
    )
    
    args = parser.parse_args()
    
    # Determine if using latent graphs (default is True now)
    use_latent = True
    if args.mode is not None and args.mode.lower() == 'simple':
        use_latent = False

    # Determine graph file path
    script_dir = Path(__file__).resolve().parent
    if use_latent:
        graph_constructed_dir = script_dir.parent / 'data' / 'graph_constructed_latent'
        config_file = 'latent_stage_(1)/graph_construction_latent.yaml'
    else:
        graph_constructed_dir = script_dir.parent / 'data' / 'graph_constructed'
        config_file = 'graph_construction.yaml'
    
    # Check if input is "<dataset> <index>" format
    if args.index is not None:
        # Format: "trainset 1" or "valset 3"
        dataset_name = args.input
        try:
            index = int(args.index)
        except ValueError:
            print(f"Error: Index must be a number, got '{args.index}'")
            sys.exit(1)
        
        dataset_dir = graph_constructed_dir / dataset_name
        if not dataset_dir.exists():
            print(f"Error: Dataset directory not found: {dataset_dir}")
            print(f"Available datasets: {[d.name for d in graph_constructed_dir.iterdir() if d.is_dir()]}")
            sys.exit(1)
        
        # Get all .pyg files sorted
        graph_files = sorted([f for f in dataset_dir.glob('*.pyg')])
        if len(graph_files) == 0:
            print(f"Error: No graph files found in {dataset_dir}")
            sys.exit(1)
        
        if index < 1 or index > len(graph_files):
            print(f"Error: Index {index} out of range. Available: 1-{len(graph_files)}")
            sys.exit(1)
        
        graph_path = graph_files[index - 1]  # Convert to 0-indexed
        print(f"Selected: {graph_path.name} (index {index} of {len(graph_files)} in {dataset_name})")
    else:
        # Format: full path to file
        graph_path = Path(args.input)
        if not graph_path.is_absolute():
            # Try relative to script directory
            graph_path = script_dir / graph_path
    
    # Validate input file
    if not graph_path.exists():
        print(f"Error: Graph file not found: {graph_path}")
        sys.exit(1)
    
    # Determine output path
    if args.output is not None:
        output_path = Path(args.output)
    else:
        # Default: save to data/visuals/gnn_data/<dataset>/ preserving dataset structure
        # e.g., data/graph_constructed_latent/trainset/event000000001-graph.pyg 
        #   -> data/visuals/gnn_data/trainset/event000000001-graph.html
        visuals_dir = script_dir.parent / 'data' / 'visuals' / 'gnn_data'
        
        # Try to preserve dataset structure (trainset/valset/testset)
        # Check if graph_path is inside data/graph_constructed_latent/ or data/graph_constructed/
        try:
            relative_path = graph_path.relative_to(graph_constructed_dir)
            # relative_path will be like: trainset/event000000001-graph.pyg
            output_path = visuals_dir / relative_path.with_suffix('.html')
        except ValueError:
            # If not in expected structure, try to extract dataset from path or use default
            # Check if any parent directory is trainset/valset/testset
            dataset_name = None
            for parent in graph_path.parents:
                if parent.name in ['trainset', 'valset', 'testset']:
                    dataset_name = parent.name
                    break
            
            if dataset_name:
                output_path = visuals_dir / dataset_name / graph_path.name.replace('.pyg', '.html')
            else:
                # Fallback: save to gnn_data root
                output_path = visuals_dir / graph_path.name.replace('.pyg', '.html')
        
        # Create visuals directory if it doesn't exist
        output_path.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Loading graph from: {graph_path}")
    try:
        graph = load_graph(graph_path)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
    
    print(f"Graph loaded successfully:")
    print(f"  • Nodes: {graph.num_nodes}")
    print(f"  • Edges: {graph.edge_index.shape[1]}")
    _te = graph.edge_index[:, graph.edge_y.bool()]
    if _te.shape[1] > 0:
        _canon = np.stack([_te.min(dim=0)[0].numpy(), _te.max(dim=0)[0].numpy()])
        _num_unique_true = np.unique(_canon.T, axis=0).shape[0]
    else:
        _num_unique_true = 0
    print(f"  • True edges: {_num_unique_true} ({100*_num_unique_true/graph.edge_index.shape[1]:.1f}%)")
    
    # Load graph construction parameters from config
    config_path = script_dir.parent / 'acorn_configs' / config_file
    r_max = None
    k_max = None
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
            if use_latent:
                # Latent config has k_max directly
                k_max = config.get('k_max')
                print(f"  • Graph construction: Learned latent space, k_max={k_max}")
            else:
                # Simple config has graph_construction section
                if 'graph_construction' in config:
                    r_max = config['graph_construction'].get('r_max')
                    k_max = config['graph_construction'].get('k_max')
                    print(f"  • Graph construction: Simple radius+KNN, r_max={r_max}, k_max={k_max}")
    except Exception as e:
        print(f"  • Warning: Could not load graph construction params from config: {e}")
    
    print("\nCreating visualization...")
    try:
        fig = create_visualization(graph, r_max=r_max, k_max=k_max)
    except Exception as e:
        print(f"Error creating visualization: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    print(f"Writing interactive HTML to: {output_path}")
    try:
        fig.write_html(str(output_path), include_plotlyjs='cdn')
    except Exception as e:
        print(f"Error writing HTML file: {e}")
        sys.exit(1)
    
    print("\n✓ Visualization complete!")
    print(f"\nOpen in browser: {output_path.absolute()}")


if __name__ == '__main__':
    main()


