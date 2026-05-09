#!/usr/bin/env python3
"""
create_visual_latent_learning_data.py <dataset> <index> [--output OUTPUT] [--color-by COLOR_BY] [--hide-edges]
create_visual_latent_learning_data.py <graph_file.pyg> [--output OUTPUT] [--color-by COLOR_BY] [--hide-edges]

Creates an interactive 3D HTML visualization of PyG graph data from feature_store (the format used by the metric learning model).
Truth edges (track_edges) are shown by default.

Usage examples:
  python create_visual_latent_learning_data.py testset 1
    # Visualizes the 1st graph in data/feature_store/testset/ (with edges shown by default)
  
  python create_visual_latent_learning_data.py trainset 1 --color-by particle
    # Colors hits by particle_id (edges still shown by default)
  
  python create_visual_latent_learning_data.py valset 3 --hide-edges
    # Hides truth edges (track_edges)
  
  python create_visual_latent_learning_data.py data/feature_store/trainset/event000000000-graph.pyg
    # Can still use full path if needed

Dependencies: plotly, torch, torch_geometric
"""

import argparse
import sys
from pathlib import Path
import numpy as np

try:
    import torch
    import plotly.graph_objects as go
    import plotly.express as px
    from torch_geometric.data import Data
except ImportError as e:
    print(f"Error: Missing dependency. Install with: pip install plotly torch torch-geometric")
    print(f"Details: {e}")
    sys.exit(1)

# Import visualization utilities
script_dir = Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir))
from visual_utils import build_hit_particle_type_map, get_particle_label


def cylindrical_to_cartesian(r, phi, z):
    """Convert cylindrical coordinates (r, phi, z) to Cartesian (x, y, z)"""
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return x, y, z


def create_visualization(graph, color_by='none', show_edges=False, max_points=None):
    """
    Create interactive 3D visualization of PyG graph

    Args:
        graph: PyG Data object
        color_by: 'none', 'particle', 'region'
        show_edges: If True, draw truth edges (track_edges) as lines
        max_points: Maximum number of points to display (for performance)

    Returns:
        plotly.graph_objects.Figure
    """
    # Build particle type mapping
    particle_type_map = build_hit_particle_type_map(graph)
    # Extract features (what the model actually uses)
    hit_r = graph.hit_r.cpu().numpy()
    hit_phi = graph.hit_phi.cpu().numpy()
    hit_z = graph.hit_z.cpu().numpy()
    
    # Convert to Cartesian for visualization
    x, y, z = cylindrical_to_cartesian(hit_r, hit_phi, hit_z)
    
    # Downsample if needed
    if max_points and len(x) > max_points:
        indices = np.random.choice(len(x), max_points, replace=False)
        x = x[indices]
        y = y[indices]
        z = z[indices]
        hit_r = hit_r[indices]
        hit_phi = hit_phi[indices]
        hit_z = hit_z[indices]
        if 'hit_particle_id' in graph.keys():
            hit_particle_id = graph.hit_particle_id.cpu().numpy()[indices]
        else:
            hit_particle_id = None
        if 'hit_region' in graph.keys():
            hit_region = graph.hit_region.cpu().numpy()[indices]
        else:
            hit_region = None
    else:
        if 'hit_particle_id' in graph.keys():
            hit_particle_id = graph.hit_particle_id.cpu().numpy()
        else:
            hit_particle_id = None
        if 'hit_region' in graph.keys():
            hit_region = graph.hit_region.cpu().numpy()
        else:
            hit_region = None
    
    # Create figure
    fig = go.Figure()
    
    # Determine coloring
    if color_by == 'particle' and hit_particle_id is not None:
        # Color by particle ID
        unique_particles = np.unique(hit_particle_id[hit_particle_id > 0])
        colors = [c for c in px.colors.qualitative.Set3 if c.lower() not in ('#ffffb3', '#ffed6f')]
        for i, pid in enumerate(unique_particles):
            mask = hit_particle_id == pid
            color = colors[i % len(colors)]

            # Get particle label with type
            particle_label = get_particle_label(int(pid), particle_type_map)

            fig.add_trace(go.Scatter3d(
                x=x[mask],
                y=y[mask],
                z=z[mask],
                mode='markers',
                name=particle_label,
                marker=dict(
                    size=2,
                    color=color,
                    opacity=0.7,
                ),
                hovertemplate=f'<b>{particle_label}</b><br>' +
                             'r=%{customdata[0]:.2f}<br>' +
                             'φ=%{customdata[1]:.3f}<br>' +
                             'z=%{customdata[2]:.2f}<br>' +
                             '<extra></extra>',
                customdata=np.column_stack([hit_r[mask], hit_phi[mask], hit_z[mask]]),
            ))
        
        # Show noise hits (particle_id == 0) separately
        if (hit_particle_id == 0).any():
            mask = hit_particle_id == 0
            fig.add_trace(go.Scatter3d(
                x=x[mask],
                y=y[mask],
                z=z[mask],
                mode='markers',
                name='Noise',
                marker=dict(
                    size=1.5,
                    color='gray',
                    opacity=0.5,
                ),
                hovertemplate='<b>Noise</b><br>' +
                             'r=%{customdata[0]:.2f}<br>' +
                             'φ=%{customdata[1]:.3f}<br>' +
                             'z=%{customdata[2]:.2f}<br>' +
                             '<extra></extra>',
                customdata=np.column_stack([hit_r[mask], hit_phi[mask], hit_z[mask]]),
            ))
    
    elif color_by == 'region' and hit_region is not None:
        # Color by region
        unique_regions = np.unique(hit_region)
        colors = [c for c in px.colors.qualitative.Set3 if c.lower() not in ('#ffffb3', '#ffed6f')]
        for i, region in enumerate(unique_regions):
            mask = hit_region == region
            color = colors[i % len(colors)]
            fig.add_trace(go.Scatter3d(
                x=x[mask],
                y=y[mask],
                z=z[mask],
                mode='markers',
                name=f'Region {int(region)}',
                marker=dict(
                    size=2,
                    color=color,
                    opacity=0.7,
                ),
                hovertemplate=f'<b>Region {int(region)}</b><br>' +
                             'r=%{customdata[0]:.2f}<br>' +
                             'φ=%{customdata[1]:.3f}<br>' +
                             'z=%{customdata[2]:.2f}<br>' +
                             '<extra></extra>',
                customdata=np.column_stack([hit_r[mask], hit_phi[mask], hit_z[mask]]),
            ))
    
    else:
        # No coloring - single color
        fig.add_trace(go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode='markers',
            name='Hits',
            marker=dict(
                size=2,
                color='blue',
                opacity=0.7,
            ),
            hovertemplate='<b>Hit</b><br>' +
                         'r=%{customdata[0]:.2f}<br>' +
                         'φ=%{customdata[1]:.3f}<br>' +
                         'z=%{customdata[2]:.2f}<br>' +
                         '<extra></extra>',
            customdata=np.column_stack([hit_r, hit_phi, hit_z]),
        ))
    
    # Add truth edges if requested
    if show_edges and 'track_edges' in graph.keys() and graph.track_edges.numel() > 0:
        track_edges = graph.track_edges.cpu().numpy()
        
        # Get particle IDs for coloring edges
        if hit_particle_id is not None:
            edge_pids = hit_particle_id[track_edges]
            # Only show edges where both nodes belong to same particle
            same_particle = edge_pids[0] == edge_pids[1]
            track_edges = track_edges[:, same_particle]
            edge_pids = edge_pids[0][same_particle]
        
        # Draw edges as lines
        edge_x = []
        edge_y = []
        edge_z = []
        edge_info = []
        
        for i in range(track_edges.shape[1]):
            src_idx = track_edges[0, i]
            tgt_idx = track_edges[1, i]
            
            # Skip if indices are out of bounds (due to downsampling)
            if src_idx >= len(x) or tgt_idx >= len(x):
                continue
            
            edge_x.extend([x[src_idx], x[tgt_idx], None])
            edge_y.extend([y[src_idx], y[tgt_idx], None])
            edge_z.extend([z[src_idx], z[tgt_idx], None])
            
            if hit_particle_id is not None:
                pid = hit_particle_id[src_idx]
                edge_info.append(f'Particle {pid}')
            else:
                edge_info.append('Edge')
        
        # Group edges by particle for coloring
        if hit_particle_id is not None and len(edge_info) > 0:
            unique_particles = np.unique([pid for pid in hit_particle_id[track_edges[0]] if pid > 0])
            colors = [c for c in px.colors.qualitative.Set3 if c.lower() not in ('#ffffb3', '#ffed6f')]
            
            for i, pid in enumerate(unique_particles):
                pid_mask = hit_particle_id[track_edges[0]] == pid
                pid_edges = track_edges[:, pid_mask]

                pid_edge_x = []
                pid_edge_y = []
                pid_edge_z = []

                for j in range(pid_edges.shape[1]):
                    src_idx = pid_edges[0, j]
                    tgt_idx = pid_edges[1, j]
                    if src_idx < len(x) and tgt_idx < len(x):
                        pid_edge_x.extend([x[src_idx], x[tgt_idx], None])
                        pid_edge_y.extend([y[src_idx], y[tgt_idx], None])
                        pid_edge_z.extend([z[src_idx], z[tgt_idx], None])

                if len(pid_edge_x) > 0:
                    # Get particle label with type for edge name
                    particle_label = get_particle_label(int(pid), particle_type_map)

                    fig.add_trace(go.Scatter3d(
                        x=pid_edge_x,
                        y=pid_edge_y,
                        z=pid_edge_z,
                        mode='lines',
                        name=f'Edges ({particle_label})',
                        line=dict(
                            color=colors[i % len(colors)],
                            width=4,
                        ),
                        showlegend=True,
                        hoverinfo='skip',
                    ))
        else:
            # Single color for all edges
            fig.add_trace(go.Scatter3d(
                x=edge_x,
                y=edge_y,
                z=edge_z,
                mode='lines',
                name='Truth Edges',
                line=dict(
                    color='red',
                    width=3,
                ),
                showlegend=True,
                hoverinfo='skip',
            ))
    
    # Update layout
    event_id = graph.event_id[0] if hasattr(graph, "event_id") else "Unknown"
    title_text = f'PyG Graph Visualization - Event {event_id}'
    
    fig.update_layout(
        title=dict(
            text=title_text,
            x=0.5,
            xanchor='center'
        ),
        scene=dict(
            xaxis=dict(title='x (mm)', backgroundcolor="white", gridcolor="lightgray"),
            yaxis=dict(title='y (mm)', backgroundcolor="white", gridcolor="lightgray"),
            zaxis=dict(title='z (mm)', backgroundcolor="white", gridcolor="lightgray"),
            aspectmode='data'
        ),
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01,
            bgcolor="rgba(255, 255, 255, 0.8)"
        ),
        margin=dict(l=0, r=0, b=0, t=40),
        hovermode='closest'
    )
    
    return fig


def main():
    parser = argparse.ArgumentParser(
        description='Create interactive 3D HTML visualization of PyG graph data from feature_store',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python create_visual_latent_learning_data.py testset 1
    # Visualizes the 1st graph in data/feature_store/testset/ (edges shown by default)
  
  python create_visual_latent_learning_data.py trainset 1 --color-by particle
    # Colors hits by particle_id (edges still shown by default)
  
  python create_visual_latent_learning_data.py valset 3 --hide-edges
    # Hides truth edges (track_edges)
  
  python create_visual_latent_learning_data.py data/feature_store/trainset/event000000000-graph.pyg
    # Can still use full path if needed
        """
    )
    
    parser.add_argument(
        'input',
        type=str,
        help='Either: (1) "<dataset> <index>" like "testset 1" or "trainset 3", or (2) full path to .pyg file'
    )
    
    parser.add_argument(
        'index',
        type=str,
        nargs='?',
        default=None,
        help='Index number (if first argument is dataset name). Ignored if first argument is a file path.'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=str,
        default=None,
        help='Output HTML file path (default: data/visuals/latent_data/<dataset>/<filename>.html)'
    )
    
    parser.add_argument(
        '--color-by',
        choices=['none', 'particle', 'region'],
        default='none',
        help='Color points by: none, particle (requires hit_particle_id), or region (default: none)'
    )
    
    parser.add_argument(
        '--hide-edges',
        action='store_false',
        dest='show_edges',
        default=True,
        help='Hide truth edges (default: edges are shown)'
    )
    
    parser.add_argument(
        '--max-points',
        type=int,
        default=None,
        help='Maximum number of points to display (for performance with large datasets)'
    )
    
    args = parser.parse_args()
    
    # Determine graph file path
    script_dir = Path(__file__).resolve().parent
    feature_store_dir = script_dir.parent / 'data' / 'feature_store'
    
    # Check if input is "<dataset> <index>" format
    if args.index is not None:
        # Format: "testset 1" or "trainset 3"
        dataset_name = args.input
        try:
            index = int(args.index)
        except ValueError:
            print(f"Error: Index must be a number, got '{args.index}'")
            sys.exit(1)
        
        dataset_dir = feature_store_dir / dataset_name
        if not dataset_dir.exists():
            print(f"Error: Dataset directory not found: {dataset_dir}")
            print(f"Available datasets: {[d.name for d in feature_store_dir.iterdir() if d.is_dir()]}")
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
    
    print(f"Loading PyG graph from: {graph_path}")
    try:
        graph = torch.load(graph_path, weights_only=False, map_location='cpu')
    except Exception as e:
        print(f"Error loading PyG file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # Check required features
    required_features = ['hit_r', 'hit_phi', 'hit_z']
    missing = [f for f in required_features if f not in graph.keys()]
    if missing:
        print(f"Error: Missing required features: {missing}")
        print(f"Available features: {list(graph.keys())}")
        sys.exit(1)
    
    print(f"Graph loaded: {len(graph.hit_r)} hits")
    if 'hit_particle_id' in graph.keys():
        unique_particles = torch.unique(graph.hit_particle_id[graph.hit_particle_id > 0])
        print(f"  Particles: {len(unique_particles)}")
    if 'track_edges' in graph.keys():
        print(f"  Truth edges: {graph.track_edges.shape[1]}")
    
    # Determine output path
    if args.output is not None:
        output_path = Path(args.output)
    else:
        # Default: save to data/visuals/latent_data/<dataset>/ preserving dataset structure
        # e.g., data/feature_store/testset/event000000000-graph.pyg 
        #   -> data/visuals/latent_data/testset/event000000000-graph.html
        visuals_dir = script_dir.parent / 'data' / 'visuals' / 'latent_data'
        
        # Try to preserve dataset structure (trainset/valset/testset)
        # Check if graph_path is inside data/feature_store/
        try:
            relative_path = graph_path.relative_to(feature_store_dir)
            # relative_path will be like: testset/event000000000-graph.pyg
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
                # Fallback: save to latent_data root
                output_path = visuals_dir / graph_path.name.replace('.pyg', '.html')
        
        # Create visuals directory if it doesn't exist
        output_path.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"\nCreating visualization (color_by={args.color_by}, show_edges={args.show_edges})...")
    try:
        fig = create_visualization(
            graph,
            color_by=args.color_by,
            show_edges=args.show_edges,
            max_points=args.max_points
        )
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
