#!/usr/bin/env python3
"""
create_visual_track_building_proposal.py - Visualize proposed track clusters

Visualize proposed track clusters from Connected Components algorithm.
Shows the actual reconstructed tracks with each cluster in a different color
for easy visual examination. Uses hit_track_labels to identify track membership.

This script loads graphs from data/track_building/ with track_edges and
hit_track_labels from the Connected Components clustering stage.

Edge colors: Each track cluster shown in a unique color from the palette
Isolated hits: Shown in gray with reduced opacity

Usage:
  python create_visual_track_building_proposal.py <dataset> <index> [--output FILE]

Examples:
  python create_visual_track_building_proposal.py trainset 1
    # Visualizes 1st graph from data/track_building/trainset/

  python create_visual_track_building_proposal.py valset 3
    # Visualizes 3rd validation graph

  python create_visual_track_building_proposal.py testset 5
    # Visualizes 5th test graph

  python create_visual_track_building_proposal.py trainset 1 --output custom/path.html
    # Custom output location

This shows the direct output of the Connected Components clustering algorithm,
with each track cluster rendered in a unique color for individual examination.
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import torch

try:
    import plotly.graph_objects as go
    import plotly.express as px
except ImportError:
    print("Error: plotly is required. Install with: pip install plotly")
    sys.exit(1)

from visual_utils import (
    validate_dataset_name,
    extract_node_coordinates,
    get_standard_scene_layout,
    create_node_hover_text
)


def load_track_building_graph(dataset_name, index, input_dir=None):
    """
    Load a graph from the track building output directory.

    Args:
        dataset_name: 'trainset', 'valset', or 'testset'
        index: Index of graph (1-based)
        input_dir: Directory containing track_building data (default: ../data/track_building)

    Returns:
        graph: PyTorch Geometric Data object with edge_index, scores, and hit_track_labels
        graph_path: Path to the loaded graph file
    """
    script_dir = Path(__file__).resolve().parent

    if input_dir is None:
        input_dir = script_dir.parent / 'data' / 'track_building'
    else:
        input_dir = Path(input_dir)

    # Validate dataset name
    validate_dataset_name(dataset_name)

    dataset_dir = input_dir / dataset_name
    if not dataset_dir.exists():
        raise FileNotFoundError(f"Dataset directory not found: {dataset_dir}")

    # Get all .pyg files sorted
    graph_files = sorted([f for f in dataset_dir.glob('*.pyg')])
    if len(graph_files) == 0:
        raise ValueError(f"No graph files found in {dataset_dir}")

    if index < 1 or index > len(graph_files):
        raise ValueError(f"Index {index} out of range. Available: 1-{len(graph_files)}")

    # Select the graph file at the specified index (convert to 0-indexed)
    graph_path = graph_files[index - 1]
    print(f"Loading {dataset_name} graph #{index} from: {graph_path.name}")

    # Load the graph file
    try:
        graph = torch.load(graph_path, map_location='cpu', weights_only=False)
    except Exception as e:
        raise RuntimeError(f"Failed to load graph file {graph_path}: {e}")

    # Check required attributes
    if not hasattr(graph, 'edge_index'):
        raise ValueError(
            f"\nGraph does not have edge_index attribute.\n"
            f"  Graph file: {graph_path}\n"
            f"  Solution: Re-run track building:\n"
            f"    python helix_segmentmatcher_walkthrough.py testset\n"
        )

    if not hasattr(graph, 'scores') and not hasattr(graph, 'edge_scores'):
        raise ValueError(
            f"\nGraph does not have scores or edge_scores attribute (GNN edge scores).\n"
            f"  Graph file: {graph_path}\n"
            f"  Solution: Re-run GNN inference and track building:\n"
            f"    python infer_gnn.py\n"
            f"    python helix_segmentmatcher_walkthrough.py testset\n"
        )

    if not hasattr(graph, 'hit_track_labels'):
        raise ValueError(
            f"\nGraph does not have hit_track_labels attribute.\n"
            f"  Graph file: {graph_path}\n"
            f"  Solution: Re-run track building:\n"
            f"    python helix_segmentmatcher_walkthrough.py testset\n"
        )

    return graph, graph_path


def visualize_track_building_graph(dataset_name, index, input_dir=None, output_path=None):
    """
    Load and visualize a graph with built tracks from track building stage.

    Args:
        dataset_name: Name of dataset ('trainset', 'valset', or 'testset')
        index: Index of graph in dataset (1-based)
        input_dir: Directory containing track_building data
        output_path: Output HTML file path
    """
    print("=" * 70)
    print(f"VISUALIZING TRACK CLUSTERS: {dataset_name} graph #{index}")
    print("=" * 70)
    print()

    # Load graph
    graph, graph_path = load_track_building_graph(dataset_name, index, input_dir)

    # Count unique tracks and isolated hits
    track_labels = graph.hit_track_labels.tolist()
    unique_tracks = set(track_labels)
    if -1 in unique_tracks:
        unique_tracks.remove(-1)
    num_tracks = len(unique_tracks)
    num_isolated = track_labels.count(-1)

    print(f"Graph loaded:")
    print(f"  Nodes: {graph.num_nodes}")
    print(f"  Track clusters: {num_tracks}")
    print(f"  Isolated hits: {num_isolated}")
    print()

    # Create visualization
    print("Creating visualization...")
    fig = create_track_visualization(graph, dataset_name, index)

    # Determine output path
    script_dir = Path(__file__).resolve().parent
    if output_path is None:
        visuals_dir = script_dir.parent / 'data' / 'visuals' / 'track_building' / dataset_name
        visuals_dir.mkdir(parents=True, exist_ok=True)
        output_path = visuals_dir / f"{dataset_name}{index:03d}_track_clusters.html"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Writing interactive HTML to: {output_path}")
    fig.write_html(str(output_path), include_plotlyjs='cdn')

    print()
    print("=" * 70)
    print("✓ Visualization complete!")
    print("=" * 70)
    print(f"\nOpen in browser: {output_path.absolute()}")

    return output_path


def create_track_visualization(graph, dataset_name, index):
    """
    Create interactive Plotly visualization showing track clusters.

    Each track cluster is shown in a unique color based on hit_track_labels.
    Shows the actual edges used by the Connected Components clustering algorithm
    (GNN-filtered edges with score > score_cut).

    Args:
        graph: PyTorch Geometric Data object with edge_index, scores, and hit_track_labels
        dataset_name: Name of dataset ('trainset', 'valset', or 'testset')
        index: Index of graph in dataset

    Returns:
        plotly.graph_objects.Figure
    """
    # Extract node positions and convert to Cartesian
    coords = extract_node_coordinates(graph, return_cylindrical=True, return_cartesian=True)
    x, y, z_cart = coords['x'], coords['y'], coords['z_cart']
    r, phi, z = coords['r'], coords['phi'], coords['z']
    num_nodes = len(x)

    # Load score_cut threshold from config
    import yaml
    script_dir = Path(__file__).resolve().parent
    config_path = script_dir.parent / 'acorn_configs' / 'track_building_stage_(3)' / 'helix_segmentmatcher_walkthrough.yaml'
    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)
        score_cut = config.get('score_cut', 0.5)
    except FileNotFoundError:
        score_cut = 0.5
        print(f"Warning: Config not found, using default score_cut={score_cut}")

    # Get the edges that Connected Components actually used
    # These are the GNN-filtered edges (score > score_cut)
    # Handle both 'scores' and 'edge_scores' attribute names
    if hasattr(graph, 'scores'):
        edge_scores = np.array(graph.scores.tolist())
    else:
        edge_scores = np.array(graph.edge_scores.tolist())
    edge_index = np.array(graph.edge_index.tolist())

    # Filter edges by score (same as Connected Components does)
    edge_mask = edge_scores > score_cut
    clustering_edges = edge_index[:, edge_mask]

    # Extract track labels
    hit_track_labels = np.array(graph.hit_track_labels.tolist())

    # Get unique track IDs (excluding isolated hits with label=-1)
    unique_tracks = sorted(set(hit_track_labels.tolist()))
    if -1 in unique_tracks:
        unique_tracks.remove(-1)

    num_tracks = len(unique_tracks)
    num_isolated = (hit_track_labels == -1).sum()

    # Create color palette (22 unique colors)
    # Combine Plotly (10) + Set3 (12) palettes for 22 distinct colors
    base_colors = px.colors.qualitative.Plotly  # 10 colors
    extended_colors = px.colors.qualitative.Set3  # 12 colors
    all_colors = base_colors + extended_colors  # 22 colors total

    # Map track_id -> color
    track_colors = {
        track_id: all_colors[i % len(all_colors)]
        for i, track_id in enumerate(unique_tracks)
    }

    # Group clustering edges by their track cluster
    # Store (src, dst, score) tuples to include score in hover text
    track_edge_groups = {track_id: [] for track_id in unique_tracks}
    mismatched = 0
    cross_cluster = 0  # Edges connecting different clusters

    # Get indices in original edge_index to retrieve scores
    filtered_edge_indices = np.where(edge_mask)[0]

    for i, original_idx in enumerate(filtered_edge_indices):
        src, dst = clustering_edges[0, i], clustering_edges[1, i]
        score = edge_scores[original_idx]
        src_label = hit_track_labels[src]
        dst_label = hit_track_labels[dst]

        # Track mismatches for diagnostics
        if src_label != dst_label:
            mismatched += 1
            # If both endpoints are in clusters (not isolated), it's a cross-cluster edge
            if src_label != -1 and dst_label != -1:
                cross_cluster += 1

        # Assign edge to source node's cluster (primary assignment)
        if src_label != -1 and src_label in track_edge_groups:
            track_edge_groups[src_label].append((src, dst, score))
        elif dst_label != -1 and dst_label in track_edge_groups:
            # Fallback: assign to destination cluster if source is isolated
            track_edge_groups[dst_label].append((src, dst, score))

    # Print statistics
    print(f"Track cluster statistics:")
    print(f"  GNN score threshold: {score_cut}")
    print(f"  Total tracks: {num_tracks}")
    print(f"  Total clustering edges: {clustering_edges.shape[1]}")
    for track_id in list(unique_tracks)[:5]:
        print(f"    Track {track_id}: {len(track_edge_groups[track_id])} edges")
    if num_tracks > 5:
        print(f"    ... ({num_tracks - 5} more tracks)")
    print(f"  Isolated hits: {num_isolated}")
    if mismatched > 0:
        print(f"  Note: {mismatched} edges with mismatched labels ({cross_cluster} cross-cluster)")
    print()

    # Create Plotly figure
    fig = go.Figure()

    # Add edges for each track cluster
    for track_id in unique_tracks:
        edges = track_edge_groups[track_id]
        if len(edges) == 0:
            continue

        # Build edge coordinates and hover text
        edge_x, edge_y, edge_z = [], [], []
        edge_hover = []
        for src, dst, score in edges:
            edge_x.extend([x[src], x[dst], None])
            edge_y.extend([y[src], y[dst], None])
            edge_z.extend([z_cart[src], z_cart[dst], None])
            # Add hover text for the edge (only for the first point, None for the second)
            edge_hover.extend([
                f'Edge: {src} → {dst}<br>Score: {score:.4f}<br>Track: {track_id}',
                '',  # No hover for the second point
                ''   # No hover for the None separator
            ])

        # Add trace for this track
        color = track_colors[track_id]
        fig.add_trace(go.Scatter3d(
            x=edge_x, y=edge_y, z=edge_z,
            mode='lines',
            line=dict(color=color, width=3),
            opacity=0.8,
            name=f'Track {track_id} ({len(edges)} edges)',
            text=edge_hover,
            hoverinfo='text',
            showlegend=True
        ))

    # Add nodes for each track (colored to match edges)
    for track_id in unique_tracks:
        mask = hit_track_labels == track_id
        if not mask.any():
            continue

        # Create hover text with track ID
        hover_text = [
            f"Node {i}<br>Track: {track_id}<br>r: {r[i]:.1f}<br>φ: {phi[i]:.2f}<br>z: {z[i]:.1f}"
            for i in np.where(mask)[0]
        ]

        color = track_colors[track_id]
        fig.add_trace(go.Scatter3d(
            x=x[mask], y=y[mask], z=z_cart[mask],
            mode='markers',
            marker=dict(
                size=4,
                color=color,
                opacity=0.7,
                line=dict(width=0.5, color='white')
            ),
            name=f'Track {track_id} nodes',
            text=hover_text,
            hoverinfo='text',
            showlegend=False  # Don't duplicate in legend
        ))

    # Add isolated hits (label=-1) in gray
    if num_isolated > 0:
        mask = hit_track_labels == -1
        hover_text = [
            f"Node {i}<br>Isolated<br>r: {r[i]:.1f}<br>φ: {phi[i]:.2f}<br>z: {z[i]:.1f}"
            for i in np.where(mask)[0]
        ]

        fig.add_trace(go.Scatter3d(
            x=x[mask], y=y[mask], z=z_cart[mask],
            mode='markers',
            marker=dict(
                size=3,
                color='gray',
                opacity=0.3,
                line=dict(width=0.5, color='darkgray')
            ),
            name=f'Isolated hits ({num_isolated})',
            text=hover_text,
            hoverinfo='text',
            showlegend=True
        ))

    # Update layout and title
    title_text = f"Track Building Clusters: {dataset_name.capitalize()} Graph #{index}<br>"
    title_text += f"Tracks: {num_tracks} | Clustering Edges (score>{score_cut}): {clustering_edges.shape[1]} | Isolated Hits: {num_isolated}"

    fig.update_layout(**get_standard_scene_layout(title_text, margin_t=60))

    return fig


def main():
    parser = argparse.ArgumentParser(
        description='Visualize proposed track clusters from track building stage',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python create_visual_track_building_proposal.py trainset 1
    # Visualize 1st graph from data/track_building/trainset/

  python create_visual_track_building_proposal.py testset 5 --output custom/path.html
    # Save to custom location

Track Cluster Visualization:
  Each reconstructed track cluster is shown in a different color.
  Edges and nodes are colored by their track assignment (hit_track_labels).
  Isolated hits (not assigned to any track) are shown in gray.

This shows the direct output of the Connected Components clustering algorithm.
        """
    )
    parser.add_argument(
        'dataset',
        type=str,
        choices=['trainset', 'valset', 'testset'],
        help='Dataset to visualize: trainset, valset, or testset'
    )
    parser.add_argument(
        'index',
        type=int,
        help='Index of graph to visualize (1-based, e.g., 1 for first graph)'
    )
    parser.add_argument(
        '--input-dir',
        type=str,
        default=None,
        help='Input directory containing track_building data (default: ../data/track_building)'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        default=None,
        help='Output HTML file path (default: ../data/visuals/track_building/<dataset>/<dataset><N>_track_clusters.html)'
    )

    args = parser.parse_args()

    try:
        visualize_track_building_graph(
            args.dataset,
            args.index,
            args.input_dir,
            args.output
        )
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
