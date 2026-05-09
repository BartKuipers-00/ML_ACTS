#!/usr/bin/env python3
"""
create_visual_latent_learning_data_connect.py <dataset> <index> [--output OUTPUT] [--color-by COLOR_BY] [--hide-edges]
create_visual_latent_learning_data_connect.py <graph_file.pyg> [--output OUTPUT] [--color-by COLOR_BY] [--hide-edges]

Same as create_visual_latent_learning_data.py but ALSO draws connector edges
between consecutive segments of the same particle (segment N -> segment N+1).
These cross-segment links are not in `track_edges` (which is built per
(particle, segment)). They are drawn dashed, in the same per-particle color.

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


line_width = 7

# Yellow tones to filter out of plotly palettes (user preference)
_YELLOW_TONES = ('#ffffb3', '#ffed6f')


def _palette():
    return [c for c in px.colors.qualitative.Set3 if c.lower() not in _YELLOW_TONES]


def cylindrical_to_cartesian(r, phi, z):
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return x, y, z


def compute_segment_connectors(hit_particle_id, hit_segment_id, hit_t):
    """
    For each particle, find pairs (last_hit_of_seg_N, first_hit_of_seg_N+1)
    where N+1 is the next segment id present for that particle.

    Returns: int array of shape (2, n_connectors) with hit indices,
    plus int array of length n_connectors with the particle id per pair.
    """
    src = []
    tgt = []
    pids = []

    if hit_particle_id is None or hit_segment_id is None or hit_t is None:
        return np.empty((2, 0), dtype=np.int64), np.empty((0,), dtype=np.int64)

    unique_pids = np.unique(hit_particle_id[hit_particle_id > 0])
    for pid in unique_pids:
        p_mask = np.where(hit_particle_id == pid)[0]
        if p_mask.size < 2:
            continue
        seg_ids = hit_segment_id[p_mask]
        times = hit_t[p_mask]

        unique_segs = np.unique(seg_ids[seg_ids > 0])
        if unique_segs.size < 2:
            continue

        # For each segment, find first and last hit (by time) -- store global idx
        first_idx = {}
        last_idx = {}
        for s in unique_segs:
            s_local = np.where(seg_ids == s)[0]
            if s_local.size == 0:
                continue
            order = np.argsort(times[s_local])
            first_idx[int(s)] = int(p_mask[s_local[order[0]]])
            last_idx[int(s)] = int(p_mask[s_local[order[-1]]])

        sorted_segs = sorted(first_idx.keys())
        for a, b in zip(sorted_segs[:-1], sorted_segs[1:]):
            src.append(last_idx[a])
            tgt.append(first_idx[b])
            pids.append(int(pid))

    if not src:
        return np.empty((2, 0), dtype=np.int64), np.empty((0,), dtype=np.int64)
    return np.array([src, tgt], dtype=np.int64), np.array(pids, dtype=np.int64)


def create_visualization(graph, color_by='none', show_edges=False, max_points=None, dashed_connectors=False, no_grid=False):
    particle_type_map = build_hit_particle_type_map(graph)

    hit_r = graph.hit_r.cpu().numpy()
    hit_phi = graph.hit_phi.cpu().numpy()
    hit_z = graph.hit_z.cpu().numpy()

    x, y, z = cylindrical_to_cartesian(hit_r, hit_phi, hit_z)

    # Optional features
    hit_particle_id = graph.hit_particle_id.cpu().numpy() if 'hit_particle_id' in graph.keys() else None
    hit_region = graph.hit_region.cpu().numpy() if 'hit_region' in graph.keys() else None
    hit_segment_id = graph.hit_segment_id.cpu().numpy() if 'hit_segment_id' in graph.keys() else None
    hit_t = graph.hit_t.cpu().numpy() if 'hit_t' in graph.keys() else None

    # Downsample if needed
    if max_points and len(x) > max_points:
        indices = np.random.choice(len(x), max_points, replace=False)
        x = x[indices]; y = y[indices]; z = z[indices]
        hit_r = hit_r[indices]; hit_phi = hit_phi[indices]; hit_z = hit_z[indices]
        if hit_particle_id is not None: hit_particle_id = hit_particle_id[indices]
        if hit_region is not None: hit_region = hit_region[indices]
        if hit_segment_id is not None: hit_segment_id = hit_segment_id[indices]
        if hit_t is not None: hit_t = hit_t[indices]

    fig = go.Figure()

    # ---- Markers ----
    if color_by == 'particle' and hit_particle_id is not None:
        unique_particles = np.unique(hit_particle_id[hit_particle_id > 0])
        colors = _palette()
        for i, pid in enumerate(unique_particles):
            mask = hit_particle_id == pid
            color = colors[i % len(colors)]
            particle_label = get_particle_label(int(pid), particle_type_map)
            fig.add_trace(go.Scatter3d(
                x=x[mask], y=y[mask], z=z[mask],
                mode='markers', name=particle_label,
                marker=dict(size=2, color=color, opacity=0.7),
                hovertemplate=f'<b>{particle_label}</b><br>r=%{{customdata[0]:.2f}}<br>'
                              'φ=%{customdata[1]:.3f}<br>z=%{customdata[2]:.2f}<br><extra></extra>',
                customdata=np.column_stack([hit_r[mask], hit_phi[mask], hit_z[mask]]),
            ))
        if (hit_particle_id == 0).any():
            mask = hit_particle_id == 0
            fig.add_trace(go.Scatter3d(
                x=x[mask], y=y[mask], z=z[mask],
                mode='markers', name='Noise',
                marker=dict(size=1.5, color='gray', opacity=0.5),
                hovertemplate='<b>Noise</b><br>r=%{customdata[0]:.2f}<br>'
                              'φ=%{customdata[1]:.3f}<br>z=%{customdata[2]:.2f}<br><extra></extra>',
                customdata=np.column_stack([hit_r[mask], hit_phi[mask], hit_z[mask]]),
            ))
    elif color_by == 'region' and hit_region is not None:
        unique_regions = np.unique(hit_region)
        colors = _palette()
        for i, region in enumerate(unique_regions):
            mask = hit_region == region
            color = colors[i % len(colors)]
            fig.add_trace(go.Scatter3d(
                x=x[mask], y=y[mask], z=z[mask],
                mode='markers', name=f'Region {int(region)}',
                marker=dict(size=2, color=color, opacity=0.7),
                hovertemplate=f'<b>Region {int(region)}</b><br>r=%{{customdata[0]:.2f}}<br>'
                              'φ=%{customdata[1]:.3f}<br>z=%{customdata[2]:.2f}<br><extra></extra>',
                customdata=np.column_stack([hit_r[mask], hit_phi[mask], hit_z[mask]]),
            ))
    else:
        fig.add_trace(go.Scatter3d(
            x=x, y=y, z=z, mode='markers', name='Hits',
            marker=dict(size=2, color='blue', opacity=0.7),
            hovertemplate='<b>Hit</b><br>r=%{customdata[0]:.2f}<br>'
                          'φ=%{customdata[1]:.3f}<br>z=%{customdata[2]:.2f}<br><extra></extra>',
            customdata=np.column_stack([hit_r, hit_phi, hit_z]),
        ))

    # ---- Combined truth edges + cross-segment connectors, per particle ----
    if show_edges:
        # Gather track_edges (within-segment) and connectors (cross-segment)
        within_edges = np.empty((2, 0), dtype=np.int64)
        if 'track_edges' in graph.keys() and graph.track_edges.numel() > 0:
            te = graph.track_edges.cpu().numpy()
            if hit_particle_id is not None:
                edge_pids = hit_particle_id[te]
                te = te[:, edge_pids[0] == edge_pids[1]]
            within_edges = te

        connectors, _ = compute_segment_connectors(hit_particle_id, hit_segment_id, hit_t)

        def _edge_xyz(edges):
            ex, ey, ez = [], [], []
            for j in range(edges.shape[1]):
                s, t = edges[0, j], edges[1, j]
                if s < len(x) and t < len(x):
                    ex.extend([x[s], x[t], None])
                    ey.extend([y[s], y[t], None])
                    ez.extend([z[s], z[t], None])
            return ex, ey, ez

        if dashed_connectors:
            # Solid trace for within-segment edges, dashed trace for connectors,
            # grouped per particle so the legend toggles them together.
            if hit_particle_id is None:
                if within_edges.shape[1] > 0:
                    ex, ey, ez = _edge_xyz(within_edges)
                    if ex:
                        fig.add_trace(go.Scatter3d(
                            x=ex, y=ey, z=ez, mode='lines', name='Edges',
                            line=dict(color='red', width=line_width),
                            legendgroup='edges', showlegend=True, hoverinfo='skip',
                        ))
                if connectors.shape[1] > 0:
                    ex, ey, ez = _edge_xyz(connectors)
                    if ex:
                        fig.add_trace(go.Scatter3d(
                            x=ex, y=ey, z=ez, mode='lines', name='Connectors',
                            line=dict(color='red', width=line_width + 2, dash='dash'),
                            legendgroup='edges', showlegend=False, hoverinfo='skip',
                        ))
            else:
                pid_pool = np.concatenate([
                    hit_particle_id[within_edges[0]] if within_edges.shape[1] else np.empty(0, dtype=np.int64),
                    hit_particle_id[connectors[0]] if connectors.shape[1] else np.empty(0, dtype=np.int64),
                ])
                unique_particles = np.unique(pid_pool[pid_pool > 0])
                colors = _palette()
                for i, pid in enumerate(unique_particles):
                    color = colors[i % len(colors)]
                    particle_label = get_particle_label(int(pid), particle_type_map)
                    group = f'edges_{int(pid)}'

                    if within_edges.shape[1]:
                        wmask = hit_particle_id[within_edges[0]] == pid
                        sub = within_edges[:, wmask]
                        ex, ey, ez = _edge_xyz(sub)
                        if ex:
                            fig.add_trace(go.Scatter3d(
                                x=ex, y=ey, z=ez, mode='lines',
                                name=f'Edges ({particle_label})',
                                line=dict(color=color, width=line_width),
                                legendgroup=group, showlegend=True, hoverinfo='skip',
                            ))

                    if connectors.shape[1]:
                        cmask = hit_particle_id[connectors[0]] == pid
                        sub = connectors[:, cmask]
                        ex, ey, ez = _edge_xyz(sub)
                        if ex:
                            fig.add_trace(go.Scatter3d(
                                x=ex, y=ey, z=ez, mode='lines',
                                name=f'Connectors ({particle_label})',
                                line=dict(color=color, width=line_width + 1, dash='dash'),
                                legendgroup=group, showlegend=False, hoverinfo='skip',
                            ))
        else:
            # Original behavior: merge into a single solid trace per particle
            all_edges = np.concatenate([within_edges, connectors], axis=1) if connectors.shape[1] else within_edges

            if all_edges.shape[1] > 0:
                if hit_particle_id is None:
                    ex, ey, ez = _edge_xyz(all_edges)
                    if ex:
                        fig.add_trace(go.Scatter3d(
                            x=ex, y=ey, z=ez, mode='lines', name='Edges',
                            line=dict(color='red', width=line_width),
                            showlegend=True, hoverinfo='skip',
                        ))
                else:
                    edge_pids_src = hit_particle_id[all_edges[0]]
                    unique_particles = np.unique(edge_pids_src[edge_pids_src > 0])
                    colors = _palette()
                    for i, pid in enumerate(unique_particles):
                        pmask = edge_pids_src == pid
                        sub = all_edges[:, pmask]
                        ex, ey, ez = _edge_xyz(sub)
                        if not ex:
                            continue
                        particle_label = get_particle_label(int(pid), particle_type_map)
                        fig.add_trace(go.Scatter3d(
                            x=ex, y=ey, z=ez, mode='lines',
                            name=f'Edges ({particle_label})',
                            line=dict(color=colors[i % len(colors)], width=line_width),
                            showlegend=True, hoverinfo='skip',
                        ))

    # ---- Layout ----
    event_id = graph.event_id[0] if hasattr(graph, "event_id") else "Unknown"
    if no_grid:
        hidden_axis = dict(
            title='', showbackground=False, showgrid=False, zeroline=False,
            showticklabels=False, showline=False, visible=False,
        )
        scene = dict(
            xaxis=hidden_axis, yaxis=hidden_axis, zaxis=hidden_axis,
            aspectmode='data',
        )
    else:
        scene = dict(
            xaxis=dict(title='x (mm)', backgroundcolor="white", gridcolor="lightgray"),
            yaxis=dict(title='y (mm)', backgroundcolor="white", gridcolor="lightgray"),
            zaxis=dict(title='z (mm)', backgroundcolor="white", gridcolor="lightgray"),
            aspectmode='data',
        )
    fig.update_layout(
        title=dict(text=f'PyG Graph + Segment Connectors - Event {event_id}',
                   x=0.5, xanchor='center'),
        scene=scene,
        showlegend=True,
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01,
                    bgcolor="rgba(255, 255, 255, 0.8)"),
        margin=dict(l=0, r=0, b=0, t=40),
        hovermode='closest',
    )
    return fig


def main():
    parser = argparse.ArgumentParser(
        description='Create interactive 3D HTML visualization of PyG graph data with cross-segment connectors',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('input', type=str,
                        help='Either: (1) "<dataset> <index>" like "testset 1", or (2) full path to .pyg file')
    parser.add_argument('index', type=str, nargs='?', default=None,
                        help='Index number (if first argument is dataset name)')
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='Output HTML file path (default: data/visuals/latent_data_connect/<dataset>/<filename>.html)')
    parser.add_argument('--color-by', choices=['none', 'particle', 'region'], default='none',
                        help='Color points by: none, particle, or region (default: none)')
    parser.add_argument('--hide-edges', action='store_false', dest='show_edges', default=True,
                        help='Hide truth edges and segment connectors (default: shown)')
    parser.add_argument('--dashed', action='store_true', default=False,
                        help='Draw cross-segment connector edges as dashed (same color, toggles with the particle)')
    parser.add_argument('--no_grid', '--no-grid', action='store_true', default=False,
                        help='Hide axes, labels, and gridlines — show only the points and edges')
    parser.add_argument('--max-points', type=int, default=None,
                        help='Maximum number of points to display')
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    feature_store_dir = script_dir.parent / 'data' / 'feature_store'

    if args.index is not None:
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

        graph_files = sorted([f for f in dataset_dir.glob('*.pyg')])
        if not graph_files:
            print(f"Error: No graph files found in {dataset_dir}")
            sys.exit(1)
        if index < 1 or index > len(graph_files):
            print(f"Error: Index {index} out of range. Available: 1-{len(graph_files)}")
            sys.exit(1)

        graph_path = graph_files[index - 1]
        print(f"Selected: {graph_path.name} (index {index} of {len(graph_files)} in {dataset_name})")
    else:
        graph_path = Path(args.input)
        if not graph_path.is_absolute():
            graph_path = script_dir / graph_path

    if not graph_path.exists():
        print(f"Error: Graph file not found: {graph_path}")
        sys.exit(1)

    print(f"Loading PyG graph from: {graph_path}")
    try:
        graph = torch.load(graph_path, weights_only=False, map_location='cpu')
    except Exception as e:
        print(f"Error loading PyG file: {e}")
        import traceback; traceback.print_exc()
        sys.exit(1)

    required_features = ['hit_r', 'hit_phi', 'hit_z']
    missing = [f for f in required_features if f not in graph.keys()]
    if missing:
        print(f"Error: Missing required features: {missing}")
        print(f"Available features: {list(graph.keys())}")
        sys.exit(1)

    print(f"Graph loaded: {len(graph.hit_r)} hits")
    if 'hit_segment_id' not in graph.keys():
        print("  Warning: hit_segment_id not present — no segment connectors will be drawn.")
    if 'hit_t' not in graph.keys():
        print("  Warning: hit_t not present — segment connectors require time to order hits.")

    if args.output is not None:
        output_path = Path(args.output)
    else:
        visuals_dir = script_dir.parent / 'data' / 'visuals' / 'latent_data_connect'
        try:
            relative_path = graph_path.relative_to(feature_store_dir)
            output_path = visuals_dir / relative_path.with_suffix('.html')
        except ValueError:
            dataset_name = None
            for parent in graph_path.parents:
                if parent.name in ['trainset', 'valset', 'testset']:
                    dataset_name = parent.name
                    break
            if dataset_name:
                output_path = visuals_dir / dataset_name / graph_path.name.replace('.pyg', '.html')
            else:
                output_path = visuals_dir / graph_path.name.replace('.pyg', '.html')
        output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"\nCreating visualization (color_by={args.color_by}, show_edges={args.show_edges}, dashed={args.dashed}, no_grid={args.no_grid})...")
    fig = create_visualization(graph, color_by=args.color_by,
                               show_edges=args.show_edges, max_points=args.max_points,
                               dashed_connectors=args.dashed, no_grid=args.no_grid)

    print(f"Writing interactive HTML to: {output_path}")
    fig.write_html(str(output_path), include_plotlyjs='cdn')
    print("\n✓ Visualization complete!")
    print(f"\nOpen in browser: {output_path.absolute()}")


if __name__ == '__main__':
    main()
