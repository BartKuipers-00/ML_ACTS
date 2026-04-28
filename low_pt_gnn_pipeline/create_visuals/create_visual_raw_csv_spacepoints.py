#!/usr/bin/env python3
"""
create_visual_raw_csv_spacepoints.py <input> [--output OUTPUT] [--color-by COLOR_BY] [--max-points MAX] [--use-simhits]

Creates an interactive 3D HTML visualization of spacepoints.

Two input modes:
  1. PyG graph file (default smeared mode):
       Pass a .pyg file from feature_store/. Uses hit_x/y/z — the correctly
       geometry-transformed smeared coordinates that the model actually trains on.

  2. CSV event prefix (truth simhits mode, --use-simhits):
       Pass a CSV event prefix (e.g. data/csv/event000000000).
       Uses tx/ty/tz exact truth positions from hits.csv.
       Also works for ML_data format (x/y/z columns).

Usage examples:
  python create_visual_raw_csv_spacepoints.py ../data/feature_store/testset/event000011000-graph.pyg
    # Smeared spacepoints — what the model trains on (default)

  python create_visual_raw_csv_spacepoints.py ../data/csv/event000000000 --use-simhits
    # Exact truth hit positions from hits.csv

  python create_visual_raw_csv_spacepoints.py ../data/feature_store/testset/event000011000-graph.pyg --color-by particle
    # Colors hits by particle_id

Dependencies: plotly, pandas, numpy, torch, torch_geometric
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd

try:
    import plotly.graph_objects as go
    import plotly.express as px
except ImportError:
    print("Error: plotly is required. Install with: pip install plotly pandas numpy")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

def read_pyg_graph(pyg_path):
    """
    Load a PyG graph from feature_store and return (df, particle_info).

    df has columns: x, y, z, particle_id, hit_id, t
    particle_info: dict of particle_id -> (pt, eta)
    """
    try:
        import torch
    except ImportError:
        raise RuntimeError("torch is required to read .pyg files")

    g = torch.load(pyg_path, map_location='cpu', weights_only=False)

    df = pd.DataFrame({
        'x': g.hit_x.numpy().astype(float),
        'y': g.hit_y.numpy().astype(float),
        'z': g.hit_z.numpy().astype(float),
        'particle_id': g.hit_particle_id.numpy().astype(int),
        'hit_id': g.hit_id.numpy().astype(int) if hasattr(g, 'hit_id') else np.arange(len(g.hit_x)),
        't': g.hit_t.numpy().astype(float) if hasattr(g, 'hit_t') else np.zeros(len(g.hit_x)),
    })

    particle_info = {}
    if hasattr(g, 'track_particle_id') and hasattr(g, 'track_particle_pt'):
        pids = g.track_particle_id.numpy()
        pts  = g.track_particle_pt.numpy()
        etas = g.track_particle_eta.numpy() if hasattr(g, 'track_particle_eta') else np.zeros(len(pids))
        for pid, pt, eta in zip(pids, pts, etas):
            particle_info[int(pid)] = (float(pt), float(eta))

    return df, particle_info


def detect_format(df_hits):
    if 'x' in df_hits.columns and 'y' in df_hits.columns and 'z' in df_hits.columns:
        return 'ml_data'
    elif 'tx' in df_hits.columns and 'ty' in df_hits.columns and 'tz' in df_hits.columns:
        return 'acts'
    else:
        raise RuntimeError("Could not detect data format. Expected either (x,y,z) or (tx,ty,tz) columns")


def reconstruct_particle_id(df):
    """
    Reconstruct a single particle_id integer from ACTS split barcode columns.

    New ACTS versions write the barcode as 5 separate columns instead of one
    packed 64-bit integer.  Bit layout (ACTS Barcode):
      bits 52-63 : primary vertex  (12 bits)
      bits 36-51 : secondary vertex (16 bits)
      bits 20-35 : particle         (16 bits)
      bits 10-19 : generation       (10 bits)
      bits  0- 9 : sub-particle     (10 bits)

    Uses .astype(object) so pandas delegates to Python's arbitrary-precision
    integers — avoids both the int64 overflow at bit 52 and the missing <<
    operator on numeric Series in older pandas.
    """
    # Plain Python ints: no int64 overflow, << always works, no pandas version issues.
    return [
        (int(pv) << 52) | (int(sv) << 36) | (int(pa) << 20) | (int(ge) << 10) | int(su)
        for pv, sv, pa, ge, su in zip(
            df['particle_id_pv'], df['particle_id_sv'], df['particle_id_part'],
            df['particle_id_gen'], df['particle_id_subpart'],
        )
    ]


def read_hits(hits_path):
    """Read hits CSV (truth simhits), return (df, format_type)."""
    try:
        df = pd.read_csv(hits_path)
        format_type = detect_format(df)

        if format_type == 'acts':
            df = df.rename(columns={'tx': 'x', 'ty': 'y', 'tz': 'z'})
            if 'hit_id' not in df.columns and 'index' in df.columns:
                df['hit_id'] = df['index']
            # New ACTS format: split barcode columns → reconstruct single particle_id
            if 'particle_id' not in df.columns and 'particle_id_part' in df.columns:
                df['particle_id'] = reconstruct_particle_id(df)

        required_cols = ['x', 'y', 'z']
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            raise RuntimeError(f"Missing required columns: {missing}")

        return df, format_type
    except Exception as e:
        raise RuntimeError(f"Error reading hits file: {e}")


def read_truth(truth_path):
    try:
        return pd.read_csv(truth_path)
    except Exception as e:
        print(f"Warning: Could not read truth file: {e}")
        return None


def read_particles(particles_path, particles_simulated_path=None):
    for path in [particles_path, particles_simulated_path]:
        if path is not None and path.exists():
            try:
                df = pd.read_csv(path)
                # New ACTS format: reconstruct single particle_id from split columns
                if 'particle_id' not in df.columns and 'particle_id_part' in df.columns:
                    df['particle_id'] = reconstruct_particle_id(df)
                return df
            except Exception as e:
                print(f"Warning: Could not read particles file {path}: {e}")
    return None


def calculate_pt_eta(px, py, pz):
    pt = np.sqrt(px**2 + py**2)
    p  = np.sqrt(px**2 + py**2 + pz**2)
    eta = np.where(p != pz, 0.5 * np.log((p + pz) / (p - pz)), 0.0)
    return pt, eta


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def create_visualization(df_hits, df_truth=None, df_particles=None,
                         particle_info=None, color_by='particle', max_points=None):
    """
    Create an interactive Plotly figure of the spacepoints.

    particle_info: optional pre-built dict of particle_id -> (pt, eta).
                   When provided, df_particles is ignored for pt/eta lookup.
    """
    # Merge truth data if available (ML_data format)
    if df_truth is not None and 'hit_id' in df_truth.columns and 'particle_id' in df_truth.columns:
        if 'hit_id' in df_hits.columns:
            df_hits = df_hits.merge(df_truth[['hit_id', 'particle_id']], on='hit_id', how='left')
            df_hits['particle_id'] = df_hits['particle_id'].fillna(0)

    # Build particle_info from df_particles if not already supplied
    if particle_info is None:
        particle_info = {}
        if df_particles is not None and 'particle_id' in df_particles.columns:
            for _, row in df_particles.iterrows():
                pid = row['particle_id']
                px_val, py_val, pz_val = row.get('px'), row.get('py'), row.get('pz')
                if px_val is not None and py_val is not None and pz_val is not None:
                    pt, eta = calculate_pt_eta(px_val, py_val, pz_val)
                    particle_info[int(pid)] = (pt, eta)

    x = df_hits['x'].values
    y = df_hits['y'].values
    z = df_hits['z'].values
    num_points = len(x)

    if max_points is not None and num_points > max_points:
        idx = np.random.default_rng(seed=0).choice(num_points, size=max_points, replace=False)
        x = x[idx]; y = y[idx]; z = z[idx]
        df_hits = df_hits.iloc[idx].reset_index(drop=True)
        info_msg = f" (downsampled from {num_points} to {max_points} points)"
    else:
        info_msg = ""

    color_col = None
    title_suffix = ""

    if 'particle_id' in df_hits.columns:
        unique_particles = len(df_hits[df_hits['particle_id'] > 0]['particle_id'].unique())
        title_suffix = f" — {unique_particles} unique particles"

    if color_by == 'particle':
        if 'particle_id' in df_hits.columns:
            color_col = df_hits['particle_id'].astype(str)
        else:
            print("Warning: particle_id not available, using no coloring")
            color_by = 'none'

    if color_by == 'volume':
        if 'volume_id' in df_hits.columns:
            color_col = df_hits['volume_id'].astype(str)
            title_suffix = f" — {df_hits['volume_id'].nunique()} volumes"
        else:
            print("Warning: volume_id not available, using no coloring")
            color_by = 'none'

    if color_by == 'layer':
        if 'layer_id' in df_hits.columns:
            color_col = df_hits['layer_id'].astype(str)
            title_suffix = f" — {df_hits['layer_id'].nunique()} layers"
        else:
            print("Warning: layer_id not available, using no coloring")
            color_by = 'none'

    # Build hover text
    hover_text_final = []
    for i in range(len(x)):
        hit_id = df_hits.iloc[i].get('hit_id', i)
        x_val, y_val, z_val = float(x[i]), float(y[i]), float(z[i])
        r_val = np.sqrt(x_val**2 + y_val**2)
        info = f"Hit {hit_id}<br>"
        info += f"x: {x_val:.2f} mm<br>y: {y_val:.2f} mm<br>z: {z_val:.2f} mm<br>r: {r_val:.2f} mm<br>"
        if 'volume_id' in df_hits.columns:
            info += f"volume_id: {df_hits.iloc[i]['volume_id']}<br>"
        if 'geometry_id' in df_hits.columns:
            info += f"geometry_id: {df_hits.iloc[i]['geometry_id']}<br>"
        if 'layer_id' in df_hits.columns:
            info += f"layer_id: {df_hits.iloc[i]['layer_id']}<br>"
        if 'module_id' in df_hits.columns:
            info += f"module_id: {df_hits.iloc[i]['module_id']}<br>"
        if 'particle_id' in df_hits.columns:
            pid = int(df_hits.iloc[i]['particle_id'])
            if pid != 0:
                info += f"particle_id: {pid}"
                if pid in particle_info:
                    pt, eta = particle_info[pid]
                    info += f" (Pt={pt:.3f} GeV, η={eta:.3f})"
                info += "<br>"
            else:
                info += "particle_id: noise<br>"
        hover_text_final.append(info)

    fig = go.Figure()

    if color_col is not None and color_by == 'particle' and 'particle_id' in df_hits.columns:
        unique_pids = sorted(df_hits['particle_id'].unique())
        colors = px.colors.qualitative.Set3
        for i, pid in enumerate(unique_pids):
            mask = (df_hits['particle_id'] == pid).values
            fig.add_trace(go.Scatter3d(
                x=x[mask], y=y[mask], z=z[mask],
                mode='markers',
                marker=dict(size=3, color='gray' if pid == 0 else colors[i % len(colors)], opacity=0.8),
                text=[hover_text_final[j] for j in range(len(x)) if mask[j]],
                hoverinfo='text',
                name='Noise' if pid == 0 else (
                    f"Particle {int(pid)} (Pt={particle_info[int(pid)][0]:.3f} GeV, η={particle_info[int(pid)][1]:.3f})"
                    if int(pid) in particle_info else f"Particle {int(pid)}"
                ),
                showlegend=True,
            ))
    elif color_col is not None:
        df_plot = pd.DataFrame({'x': x, 'y': y, 'z': z, 'color': color_col})
        temp_fig = px.scatter_3d(df_plot, x='x', y='y', z='z', color='color')
        for trace in temp_fig.data:
            trace.text = hover_text_final
            trace.hoverinfo = 'text'
            fig.add_trace(trace)
    else:
        fig.add_trace(go.Scatter3d(
            x=x, y=y, z=z,
            mode='markers',
            marker=dict(size=3, color='steelblue', opacity=0.8),
            text=hover_text_final,
            hoverinfo='text',
            name=f'Hits ({len(x)})',
        ))

    fig.update_layout(
        title=dict(text=f"Spacepoints ({len(x)} hits){title_suffix}{info_msg}", x=0.5, xanchor='center'),
        scene=dict(
            xaxis=dict(title='x (mm)', backgroundcolor="white", gridcolor="lightgray"),
            yaxis=dict(title='y (mm)', backgroundcolor="white", gridcolor="lightgray"),
            zaxis=dict(title='z (mm)', backgroundcolor="white", gridcolor="lightgray"),
            aspectmode='data',
        ),
        showlegend=True,
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01, bgcolor="rgba(255,255,255,0.8)"),
        margin=dict(l=0, r=0, b=0, t=40),
        hovermode='closest',
    )
    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Create interactive 3D HTML visualization of spacepoints',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python create_visual_raw_csv_spacepoints.py ../data/feature_store/testset/event000011000-graph.pyg
    # Smeared spacepoints from PyG graph — what the model trains on

  python create_visual_raw_csv_spacepoints.py ../data/csv/event000000000 --use-simhits
    # Exact truth hit positions from hits.csv

  python create_visual_raw_csv_spacepoints.py ../data/feature_store/testset/event000011000-graph.pyg --color-by particle
    # Color by particle_id with pt/eta in legend
        """
    )
    parser.add_argument(
        'input',
        type=str,
        help='Path to a .pyg graph file (smeared coords) or a CSV event prefix with --use-simhits'
    )
    parser.add_argument('--output', '-o', type=str, default=None,
                        help='Output HTML file path')
    parser.add_argument('--color-by', type=str, choices=['none', 'particle', 'volume', 'layer'],
                        default='particle', help='Color points by (default: particle)')
    parser.add_argument('--max-points', type=int, default=None,
                        help='Maximum number of points to display')
    parser.add_argument('--use-simhits', action='store_true', default=False,
                        help='Use exact truth hit positions from hits.csv (pass a CSV event prefix as input)')

    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.is_absolute():
        input_path = Path.cwd() / input_path

    df_hits = None
    df_truth = None
    df_particles = None
    particle_info = None
    format_type = None

    if args.use_simhits:
        # CSV truth hits mode
        if input_path.suffix in ('.csv', '.pyg'):
            parent_dir, base_name = input_path.parent, input_path.stem
        else:
            parent_dir, base_name = input_path.parent, input_path.name

        hits_path = parent_dir / f"{base_name}-hits.csv"
        if not hits_path.exists():
            print(f"Error: Hits file not found: {hits_path}")
            sys.exit(1)

        print(f"Loading truth simhits from: {hits_path}")
        try:
            df_hits, format_type = read_hits(hits_path)
        except Exception as e:
            print(f"Error: {e}"); sys.exit(1)
        print(f"Simhits loaded: {len(df_hits)} points (format: {format_type})")

        truth_path             = parent_dir / f"{base_name}-truth.csv"
        particles_path         = parent_dir / f"{base_name}-particles.csv"
        particles_simulated_path = parent_dir / f"{base_name}-particles_simulated.csv"

        if format_type == 'ml_data' and truth_path.exists():
            df_truth = read_truth(truth_path)

        if particles_path.exists() or particles_simulated_path.exists():
            df_particles = read_particles(particles_path, particles_simulated_path)

    elif input_path.suffix == '.pyg':
        # PyG graph mode — correctly smeared coordinates
        if not input_path.exists():
            print(f"Error: PyG file not found: {input_path}")
            sys.exit(1)
        print(f"Loading smeared spacepoints from PyG graph: {input_path}")
        try:
            df_hits, particle_info = read_pyg_graph(input_path)
        except Exception as e:
            print(f"Error: {e}"); sys.exit(1)
        print(f"Graph loaded: {len(df_hits)} hits, {len(particle_info)} particles with pt/eta")
        format_type = 'acts'
        base_name = input_path.stem.replace('-graph', '')
        parent_dir = input_path.parent

    else:
        # CSV prefix without --use-simhits: ML_data fallback or helpful error
        if input_path.suffix in ('', '.csv'):
            parent_dir = input_path.parent
            base_name  = input_path.name if input_path.suffix == '' else input_path.stem
        else:
            parent_dir, base_name = input_path.parent, input_path.stem

        hits_path = parent_dir / f"{base_name}-hits.csv"
        if hits_path.exists():
            try:
                df_test, fmt = read_hits(hits_path)
                if fmt == 'ml_data':
                    # ML_data format has no geometry transform issue
                    print(f"Loading hits from: {hits_path}")
                    df_hits, format_type = df_test, fmt
                    truth_path = parent_dir / f"{base_name}-truth.csv"
                    particles_path = parent_dir / f"{base_name}-particles.csv"
                    if truth_path.exists():
                        df_truth = read_truth(truth_path)
                    if particles_path.exists():
                        df_particles = read_particles(particles_path)
                else:
                    print("Error: For ACTS format, pass a .pyg file from feature_store/ for correct smeared coordinates.")
                    print("       Use --use-simhits to visualize truth hit positions from hits.csv instead.")
                    print(f"       Example: python create_visual_raw_csv_spacepoints.py data/feature_store/testset/event000011000-graph.pyg")
                    sys.exit(1)
            except Exception as e:
                print(f"Error: {e}"); sys.exit(1)
        else:
            print(f"Error: Input not recognised. Pass a .pyg file or a CSV event prefix with --use-simhits.")
            sys.exit(1)

    # Determine output path
    if args.output is not None:
        output_path = Path(args.output)
    else:
        script_dir = Path(__file__).resolve().parent
        visuals_dir = script_dir.parent / 'data' / 'visuals' / 'raw_spacepoints'

        dataset_name = None
        for p in input_path.parents:
            if p.name in ('trainset', 'valset', 'testset'):
                dataset_name = p.name
                break

        if dataset_name:
            output_path = visuals_dir / dataset_name / f"{base_name}-spacepoints.html"
        else:
            output_path = visuals_dir / f"{base_name}-spacepoints.html"

    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"\nCreating visualization (color_by={args.color_by})...")
    try:
        fig = create_visualization(
            df_hits,
            df_truth=df_truth,
            df_particles=df_particles,
            particle_info=particle_info,
            color_by=args.color_by,
            max_points=args.max_points,
        )
    except Exception as e:
        print(f"Error creating visualization: {e}")
        import traceback; traceback.print_exc()
        sys.exit(1)

    print(f"Writing interactive HTML to: {output_path}")
    try:
        fig.write_html(str(output_path), include_plotlyjs='cdn')
    except Exception as e:
        print(f"Error writing HTML file: {e}"); sys.exit(1)

    print("\n✓ Visualization complete!")
    print(f"\nOpen in browser: {output_path.absolute()}")


if __name__ == '__main__':
    main()
