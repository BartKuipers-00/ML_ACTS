#!/usr/bin/env python3
"""
minimal_csv_spacepoints.py <path_to_hits.csv>

Minimal 3D HTML visualization of truth spacepoints from an ACTS CSV event.

Writes <event>-spacepoints.html in the current working directory.

Usage:

  python minimal_csv_spacepoints.py /abs/path/to/event000000000-hits.csv

Dependencies: plotly, pandas, numpy
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px


def reconstruct_particle_id(df):
    """Rebuild a single 64-bit ACTS barcode from the 5 split columns."""
    return [
        (int(pv) << 52) | (int(sv) << 36) | (int(pa) << 20) | (int(ge) << 10) | int(su)
        for pv, sv, pa, ge, su in zip(
            df['particle_id_pv'], df['particle_id_sv'], df['particle_id_part'],
            df['particle_id_gen'], df['particle_id_subpart'],
        )
    ]


def main():
    if len(sys.argv) != 2:
        print("Usage: python minimal_csv_spacepoints.py <path_to_hits.csv>")
        sys.exit(1)

    hits_path = Path(sys.argv[1])
    base_name = hits_path.name.replace('-hits.csv', '')
    if not hits_path.exists():
        print(f"Error: hits file not found: {hits_path}")
        sys.exit(1)

    print(f"Loading: {hits_path}")
    df = pd.read_csv(hits_path)

    # Truth positions: tx/ty/tz -> x/y/z
    df = df.rename(columns={'tx': 'x', 'ty': 'y', 'tz': 'z'})
    if not {'x', 'y', 'z'}.issubset(df.columns):
        print("Error: expected tx/ty/tz columns in hits CSV")
        sys.exit(1)

    # Reconstruct particle_id from split barcode columns if needed
    if 'particle_id' not in df.columns and 'particle_id_part' in df.columns:
        df['particle_id'] = reconstruct_particle_id(df)
    if 'particle_id' not in df.columns:
        df['particle_id'] = 0

    print(f"Loaded {len(df)} hits, "
          f"{df[df['particle_id'] > 0]['particle_id'].nunique()} particles")

    x, y, z = df['x'].values, df['y'].values, df['z'].values

    fig = go.Figure()
    colors = px.colors.qualitative.Set3
    for i, pid in enumerate(sorted(df['particle_id'].unique())):
        mask = (df['particle_id'] == pid).values
        fig.add_trace(go.Scatter3d(
            x=x[mask], y=y[mask], z=z[mask],
            mode='markers',
            marker=dict(size=3, color='gray' if pid == 0 else colors[i % len(colors)], opacity=0.8),
            name='Noise' if pid == 0 else f"Particle {int(pid)}",
        ))

    fig.update_layout(
        title=f"Spacepoints ({len(df)} hits)",
        scene=dict(
            xaxis=dict(title='x (mm)'),
            yaxis=dict(title='y (mm)'),
            zaxis=dict(title='z (mm)'),
            aspectmode='data',
        ),
        margin=dict(l=0, r=0, b=0, t=40),
    )

    output_path = Path.cwd() / f"{base_name}-spacepoints.html"
    fig.write_html(str(output_path), include_plotlyjs='cdn')
    print(f"✓ Wrote {output_path}")


if __name__ == '__main__':
    main()
