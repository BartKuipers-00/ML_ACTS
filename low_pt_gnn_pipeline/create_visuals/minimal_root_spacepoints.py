"""
minimal_root_spacepoints.py <path_to_hits.root> [event_id]

Minimal 3D HTML visualization of truth spacepoints from an ACTS ROOT hits file.

Pass the full path to the hits.root file, and optionally an event_id (default 0).
pt/eta are taken from a sibling particles_simulation.root if present.
Writes <stem>-event<N>-spacepoints.html in the current working directory.

Usage (for 5th event):
  python minimal_root_spacepoints.py /abs/path/to/hits.root 5

Dependencies: uproot, plotly, pandas, numpy
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import uproot
import plotly.graph_objects as go
import plotly.express as px


def main():
    if len(sys.argv) not in (2, 3):
        print("Usage: python minimal_root_spacepoints.py <path_to_hits.root> [event_id]")
        sys.exit(1)

    hits_path = Path(sys.argv[1])
    event_id = int(sys.argv[2]) if len(sys.argv) == 3 else 0
    if not hits_path.exists():
        print(f"Error: hits file not found: {hits_path}")
        sys.exit(1)

    print(f"Loading: {hits_path} (event {event_id})")
    tree = uproot.open(hits_path)['hits']
    arr = tree.arrays(['event_id', 'tx', 'ty', 'tz', 'particle_id'], library='np')

    mask = arr['event_id'] == event_id
    if not mask.any():
        print(f"Error: no hits for event_id {event_id}")
        sys.exit(1)

    df = pd.DataFrame({
        'x': arr['tx'][mask],
        'y': arr['ty'][mask],
        'z': arr['tz'][mask],
        'particle_id': arr['particle_id'][mask],
    })

    print(f"Loaded {len(df)} hits, "
          f"{df[df['particle_id'] > 0]['particle_id'].nunique()} particles")

    # Optional: read particles_simulation.root for pt/eta per particle (this event).
    # particles.root has one entry per event; per-particle branches are jagged.
    particle_info = {}
    p_path = hits_path.parent / 'particles_simulation.root'
    if p_path.exists():
        pt_tree = uproot.open(p_path)['particles']
        parr = pt_tree.arrays(['event_id', 'particle_id', 'pt', 'eta'], library='np')
        entry = np.where(parr['event_id'] == event_id)[0]
        if len(entry):
            e = entry[0]
            particle_info = {int(pid): (float(t), float(et))
                             for pid, t, et in zip(parr['particle_id'][e],
                                                   parr['pt'][e], parr['eta'][e])}

    x, y, z = df['x'].values, df['y'].values, df['z'].values

    fig = go.Figure()
    colors = px.colors.qualitative.Set3
    for i, pid in enumerate(sorted(df['particle_id'].unique())):
        m = (df['particle_id'] == pid).values
        if pid == 0:
            name = 'Noise'
        elif int(pid) in particle_info:
            pt, eta = particle_info[int(pid)]
            name = f"Particle {int(pid)} (Pt={pt:.3f} GeV, eta={eta:.2f})"
        else:
            name = f"Particle {int(pid)}"
        fig.add_trace(go.Scatter3d(
            x=x[m], y=y[m], z=z[m],
            mode='markers',
            marker=dict(size=3, color='gray' if pid == 0 else colors[i % len(colors)], opacity=0.8),
            name=name,
        ))

    fig.update_layout(
        title=f"Spacepoints ({len(df)} hits, event {event_id})",
        scene=dict(
            xaxis=dict(title='x (mm)'),
            yaxis=dict(title='y (mm)'),
            zaxis=dict(title='z (mm)'),
            aspectmode='data',
        ),
        margin=dict(l=0, r=0, b=0, t=40),
    )

    output_path = Path.cwd() / f"{hits_path.stem}-event{event_id}-spacepoints.html"
    fig.write_html(str(output_path), include_plotlyjs='cdn')
    print(f"OK Wrote {output_path}")


if __name__ == '__main__':
    main()
