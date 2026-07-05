"""
minimal_root_spacepoints_barcode.py <path_to_hits.root> [event]

Like minimal_root_spacepoints.py, but for the newer ACTS ROOT schema where the
hits tree has no single `particle_id` branch. Instead the particle barcode is
split into 5 columns (barcode_vertex_primary/secondary/particle/generation/
sub_particle). These are packed into one composite particle id for coloring.

pt/eta are taken from a sibling particles_simulation.root if present, matched on
the same 5 barcode components (its `particle_hash` is a hash, NOT the packed
barcode, so we match on the components instead).

Event numbering: the `event` argument is the literal event_id (e.g. 98000). For
convenience, a small value is also accepted as an index into the sorted list of
event_ids present in the file (so `1` = second event).

Usage:
  python minimal_root_spacepoints_barcode.py /abs/path/to/hits.root 98000
  python minimal_root_spacepoints_barcode.py /abs/path/to/hits.root 1   # index

Dependencies: uproot, plotly, pandas, numpy
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import uproot
import plotly.graph_objects as go
import plotly.express as px


# ACTS barcode components -> one composite uint64 id (unique within an event).
# Bit layout is arbitrary here: only needs to be a consistent, collision-free key.
def pack_barcode(vp, vs, particle, generation, sub):
    vp = np.asarray(vp, dtype=np.uint64)
    vs = np.asarray(vs, dtype=np.uint64)
    particle = np.asarray(particle, dtype=np.uint64)
    generation = np.asarray(generation, dtype=np.uint64)
    sub = np.asarray(sub, dtype=np.uint64)
    return ((vp << np.uint64(52)) | (vs << np.uint64(40)) |
            (particle << np.uint64(24)) | (generation << np.uint64(16)) | sub)


def main():
    if len(sys.argv) not in (2, 3):
        print("Usage: python minimal_root_spacepoints_barcode.py <path_to_hits.root> [event]")
        sys.exit(1)

    hits_path = Path(sys.argv[1])
    event_arg = int(sys.argv[2]) if len(sys.argv) == 3 else 0
    if not hits_path.exists():
        print(f"Error: hits file not found: {hits_path}")
        sys.exit(1)

    print(f"Loading: {hits_path} (event arg {event_arg})")
    tree = uproot.open(hits_path)['hits']
    arr = tree.arrays(
        ['event_id', 'tx', 'ty', 'tz',
         'barcode_vertex_primary', 'barcode_vertex_secondary',
         'barcode_particle', 'barcode_generation', 'barcode_sub_particle'],
        library='np')

    # Resolve event: literal event_id if present, else treat as an index.
    event_ids = np.unique(arr['event_id'])
    if event_arg in event_ids:
        event_id = event_arg
    elif 0 <= event_arg < len(event_ids):
        event_id = int(event_ids[event_arg])
        print(f"  (interpreting {event_arg} as index -> event_id {event_id})")
    else:
        print(f"Error: event {event_arg} not found. Available event_ids: "
              f"{event_ids.min()}..{event_ids.max()} ({len(event_ids)} events)")
        sys.exit(1)

    mask = arr['event_id'] == event_id
    particle_id = pack_barcode(
        arr['barcode_vertex_primary'][mask], arr['barcode_vertex_secondary'][mask],
        arr['barcode_particle'][mask], arr['barcode_generation'][mask],
        arr['barcode_sub_particle'][mask])

    df = pd.DataFrame({
        'x': arr['tx'][mask],
        'y': arr['ty'][mask],
        'z': arr['tz'][mask],
        'particle_id': particle_id,
    })

    print(f"Loaded {len(df)} hits, {df['particle_id'].nunique()} particles")

    # Optional: read particles_simulation.root for pt/eta, matched on barcode.
    particle_info = {}
    p_path = hits_path.parent / 'particles_simulation.root'
    if p_path.exists():
        p_tree = uproot.open(p_path)['particles']
        parr = p_tree.arrays(
            ['event_id', 'vertex_primary', 'vertex_secondary', 'particle',
             'generation', 'sub_particle', 'pt', 'eta'], library='np')
        entry = np.where(parr['event_id'] == event_id)[0]
        if len(entry):
            e = entry[0]
            pkey = pack_barcode(parr['vertex_primary'][e], parr['vertex_secondary'][e],
                                parr['particle'][e], parr['generation'][e],
                                parr['sub_particle'][e])
            particle_info = {int(k): (float(t), float(et))
                             for k, t, et in zip(pkey, parr['pt'][e], parr['eta'][e])}

    x, y, z = df['x'].values, df['y'].values, df['z'].values

    fig = go.Figure()
    colors = px.colors.qualitative.Set3
    for i, pid in enumerate(sorted(df['particle_id'].unique())):
        m = (df['particle_id'] == pid).values
        if int(pid) in particle_info:
            pt, eta = particle_info[int(pid)]
            name = f"Particle {int(pid)} (Pt={pt:.3f} GeV, eta={eta:.2f})"
        else:
            name = f"Particle {int(pid)}"
        fig.add_trace(go.Scatter3d(
            x=x[m], y=y[m], z=z[m],
            mode='markers',
            marker=dict(size=3, color=colors[i % len(colors)], opacity=0.8),
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
