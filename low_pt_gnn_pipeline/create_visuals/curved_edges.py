"""
Helpers to draw edges as short helical arcs instead of straight chords.

Used by visualization scripts to make track-edge lines bow along the fitted
helix circle rather than zigzag from hit to hit. Per-segment helix is fit once
with the existing `fit_helix_to_segment` utility; each edge is sampled along
the source segment's circle, with z linearly interpolated between endpoints.

A linear-ramp correction snaps the arc's endpoints to the actual hit positions,
so consecutive edges join without visible breaks even when the source and
target lie on different segment fits (e.g. cross-segment connectors).
"""

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from low_pt_custom_utils.helix_fitting import fit_helix_to_segment, HelixParams


def fit_helices_per_segment(hit_particle_id, hit_segment_id, x, y, z):
    """Fit one helix per (particle_id, segment_id). Returns dict keyed by (pid, seg)."""
    helices = {}
    if hit_particle_id is None or hit_segment_id is None:
        return helices

    unique_pids = np.unique(hit_particle_id[hit_particle_id > 0])
    for pid in unique_pids:
        p_mask = hit_particle_id == pid
        seg_ids = hit_segment_id[p_mask]
        unique_segs = np.unique(seg_ids[seg_ids > 0])
        for seg in unique_segs:
            seg_mask = p_mask & (hit_segment_id == seg)
            if seg_mask.sum() < 3:
                continue
            params = fit_helix_to_segment(x[seg_mask], y[seg_mask], z[seg_mask])
            if params.fit_quality == "good":
                helices[(int(pid), int(seg))] = params
    return helices


def _arc_xy(helix: HelixParams, x_s, y_s, x_t, y_t, n_points: int):
    """
    Sample (x, y) along the helix circle from source to target via shortest arc.

    A linear-ramp correction snaps the endpoints to the actual hit positions,
    so consecutive edges join without visible breaks even when the source and
    target lie on different fits (e.g. cross-segment connectors).
    """
    xc, yc, R = helix.xc, helix.yc, helix.R
    theta_s = np.arctan2(y_s - yc, x_s - xc)
    theta_t = np.arctan2(y_t - yc, x_t - xc)
    dtheta = (theta_t - theta_s + np.pi) % (2 * np.pi) - np.pi
    ts = np.linspace(0.0, 1.0, n_points)
    theta = theta_s + ts * dtheta
    xs = xc + R * np.cos(theta)
    ys = yc + R * np.sin(theta)

    dx_s, dy_s = x_s - xs[0], y_s - ys[0]
    dx_t, dy_t = x_t - xs[-1], y_t - ys[-1]
    xs = xs + (1.0 - ts) * dx_s + ts * dx_t
    ys = ys + (1.0 - ts) * dy_s + ts * dy_t
    return xs, ys


def build_curved_edge_xyz(
    edges,
    x,
    y,
    z,
    helices_by_seg,
    hit_particle_id,
    hit_segment_id,
    n_points: int = 20,
):
    """
    Like the straight-line _edge_xyz, but sample along the source segment's helix.

    Falls back to a straight line when no good helix is available (segment too
    short, or particle/segment id unknown).
    """
    ex, ey, ez = [], [], []
    n = len(x)
    for j in range(edges.shape[1]):
        s, t = int(edges[0, j]), int(edges[1, j])
        if s >= n or t >= n:
            continue

        helix = None
        if hit_particle_id is not None and hit_segment_id is not None:
            key = (int(hit_particle_id[s]), int(hit_segment_id[s]))
            helix = helices_by_seg.get(key)

        if helix is not None:
            xs, ys = _arc_xy(helix, x[s], y[s], x[t], y[t], n_points)
            zs = np.linspace(z[s], z[t], n_points)
            ex.extend(xs.tolist())
            ey.extend(ys.tolist())
            ez.extend(zs.tolist())
        else:
            ex.extend([x[s], x[t]])
            ey.extend([y[s], y[t]])
            ez.extend([z[s], z[t]])

        ex.append(None)
        ey.append(None)
        ez.append(None)
    return ex, ey, ez
