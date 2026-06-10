"""
Segment matching logic for helix-based track building.

Extracts segments from GNN output (via CC clustering or ground truth),
fits helices to each segment, and matches segments from the same particle
based on helix parameter compatibility (circle center, radius, pitch).

The key physics insight: segments from the same looping particle share
the same circle center and similar radius in the x-y plane.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional

import numpy as np
import torch
import scipy.sparse as sps
from torch_geometric.utils import remove_isolated_nodes, to_scipy_sparse_matrix

from low_pt_custom_utils.helix_fitting import HelixParams, fit_helix_to_segment


@dataclass
class SegmentInfo:
    """A segment: a group of hits belonging to one radial arc of a particle."""
    hits: List[int]                    # Node indices (time-ordered)
    helix: Optional[HelixParams]       # Fitted helix parameters (None if not fitted)
    inner_endpoint: Tuple[float, float, float]   # (x, y, z) of innermost-r hit
    outer_endpoint: Tuple[float, float, float]   # (x, y, z) of outermost-r hit
    inner_r: float                     # Radial coordinate of innermost hit
    outer_r: float                     # Radial coordinate of outermost hit
    cluster_id: int = -1               # CC cluster or GT segment identifier


# ─── Segment Extraction ────────────────────────────────────────────────────


def extract_segments_from_ground_truth(graph) -> List[SegmentInfo]:
    """
    Extract segments using ground truth hit_segment_id and hit_particle_id.

    Each unique (particle_id, segment_id) pair with particle_id > 0 forms a segment.
    Hits are ordered by radial position (hit_r) outward.

    Args:
        graph: PyG Data object with hit_segment_id, hit_particle_id,
               hit_x, hit_y, hit_z, hit_r.

    Returns:
        List of SegmentInfo objects.
    """

    # particle_ids = np.asarray(graph.hit_particle_id.cpu().numpy(), dtype=np.int64)
    # segment_ids = np.asarray(graph.hit_segment_id.cpu().numpy(), dtype=np.int64)
    # times = np.asarray(graph.hit_t.cpu().numpy(), dtype=np.float64)
    # x = np.asarray(graph.hit_x.cpu().numpy(), dtype=np.float64)
    # y = np.asarray(graph.hit_y.cpu().numpy(), dtype=np.float64)
    # z = np.asarray(graph.hit_z.cpu().numpy(), dtype=np.float64)
    # r = np.asarray(graph.hit_r.cpu().numpy(), dtype=np.float64)

    segments = []

    particle_ids = np.asarray(graph.hit_particle_id.cpu().numpy(), dtype=np.int64)
    segment_ids = np.asarray(graph.hit_segment_id.cpu().numpy(), dtype=np.int64)
    x = np.asarray(graph.hit_x.cpu().numpy(), dtype=np.float64)
    y = np.asarray(graph.hit_y.cpu().numpy(), dtype=np.float64)
    z = np.asarray(graph.hit_z.cpu().numpy(), dtype=np.float64)
    r = np.asarray(graph.hit_r.cpu().numpy(), dtype=np.float64)

    segments = []

    # Find unique (particle_id, segment_id) pairs, skip noise (pid=0)
    signal_mask = particle_ids > 0
    signal_pids = particle_ids[signal_mask]
    signal_sids = segment_ids[signal_mask]
    signal_indices = np.where(signal_mask)[0]

    # Group by (particle_id, segment_id)
    composite_ids = signal_pids * 10000 + signal_sids
    unique_composites = np.unique(composite_ids)

    for cid in unique_composites.tolist():
        mask = composite_ids == cid
        hit_indices = signal_indices[mask]

        # Order by radial position outward
        r_order = np.argsort(r[hit_indices])
        hit_indices = hit_indices[r_order]

        seg = _build_segment_info(hit_indices.tolist(), x, y, z, r, cluster_id=int(cid))
        segments.append(seg)

    return segments


def extract_segments_from_cc(graph, score_cut: float) -> List[SegmentInfo]:
    """
    Extract segments via Connected Components clustering on GNN edge scores.

    Each CC cluster is treated as one segment. Single isolated nodes are skipped.
    Hits within each segment are ordered by radial position (hit_r) outward.

    Args:
        graph: PyG Data object with edge_index, edge_scores,
               hit_x, hit_y, hit_z, hit_r.
        score_cut: Threshold for edge score filtering.

    Returns:
        List of SegmentInfo objects.
    """
    edge_mask = graph.edge_scores > score_cut
    edges = graph.edge_index[:, edge_mask]

    num_nodes = graph.hit_x.size(0)

    edges_clean, _, node_mask = remove_isolated_nodes(edges, num_nodes=num_nodes)
    n_connected = node_mask.sum().item()

    if n_connected == 0:
        return []

    sparse_edges = to_scipy_sparse_matrix(edges_clean, num_nodes=n_connected)
    n_components, candidate_labels = sps.csgraph.connected_components(
        sparse_edges, directed=False, return_labels=True
    )

    # Map back to original node indices
    labels = (torch.ones(num_nodes, dtype=torch.long) * -1)
    labels[node_mask] = torch.tensor(candidate_labels, dtype=torch.long)
    labels = labels.numpy()

    x = np.asarray(graph.hit_x.cpu().numpy(), dtype=np.float64)
    y = np.asarray(graph.hit_y.cpu().numpy(), dtype=np.float64)
    z = np.asarray(graph.hit_z.cpu().numpy(), dtype=np.float64)
    r = np.asarray(graph.hit_r.cpu().numpy(), dtype=np.float64)

    segments = []
    for cluster_id in range(n_components):
        hit_indices = np.where(labels == cluster_id)[0]

        if len(hit_indices) < 1:
            continue

        # Order by radial position outward
        r_order = np.argsort(r[hit_indices])
        hit_indices = hit_indices[r_order]

        seg = _build_segment_info(hit_indices.tolist(), x, y, z, r, cluster_id=cluster_id)
        segments.append(seg)

    return segments


def _build_segment_info(
    hit_indices: List[int],
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    r: np.ndarray,
    cluster_id: int = -1,
) -> SegmentInfo:
    """Build a SegmentInfo from hit indices and coordinate arrays."""
    hits_r = r[hit_indices]
    inner_idx = hit_indices[np.argmin(hits_r)]
    outer_idx = hit_indices[np.argmax(hits_r)]

    return SegmentInfo(
        hits=hit_indices,
        helix=None,  # Will be fitted in a separate step
        inner_endpoint=(float(x[inner_idx]), float(y[inner_idx]), float(z[inner_idx])),
        outer_endpoint=(float(x[outer_idx]), float(y[outer_idx]), float(z[outer_idx])),
        inner_r=float(r[inner_idx]),
        outer_r=float(r[outer_idx]),
        cluster_id=cluster_id,
    )


# ─── Helix Fitting ─────────────────────────────────────────────────────────


def fit_helices_to_segments(
    segments: List[SegmentInfo],
    graph,
    B_field: float = 2.0,
    outlier_rejection: bool = False,
) -> List[SegmentInfo]:
    """
    Fit helix parameters to each segment using ALL hits.

    Args:
        segments: List of SegmentInfo objects.
        graph: PyG Data object with hit_x, hit_y, hit_z.
        B_field: Magnetic field strength (Tesla).
        outlier_rejection: Whether to apply outlier rejection.

    Returns:
        Same list with helix field populated.

    Batched (bincount sums + batched 3x3 solve); n<3 / singular / outlier-rejection fall back
    to the per-segment routine.
    """
    x = graph.hit_x.cpu().numpy().astype(np.float64, copy=False)
    y = graph.hit_y.cpu().numpy().astype(np.float64, copy=False)
    z = graph.hit_z.cpu().numpy().astype(np.float64, copy=False)
    N = len(segments)
    if N == 0:
        return segments

    if outlier_rejection:   # iterative path — keep the exact per-segment routine
        for seg in segments:
            seg.helix = fit_helix_to_segment(x[seg.hits], y[seg.hits], z[seg.hits],
                                             B_field=B_field, outlier_rejection=True)
        return segments

    nhits = np.fromiter((len(s.hits) for s in segments), dtype=np.int64, count=N)
    flat = np.concatenate([np.asarray(s.hits, dtype=np.int64) for s in segments])
    seg_ids = np.repeat(np.arange(N), nhits)
    xs, ys, zs = x[flat], y[flat], z[flat]
    nf = nhits.astype(np.float64)
    bc = lambda w: np.bincount(seg_ids, weights=w, minlength=N)

    # --- batched Kasa circle fit via per-segment normal equations  A=[2x,2y,1], b=x^2+y^2 ---
    bb = xs * xs + ys * ys
    Sx, Sy = bc(xs), bc(ys)
    Sxx, Syy, Sxy = bc(xs * xs), bc(ys * ys), bc(xs * ys)
    Sb, Sxb, Syb = bc(bb), bc(xs * bb), bc(ys * bb)
    AtA = np.zeros((N, 3, 3))
    AtA[:, 0, 0] = 4 * Sxx; AtA[:, 0, 1] = 4 * Sxy; AtA[:, 0, 2] = 2 * Sx
    AtA[:, 1, 0] = 4 * Sxy; AtA[:, 1, 1] = 4 * Syy; AtA[:, 1, 2] = 2 * Sy
    AtA[:, 2, 0] = 2 * Sx;  AtA[:, 2, 1] = 2 * Sy;  AtA[:, 2, 2] = nf
    Atb = np.stack([2 * Sxb, 2 * Syb, Sb], axis=1)

    xc = np.zeros(N); yc = np.zeros(N); R = np.full(N, np.nan)
    ok = nhits >= 3
    gi = np.where(ok)[0]
    if gi.size:
        sol = np.full((gi.size, 3), np.nan)
        try:
            sol = np.linalg.solve(AtA[gi], Atb[gi])
        except np.linalg.LinAlgError:               # a singular segment in the batch
            for j, k in enumerate(gi):
                try: sol[j] = np.linalg.solve(AtA[k], Atb[k])
                except np.linalg.LinAlgError: sol[j] = np.nan
        xc[gi], yc[gi] = sol[:, 0], sol[:, 1]
        R2 = sol[:, 2] + sol[:, 0] ** 2 + sol[:, 1] ** 2
        R[gi] = np.sqrt(np.where(R2 > 0, R2, np.nan))
    valid = ok & np.isfinite(R) & (R > 0)
    Rcap = np.minimum(R, 10000.0)                    # cap (used for arc/pT/stored R; residuals use raw R)

    # residual_rms per segment (uses uncapped R, matching the scalar routine)
    dist = np.sqrt((xs - xc[seg_ids]) ** 2 + (ys - yc[seg_ids]) ** 2)
    with np.errstate(invalid="ignore"):
        residual_rms = np.sqrt(bc((dist - R[seg_ids]) ** 2) / nf)

    # arc length per hit (cumulative within segment), then closed-form pitch regression
    theta = np.arctan2(ys - yc[seg_ids], xs - xc[seg_ids])
    starts = np.zeros(N, dtype=np.int64); starts[1:] = np.cumsum(nhits)[:-1]
    dtheta = np.zeros(flat.size); dtheta[1:] = theta[1:] - theta[:-1]
    dtheta = (dtheta + np.pi) % (2 * np.pi) - np.pi
    dtheta[starts] = 0.0                             # no step at a segment's first hit
    contrib = np.abs(dtheta) * Rcap[seg_ids]
    contrib = np.where(valid[seg_ids], contrib, 0.0)  # don't let an invalid segment poison the cumsum
    cum = np.cumsum(contrib)
    offset = np.where(starts > 0, cum[np.maximum(starts - 1, 0)], 0.0)
    arc = cum - offset[seg_ids]
    arc_span = arc[starts + nhits - 1]              # arc is monotonic within a segment, first = 0

    Ss, Sz, Sss, Ssz = bc(arc), bc(zs), bc(arc * arc), bc(arc * zs)
    denom = nf * Sss - Ss * Ss
    with np.errstate(invalid="ignore", divide="ignore"):
        pitch_arr = (nf * Ssz - Ss * Sz) / denom
        z0_lin = (Sz - pitch_arr * Ss) / nf
    z_mean = Sz / nf
    pT = 0.3 * B_field * Rcap / 1000.0
    phi = np.arctan2(yc, xc)

    # Build HelixParams — math is batched; only the cheap object construction loops.
    for i, seg in enumerate(segments):
        if not valid[i]:                            # n<3 / singular / R^2<=0 -> exact scalar routine
            seg.helix = fit_helix_to_segment(x[seg.hits], y[seg.hits], z[seg.hits],
                                             B_field=B_field, outlier_rejection=False)
            continue
        has_pitch = arc_span[i] > 1e-6
        seg.helix = HelixParams(
            xc=float(xc[i]), yc=float(yc[i]), R=float(Rcap[i]),
            pitch=float(pitch_arr[i]) if has_pitch else None,
            z0=float(z0_lin[i]) if has_pitch else float(z_mean[i]),
            pT=float(pT[i]), phi_center=float(phi[i]),
            fit_quality="good", residual_rms=float(residual_rms[i]), nhits=int(nhits[i]),
        )
    return segments


# ─── Matching Score ────────────────────────────────────────────────────────


def compute_matching_score(
    seg_a: SegmentInfo,
    seg_b: SegmentInfo,
    config: dict,
) -> Tuple[float, bool]:
    """
    Compute helix-based matching score between two segments.

    Uses circle center distance, radius ratio, and pitch difference.

    Args:
        seg_a, seg_b: Segments with fitted helix parameters.
        config: Dict with matching parameters:
            - max_center_distance (mm)
            - min_R_ratio
            - sigma_center (mm)
            - sigma_R
            - sigma_pitch
            - weight_center, weight_R, weight_pitch

    Returns:
        (score, is_compatible): score in [0, 1], is_compatible is bool.
    """
    helix_a = seg_a.helix
    helix_b = seg_b.helix

    # Both must have valid circle fits
    if helix_a is None or helix_b is None:
        return 0.0, False
    if helix_a.fit_quality != "good" or helix_b.fit_quality != "good":
        return 0.0, False

    # Config parameters
    max_center_dist = config.get("max_center_distance", 100.0)
    min_R_ratio = config.get("min_R_ratio", 0.5)
    sigma_center = config.get("sigma_center", 30.0)
    sigma_R = config.get("sigma_R", 0.1)
    sigma_pitch = config.get("sigma_pitch", 0.01)
    w_center = config.get("weight_center", 0.5)
    w_R = config.get("weight_R", 0.3)
    w_pitch = config.get("weight_pitch", 0.2)

    # Center distance
    center_dist = np.sqrt((helix_a.xc - helix_b.xc)**2 + (helix_a.yc - helix_b.yc)**2)

    # Hard cut on center distance
    if center_dist > max_center_dist:
        return 0.0, False

    # Radius ratio (always <= 1)
    R_max = max(helix_a.R, helix_b.R)
    R_min = min(helix_a.R, helix_b.R)
    if R_max < 1e-6:
        return 0.0, False
    R_ratio = R_min / R_max

    # Hard cut on radius ratio
    if R_ratio < min_R_ratio:
        return 0.0, False

    # Soft scores (Gaussian-weighted)
    score_center = np.exp(-center_dist**2 / (2 * sigma_center**2))
    score_R = np.exp(-(R_ratio - 1.0)**2 / (2 * sigma_R**2))

    # Pitch score (if both have valid pitch fits).
    # Compare |pitch| magnitudes: the wrangler walk reverses arc2's hit order
    # (radially outward), making its dz/ds negative while arc1's is positive.
    # Both arcs of the same helix have the same |pitch|, so magnitude comparison
    # is correct regardless of extraction method.
    if helix_a.pitch is not None and helix_b.pitch is not None:
        dpitch = abs(abs(helix_a.pitch) - abs(helix_b.pitch))
        score_pitch = np.exp(-dpitch**2 / (2 * sigma_pitch**2))
    else:
        score_pitch = 1.0  # Neutral if pitch unavailable
        w_pitch = 0.0  # Don't weight pitch if unavailable

    # Normalize weights
    total_weight = w_center + w_R + w_pitch
    if total_weight < 1e-6:
        return 0.0, False

    score = (w_center * score_center + w_R * score_R + w_pitch * score_pitch) / total_weight

    return float(score), True


# ─── Greedy Matching ───────────────────────────────────────────────────────


def match_segments(
    segments: List[SegmentInfo],
    config: dict,
    timing: dict = None,
) -> Tuple[List[List[SegmentInfo]], List[SegmentInfo]]:
    """
    Greedy matching of segments based on helix parameter compatibility.

    First, segments reaching the outer detector radius are marked as complete
    (straight-through) tracks and excluded from matching. Then, pairwise
    matching scores are computed between remaining segments and assigned greedily.

    Args:
        segments: List of SegmentInfo objects with fitted helices.
        config: Dict with matching parameters (see compute_matching_score).
            Also uses 'outer_r_threshold' (default 1000mm) to identify
            complete tracks that don't need matching.

    Returns:
        matched_tracks: List of lists of SegmentInfo (each list = one track).
        unmatched: List of SegmentInfo objects that were not matched.
    """
    from time import perf_counter
    t_score0 = perf_counter()
    outer_r_threshold = config.get("outer_r_threshold", 1000.0)

    # Step 1: Separate complete (straight-through) segments from matching candidates.
    # Segments that reach the outer detector boundary are already complete tracks.
    complete = set()  # Indices of segments that are complete tracks
    for i, seg in enumerate(segments):
        if seg.outer_r >= outer_r_threshold:
            complete.add(i)

    # Step 2: Among non-complete segments, separate by fit quality
    fittable = [
        (i, seg) for i, seg in enumerate(segments)
        if i not in complete and seg.helix and seg.helix.fit_quality == "good"
    ]

    # Step 3: vectorised pairwise scoring (same candidates as compute_matching_score, no loop).
    cand_a, cand_b = [], []
    if len(fittable) >= 2:
        idxs = np.array([i for i, _ in fittable])
        hel = [s.helix for _, s in fittable]
        xc = np.array([h.xc for h in hel]); yc = np.array([h.yc for h in hel])
        R = np.array([h.R for h in hel])
        pit = np.array([abs(h.pitch) if h.pitch is not None else np.nan for h in hel])

        max_cd = config.get("max_center_distance", 100.0)
        min_Rr = config.get("min_R_ratio", 0.5)
        s_c = config.get("sigma_center", 30.0); s_R = config.get("sigma_R", 0.1); s_p = config.get("sigma_pitch", 0.01)
        w_c = config.get("weight_center", 0.5); w_R = config.get("weight_R", 0.3); w_p0 = config.get("weight_pitch", 0.2)

        ii, jj = np.triu_indices(len(fittable), k=1)        # i<j, row-major (= original loop order)
        cd = np.hypot(xc[ii] - xc[jj], yc[ii] - yc[jj])
        Rmx = np.maximum(R[ii], R[jj]); Rmn = np.minimum(R[ii], R[jj])
        with np.errstate(divide="ignore", invalid="ignore"):
            Rr = np.where(Rmx >= 1e-6, Rmn / Rmx, 0.0)
        ks = np.nonzero((cd <= max_cd) & (Rmx >= 1e-6) & (Rr >= min_Rr))[0]   # hard cuts first

        # soft Gaussian score on survivors only (avoid exp over all N^2 pairs)
        cdk, Rrk = cd[ks], Rr[ks]; pik, pjk = pit[ii[ks]], pit[jj[ks]]
        has_p = ~np.isnan(pik) & ~np.isnan(pjk)
        s_pitch = np.where(has_p, np.exp(-np.abs(pik - pjk) ** 2 / (2 * s_p ** 2)), 1.0)
        w_p = np.where(has_p, w_p0, 0.0); tw = w_c + w_R + w_p
        score = (w_c * np.exp(-cdk ** 2 / (2 * s_c ** 2))
                 + w_R * np.exp(-(Rrk - 1.0) ** 2 / (2 * s_R ** 2)) + w_p * s_pitch) / tw

        keep2 = score > 0
        sel = ks[keep2]
        order = sel[np.argsort(-score[keep2], kind="stable")]           # stable -> same tie-break
        cand_a = idxs[ii[order]].tolist(); cand_b = idxs[jj[order]].tolist()

    # Greedy matching
    t_score = perf_counter() - t_score0
    t_g0 = perf_counter()
    matched = set()
    matched_tracks = []

    for idx_a, idx_b in zip(cand_a, cand_b):
        if idx_a in matched or idx_b in matched:
            continue
        matched_tracks.append([segments[idx_a], segments[idx_b]])
        matched.add(idx_a)
        matched.add(idx_b)

    # Collect unmatched segments (including short ones)
    unmatched = []
    for i, seg in enumerate(segments):
        if i not in matched:
            unmatched.append(seg)

    if timing is not None:
        timing["score"] = t_score
        timing["greedy"] = perf_counter() - t_g0
    return matched_tracks, unmatched


# ─── Geometric Pair Cut ────────────────────────────────────────────────────


def precompute_high_r_hits(
    segments: List["SegmentInfo"],
    graph,
    r_high: float,
) -> List[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """
    Pre-sort each segment's hits by r and retain only those above r_high.

    Call once per event before the pair loop so each segment is processed O(1) times
    rather than O(n) times (once per pair).

    Returns:
        List of (rs, xs, ys, zs) arrays, one entry per segment.
        Arrays are empty if the segment has no hits above r_high.
    """
    hr = graph.hit_r.numpy()
    hx = graph.hit_x.numpy()
    hy = graph.hit_y.numpy()
    hz = graph.hit_z.numpy()

    result = []
    for seg in segments:
        idx   = np.array(seg.hits)
        order = np.argsort(hr[idx])
        rs    = hr[idx][order]
        mask  = rs > r_high
        result.append((rs[mask], hx[idx][order][mask],
                        hy[idx][order][mask], hz[idx][order][mask]))
    return result


def passes_geometric_cut(
    hi: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    hj: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    r_tol: float,
    cut_mm: float,
) -> bool:
    """
    Return True if the segment pair is geometrically compatible (should be kept).

    Computes the minimum 3D distance between nearest-r matched hits at r > r_high
    (pre-filtered by ``precompute_high_r_hits``).  A pair is vetoed (returns False)
    when at least one overlapping layer exists AND the minimum 3D separation there
    is smaller than ``cut_mm``.

    Conservative by design: if either segment has no high-r hits, or there are no
    overlapping layers, the pair passes (True) — no information means no veto.

    Args:
        hi, hj:  (rs, xs, ys, zs) arrays from ``precompute_high_r_hits``.
        r_tol:   Max |r_A - r_B| to consider two hits co-layer (mm).
        cut_mm:  Minimum required 3D separation at overlapping layers (mm).
    """
    ri, xi, yi, zi = hi
    rj, xj, yj, zj = hj
    if len(ri) == 0 or len(rj) == 0:
        return True   # no high-r hits on one side — can't veto

    dr      = np.abs(ri[:, None] - rj[None, :])   # (Hi, Hj)
    nearest = np.argmin(dr, axis=1)                # (Hi,)
    mask    = dr[np.arange(len(ri)), nearest] < r_tol
    if not mask.any():
        return True   # no overlapping layers — can't veto

    n_idx = nearest[mask]
    dx = xi[mask] - xj[n_idx]
    dy = yi[mask] - yj[n_idx]
    dz = zi[mask] - zj[n_idx]
    return float(np.min(np.sqrt(dx**2 + dy**2 + dz**2))) >= cut_mm


def geometric_cut_batch(high_r_hits, cand_i, cand_j, r_tol, cut_mm):
    """Vectorised passes_geometric_cut over all candidate pairs at once (same result, no python loop)."""
    ci = np.asarray(cand_i); cj = np.asarray(cand_j)
    P = ci.shape[0]
    if P == 0:
        return np.ones(0, dtype=bool)
    lens = np.array([len(h[0]) for h in high_r_hits], dtype=np.int64)
    H = int(lens.max()) if lens.size else 0
    if H == 0:
        return np.ones(P, dtype=bool)
    Nseg = len(high_r_hits)
    R = np.full((Nseg, H), np.inf); X = np.zeros((Nseg, H)); Y = np.zeros((Nseg, H)); Z = np.zeros((Nseg, H))
    M = np.zeros((Nseg, H), dtype=bool)
    for k, (rs, xs, ys, zs) in enumerate(high_r_hits):
        L = len(rs)
        if L:
            R[k, :L] = rs; X[k, :L] = xs; Y[k, :L] = ys; Z[k, :L] = zs; M[k, :L] = True

    ri, Mi = R[ci], M[ci]
    rj, Mj = R[cj], M[cj]
    with np.errstate(invalid="ignore"):
        dr = np.abs(ri[:, :, None] - rj[:, None, :])
    dr = np.where(Mj[:, None, :], dr, np.inf)               # nearest valid j-hit per i-hit
    nearest = np.argmin(dr, axis=2)
    nearest_dr = np.take_along_axis(dr, nearest[:, :, None], axis=2)[:, :, 0]
    co_layer = Mi & (nearest_dr < r_tol)
    xj_n = np.take_along_axis(X[cj], nearest, axis=1)
    yj_n = np.take_along_axis(Y[cj], nearest, axis=1)
    zj_n = np.take_along_axis(Z[cj], nearest, axis=1)
    d3 = np.sqrt((X[ci] - xj_n) ** 2 + (Y[ci] - yj_n) ** 2 + (Z[ci] - zj_n) ** 2)
    min_d3 = np.where(co_layer, d3, np.inf).min(axis=1)
    empty = (lens[ci] == 0) | (lens[cj] == 0)
    return empty | ~np.isfinite(min_d3) | (min_d3 >= cut_mm)


# ─── Track Assembly ────────────────────────────────────────────────────────


def segments_to_track_labels(
    matched_tracks: List[List[SegmentInfo]],
    unmatched: List[SegmentInfo],
    num_nodes: int,
    hit_t: Optional[np.ndarray] = None,
) -> torch.Tensor:
    """
    Convert matched tracks and unmatched segments to hit_track_labels tensor.

    Hits within each track are combined; within matched tracks, ordering is
    determined by hit_t if available.

    Args:
        matched_tracks: List of matched track groups.
        unmatched: List of unmatched segments (each becomes a standalone track).
        num_nodes: Total number of nodes in the graph.
        hit_t: Optional time array for ordering hits within tracks.

    Returns:
        hit_track_labels: Tensor of shape (num_nodes,), -1 for unassigned.
    """
    labels = torch.full((num_nodes,), -1, dtype=torch.long)
    track_id = 0
    idx_all, lab_all = [], []   # accumulate (hit, track_id) for one vectorised scatter at the end

    # Matched tracks: each group of segments = one track. (hit_t is unused — ordering hits
    # within a track does not change their shared label, so the old per-hit argsort was a no-op.)
    for track_segments in matched_tracks:
        hits = set()
        for seg in track_segments:
            hits.update(seg.hits)   # dedup across segments of the same track
        idx_all.extend(hits); lab_all += [track_id] * len(hits)
        track_id += 1

    # Unmatched segments: each becomes a standalone track
    for seg in unmatched:
        if len(seg.hits) == 0:
            continue
        idx_all.extend(seg.hits); lab_all += [track_id] * len(seg.hits)
        track_id += 1

    if idx_all:   # one assignment instead of ~N_hits torch scalar sets
        labels[torch.as_tensor(idx_all, dtype=torch.long)] = torch.as_tensor(lab_all, dtype=torch.long)

    return labels


# ─── High-Level Entry Point ────────────────────────────────────────────────


def build_tracks_for_event(
    graph,
    config: dict,
    return_matching_info: bool = False,
) -> Tuple[torch.Tensor, dict]:
    """
    Build tracks for a single event using helix-based segment matching.

    Steps:
        1. Extract segments (CC or ground truth)
        2. Fit helix to each segment
        3. Match segments by helix parameter compatibility
        4. Convert to hit_track_labels (short/unfitted segments become standalone)

    Args:
        graph: PyG Data object from data/gnn_stage/.
        config: Configuration dict with all parameters.
        return_matching_info: If True, return a third element — a dict with
            'segments', 'matched_tracks', and 'unmatched' — for downstream
            matching-efficiency evaluation.

    Returns:
        hit_track_labels: Tensor of shape (num_nodes,).
        stats: Dict with statistics about the matching.
        matching_info (only if return_matching_info=True): dict with keys
            'segments', 'matched_tracks', 'unmatched'.
    """
    use_gt = config.get("use_gt_segments", False)
    score_cut = config.get("score_cut", 0.5)
    use_wrangler = config.get("use_wrangler", False)
    B_field = config.get("B_field", 2.0)
    outlier_rejection = config.get("outlier_rejection", False)
    matching_config = config.get("matching", {})

    from time import perf_counter

    # Step 1: Extract segments  (CC / Wrangler on edge_scores — shared with the GNN matcher)
    t0 = perf_counter()
    n_wrangler_splits = 0
    if use_gt:
        segments = extract_segments_from_ground_truth(graph)
    elif use_wrangler:
        from low_pt_custom_utils.wrangler_utils import extract_segments_with_wrangler
        segments, n_wrangler_splits = extract_segments_with_wrangler(graph, score_cut)
    else:
        segments = extract_segments_from_cc(graph, score_cut)
    t_extract = perf_counter() - t0

    n_segments = len(segments)

    # Step 2: Fit helices (Kasa circle fit + pitch, per segment)
    t0 = perf_counter()
    segments = fit_helices_to_segments(segments, graph, B_field=B_field, outlier_rejection=outlier_rejection)
    t_fit = perf_counter() - t0

    n_good_fits = sum(1 for s in segments if s.helix and s.helix.fit_quality == "good")
    n_poor_fits = sum(1 for s in segments if s.helix and s.helix.fit_quality == "poor")
    n_no_fits = sum(1 for s in segments if s.helix and s.helix.fit_quality == "none")

    # Step 3: Match segments by helix parameters
    t0 = perf_counter()
    mtiming = {}
    matched_tracks, unmatched = match_segments(segments, matching_config, timing=mtiming)
    t_match = perf_counter() - t0

    # Step 4: Convert to labels (short segments become standalone tracks)
    num_nodes = graph.hit_x.size(0)
    hit_t = graph.hit_t.cpu().numpy() if hasattr(graph, 'hit_t') else None
    labels = segments_to_track_labels(matched_tracks, unmatched, num_nodes, hit_t=hit_t)

    # Statistics
    outer_r_threshold = matching_config.get("outer_r_threshold", 1000.0)
    n_matched_pairs = len(matched_tracks)
    n_complete_tracks = sum(1 for s in unmatched if s.outer_r >= outer_r_threshold)
    n_no_match = len(unmatched) - n_complete_tracks
    n_standalone = len(unmatched)
    n_total_tracks = n_matched_pairs + n_standalone
    n_assigned = (labels >= 0).sum().item()
    n_unassigned = num_nodes - n_assigned

    stats = {
        "n_segments": n_segments,
        "n_wrangler_splits": n_wrangler_splits,
        "n_good_fits": n_good_fits,
        "n_poor_fits": n_poor_fits,
        "n_no_fits": n_no_fits,
        "n_matched_pairs": n_matched_pairs,
        "n_complete_tracks": n_complete_tracks,
        "n_no_match": n_no_match,
        "n_standalone": n_standalone,
        "n_total_tracks": n_total_tracks,
        "n_assigned_hits": n_assigned,
        "n_unassigned_hits": n_unassigned,
        "time_extract_s": t_extract,
        "time_fit_s": t_fit,
        "time_match_s": t_match,
        "time_score_s": mtiming.get("score", 0.0),    # vectorised pairwise scoring
        "time_greedy_s": mtiming.get("greedy", 0.0),  # greedy assignment loop
    }

    if return_matching_info:
        matching_info = {
            "segments": segments,
            "matched_tracks": matched_tracks,
            "unmatched": unmatched,
        }
        return labels, stats, matching_info

    return labels, stats
