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
    """
    x = graph.hit_x.cpu().numpy()
    y = graph.hit_y.cpu().numpy()
    z = graph.hit_z.cpu().numpy()

    for seg in segments:
        seg_x = x[seg.hits]
        seg_y = y[seg.hits]
        seg_z = z[seg.hits]

        seg.helix = fit_helix_to_segment(
            seg_x, seg_y, seg_z,
            B_field=B_field,
            outlier_rejection=outlier_rejection,
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

    # Step 3: Compute all pairwise scores between fittable, non-complete segments
    candidates = []
    for ai in range(len(fittable)):
        for bi in range(ai + 1, len(fittable)):
            idx_a, seg_a = fittable[ai]
            idx_b, seg_b = fittable[bi]
            score, compatible = compute_matching_score(seg_a, seg_b, config)
            if compatible and score > 0:
                candidates.append((score, idx_a, idx_b))

    # Sort by score descending
    candidates.sort(key=lambda x: x[0], reverse=True)

    # Greedy matching
    matched = set()
    matched_tracks = []

    for score, idx_a, idx_b in candidates:
        if idx_a in matched or idx_b in matched:
            continue
        # Match these two segments into one track
        matched_tracks.append([segments[idx_a], segments[idx_b]])
        matched.add(idx_a)
        matched.add(idx_b)

    # Collect unmatched segments (including short ones)
    unmatched = []
    for i, seg in enumerate(segments):
        if i not in matched:
            unmatched.append(seg)

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
    labels = torch.ones(num_nodes, dtype=torch.long) * -1
    track_id = 0

    # Matched tracks: each group of segments = one track
    for track_segments in matched_tracks:
        # Collect all hits from all segments in this track
        all_hits = []
        for seg in track_segments:
            all_hits.extend(seg.hits)

        # Remove duplicates (shouldn't happen but be safe)
        all_hits = list(set(all_hits))

        # Optionally order by time
        if hit_t is not None and len(all_hits) > 0:
            hit_times = np.asarray(hit_t[all_hits], dtype=np.float64)
            order = np.argsort(hit_times)
            all_hits = [all_hits[i] for i in order]

        for hit_idx in all_hits:
            labels[hit_idx] = track_id
        track_id += 1

    # Unmatched segments: each becomes a standalone track
    for seg in unmatched:
        if len(seg.hits) == 0:
            continue
        for hit_idx in seg.hits:
            labels[hit_idx] = track_id
        track_id += 1

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

    # Step 1: Extract segments
    n_wrangler_splits = 0
    if use_gt:
        segments = extract_segments_from_ground_truth(graph)
    elif use_wrangler:
        from low_pt_custom_utils.wrangler_utils import extract_segments_with_wrangler
        segments, n_wrangler_splits = extract_segments_with_wrangler(graph, score_cut)
    else:
        segments = extract_segments_from_cc(graph, score_cut)

    n_segments = len(segments)

    # Step 2: Fit helices
    segments = fit_helices_to_segments(segments, graph, B_field=B_field, outlier_rejection=outlier_rejection)

    n_good_fits = sum(1 for s in segments if s.helix and s.helix.fit_quality == "good")
    n_poor_fits = sum(1 for s in segments if s.helix and s.helix.fit_quality == "poor")
    n_no_fits = sum(1 for s in segments if s.helix and s.helix.fit_quality == "none")

    # Step 3: Match segments
    matched_tracks, unmatched = match_segments(segments, matching_config)

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
    }

    if return_matching_info:
        matching_info = {
            "segments": segments,
            "matched_tracks": matched_tracks,
            "unmatched": unmatched,
        }
        return labels, stats, matching_info

    return labels, stats
