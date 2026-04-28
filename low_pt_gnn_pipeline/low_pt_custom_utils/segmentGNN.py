"""
Shared GNN-based segment matching utilities.

Used by:
  - track_building_stage_(3)/GNN_segment_matching_track_builder.py
  - track_building_stage_(3)/GNN_segmentmatcher_walkthrough.py
"""

import torch
from torch_geometric.data import Batch

from low_pt_custom_utils.segment_matching import (
    extract_segments_from_cc,
    extract_segments_from_ground_truth,
    segments_to_track_labels,
    precompute_high_r_hits,
    passes_geometric_cut,
)
from low_pt_custom_utils.mini_gnn_segment_embedding import segment_to_pyg
from low_pt_custom_utils.wrangler_utils import extract_segments_with_wrangler


def match_segments_gnn(segments, graph, model, config, device="cpu"):
    """
    Match segments using cosine similarity of GNN embeddings.

    Args:
        segments:  List[SegmentInfo] from Wrangler/CC/GT extraction.
        graph:     PyG Data object for the full event (provides hit coords).
        model:     Pretrained SegmentGNN in eval mode (has .node_scales attribute).
        config:    Dict with matching thresholds and optional geometric cut params.
        device:    Torch device string.

    Returns:
        matched_tracks: List of [SegmentInfo, SegmentInfo] pairs.
        unmatched:      List of SegmentInfo not in any matched pair
                        (includes complete tracks).
    """
    cos_sim_threshold = config["cos_sim_threshold"]
    outer_r_threshold = config["outer_r_threshold"]
    node_scales       = model.node_scales

    geo_cut_enabled = config["geo_cut_enabled"]
    geo_r_high      = config["geo_cut_r_high_mm"]
    geo_r_tol       = config["geo_cut_r_tol_mm"]
    geo_min_3d      = config["geo_cut_min_3d_mm"]

    # Separate complete tracks (reach outer detector)
    to_match, complete = [], []
    for seg in segments:
        (complete if seg.outer_r >= outer_r_threshold else to_match).append(seg)

    if not to_match:
        return [], complete

    # Embed all segments in one batched forward pass
    seg_graphs = [segment_to_pyg(seg, graph, node_scales=tuple(node_scales))
                  for seg in to_match]
    batch_data = Batch.from_data_list(seg_graphs).to(device)
    with torch.no_grad():
        embeddings = model(batch_data.x, batch_data.edge_index, batch_data.batch)

    sim_matrix = (embeddings @ embeddings.T).cpu()  # (N, N), cosine similarity

    if geo_cut_enabled:
        high_r_hits = precompute_high_r_hits(to_match, graph, geo_r_high)

    n = len(to_match)
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            s = sim_matrix[i, j].item()
            if s < cos_sim_threshold:
                continue
            if geo_cut_enabled and not passes_geometric_cut(
                    high_r_hits[i], high_r_hits[j], geo_r_tol, geo_min_3d):
                continue
            pairs.append((s, i, j))

    pairs.sort(key=lambda x: x[0], reverse=True)

    matched_tracks, matched_set = [], set()
    for score, i, j in pairs:
        if i not in matched_set and j not in matched_set:
            matched_tracks.append([to_match[i], to_match[j]])
            matched_set |= {i, j}

    unmatched = complete + [to_match[i] for i in range(n) if i not in matched_set]
    return matched_tracks, unmatched


def build_tracks_for_event_gnn(graph, config, model, device, return_matching_info=False):
    """
    Build tracks for one event using GNN-scored segment matching.

    Segment extraction order:
      - use_gt_segments=True  → ground truth segments (skip CC + Wrangler)
      - use_wrangler=True     → CC + Wrangler disentanglement
      - use_wrangler=False    → plain CC (default)

    Returns:
        hit_track_labels: Tensor of shape (num_nodes,), -1 for unassigned hits.
        stats: Dict with per-event counts (always includes n_wrangler_splits,
               which is 0 when wrangler is not used).
        matching_info (only if return_matching_info=True): Dict with keys
            'segments', 'matched_tracks', 'unmatched'.
    """
    use_gt       = config["use_gt_segments"]
    score_cut    = config["score_cut"]
    use_wrangler = config.get("use_wrangler", False)
    gnn_config   = config["gnn"]

    n_wrangler_splits = 0

    # Step 1: Extract segments
    if use_gt:
        segments = extract_segments_from_ground_truth(graph)
    elif use_wrangler:
        segments, n_wrangler_splits = extract_segments_with_wrangler(graph, score_cut)
    else:
        segments = extract_segments_from_cc(graph, score_cut)

    n_segments = len(segments)

    # Step 2+3: GNN-based matching (includes complete track separation)
    matched_tracks, unmatched = match_segments_gnn(
        segments, graph, model, gnn_config, device=device,
    )

    # Step 4: Convert to hit labels
    num_nodes = graph.hit_x.size(0)
    hit_t = graph.hit_t.cpu().numpy() if hasattr(graph, "hit_t") else None
    labels = segments_to_track_labels(matched_tracks, unmatched, num_nodes, hit_t=hit_t)

    # Statistics
    outer_r_threshold = gnn_config["outer_r_threshold"]
    n_matched_pairs   = len(matched_tracks)
    n_complete_tracks = sum(1 for s in unmatched if s.outer_r >= outer_r_threshold)
    n_no_match        = len(unmatched) - n_complete_tracks
    n_standalone      = len(unmatched)
    n_total_tracks    = n_matched_pairs + n_standalone
    n_assigned        = (labels >= 0).sum().item()
    n_unassigned      = num_nodes - n_assigned

    stats = {
        "n_segments": n_segments,
        "n_wrangler_splits": n_wrangler_splits,
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
