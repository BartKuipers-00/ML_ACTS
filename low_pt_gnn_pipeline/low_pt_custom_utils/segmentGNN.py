"""
Shared GNN-based segment matching utilities.

Used by:
  - track_building_stage_(3)/GNN_segment_matching_track_builder.py
  - track_building_stage_(3)/GNN_segmentmatcher_walkthrough.py
"""

from time import perf_counter

import numpy as np
import torch

from low_pt_custom_utils.segment_matching import (
    extract_segments_from_cc,
    extract_segments_from_ground_truth,
    segments_to_track_labels,
    precompute_high_r_hits,
    geometric_cut_batch,
)
from low_pt_custom_utils.mini_gnn_segment_embedding import build_segment_batch
from low_pt_custom_utils.wrangler_utils import extract_segments_with_wrangler


def match_segments_gnn(segments, graph, model, config, device="cpu", timing=None):
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

    _cuda = str(device).startswith("cuda")
    def _sync():
        if _cuda:
            torch.cuda.synchronize()

    # Embed all segments in one batched forward pass. Finer split: CPU-side batch assembly,
    # host->device transfer, and the actual GNN forward (the only GPU-accelerable part).
    t0 = perf_counter()
    x, edge_index, batch = build_segment_batch(
        to_match, graph, node_scales=tuple(node_scales))
    t_assemble = perf_counter() - t0

    t0 = perf_counter()
    x = x.to(device)
    edge_index = edge_index.to(device)
    batch = batch.to(device)
    _sync()
    t_to_device = perf_counter() - t0

    t0 = perf_counter()
    with torch.no_grad():
        embeddings = model(x, edge_index, batch)
    _sync()
    t_forward = perf_counter() - t0
    t_embed = t_assemble + t_to_device + t_forward

    # Dot-product cosine matching (GPU matmul + copy back) then threshold/greedy (CPU).
    t0 = perf_counter()
    sim_matrix = embeddings @ embeddings.T
    _sync()
    sim_matrix = sim_matrix.cpu()  # (N, N), cosine similarity
    t_simil = perf_counter() - t0
    t_match0 = perf_counter()

    if geo_cut_enabled:
        high_r_hits = precompute_high_r_hits(to_match, graph, geo_r_high)

    n = len(to_match)
    # Threshold the upper triangle in one op, then apply the geometric veto VECTORISED over
    # all survivors (instead of a python loop of passes_geometric_cut per pair).
    triu_i, triu_j = torch.triu_indices(n, n, offset=1)
    cand_scores = sim_matrix[triu_i, triu_j]
    keep = cand_scores >= cos_sim_threshold
    cand_i = triu_i[keep].cpu().numpy()
    cand_j = triu_j[keep].cpu().numpy()
    cand_s = cand_scores[keep].cpu().numpy()

    if geo_cut_enabled and cand_i.shape[0]:
        gkeep = geometric_cut_batch(high_r_hits, cand_i, cand_j, geo_r_tol, geo_min_3d)
        cand_i, cand_j, cand_s = cand_i[gkeep], cand_j[gkeep], cand_s[gkeep]

    pairs = [(float(s), int(i), int(j)) for s, i, j in zip(cand_s, cand_i, cand_j)]

    pairs.sort(key=lambda x: x[0], reverse=True)

    matched_tracks, matched_set = [], set()
    for score, i, j in pairs:
        if i not in matched_set and j not in matched_set:
            matched_tracks.append([to_match[i], to_match[j]])
            matched_set |= {i, j}

    unmatched = complete + [to_match[i] for i in range(n) if i not in matched_set]
    if timing is not None:
        t_greedy = perf_counter() - t_match0
        timing.update(assemble=t_assemble, to_device=t_to_device, forward=t_forward,
                      simil=t_simil, greedy=t_greedy,
                      embed=t_embed, dotmatch=t_simil + t_greedy)
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

    # Step 1: Extract segments  ("track building": CC score cut + Wrangler disentangle)
    t0 = perf_counter()
    if use_gt:
        segments = extract_segments_from_ground_truth(graph)
    elif use_wrangler:
        segments, n_wrangler_splits = extract_segments_with_wrangler(graph, score_cut)
    else:
        segments = extract_segments_from_cc(graph, score_cut)
    t_extract = perf_counter() - t0

    n_segments = len(segments)

    # Step 2+3: GNN-based matching  ("segment matching": miniGNN embed + cosine match)
    mtiming = {}
    t0 = perf_counter()
    matched_tracks, unmatched = match_segments_gnn(
        segments, graph, model, gnn_config, device=device, timing=mtiming,
    )
    t_match = perf_counter() - t0

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
        # Per-event sub-stage timings (seconds): segment extraction (CC + Wrangler)
        # vs miniGNN segment matching. Surfaced so the benchmark can split stage 4.
        "time_extract_s": t_extract,
        "time_match_s": t_match,
        "time_embed_s": mtiming.get("embed", 0.0),       # miniGNN forward pass (embedding)
        "time_dotmatch_s": mtiming.get("dotmatch", 0.0),  # dot-product cosine match + greedy
        "time_assemble_s": mtiming.get("assemble", 0.0),  # CPU-side batch assembly
        "time_todevice_s": mtiming.get("to_device", 0.0),  # host->device transfer
        "time_forward_s": mtiming.get("forward", 0.0),     # GNN forward only (GPU-accelerable)
        "time_simil_s": mtiming.get("simil", 0.0),         # NxN matmul + copy back
        "time_greedy_s": mtiming.get("greedy", 0.0),       # threshold + greedy assignment
    }

    if return_matching_info:
        matching_info = {
            "segments": segments,
            "matched_tracks": matched_tracks,
            "unmatched": unmatched,
        }
        return labels, stats, matching_info

    return labels, stats
