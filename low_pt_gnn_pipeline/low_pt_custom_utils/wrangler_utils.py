"""
Wrangler (walk-through) disentanglement utilities.

Implements the Wrangler algorithm for resolving ambiguous CC subgraphs:
  1. Build a directed graph oriented radially outward (R = r² + z²).
  2. Find source nodes (in-degree == 0).
  3. From each source, greedily walk following the highest-scoring edge.
  4. Keep the longest path per source; assign hits greedily across paths.

Reference: matches the ACORN CCandWalk / cc_and_walk_utils convention for
edge orientation (3D distance from origin, not transverse r alone).
"""

from collections import deque
from typing import List, Tuple

import numpy as np
import scipy.sparse as sps
import torch
from torch_geometric.utils import remove_isolated_nodes, to_scipy_sparse_matrix


from low_pt_custom_utils.segment_matching import SegmentInfo, _build_segment_info


# ─── Directed Graph Construction ─────────────────────────────────────────────


def _build_directed_adj(hit_set, src_arr, dst_arr, scores_arr, R3d):
    """
    Build a directed adjacency dict for a CC cluster, oriented radially outward.

    For each edge (u, v) within the cluster: direction is from smaller 3D distance
    to larger (R = r² + z²), matching the ACORN remove_cycles convention.
    Ties broken by smaller node index first. Deduplicates symmetric edges
    (common in GNN edge_index), keeping max score.

    Returns:
        adj: dict {u: {v: score}} for all nodes in hit_set.
    """
    adj = {n: {} for n in hit_set}
    for s, d, sc in zip(src_arr, dst_arr, scores_arr):
        s, d = int(s), int(d)
        if s not in hit_set or d not in hit_set or s == d:
            continue
        # Orient outward by 3D distance from origin
        if R3d[s] < R3d[d] or (R3d[s] == R3d[d] and s < d):
            u, v = s, d
        else:
            u, v = d, s
        # Keep max score for this directed edge
        if v not in adj[u] or sc > adj[u][v]:
            adj[u][v] = float(sc)
    return adj


# ─── Walk ─────────────────────────────────────────────────────────────────────


def _wrangler_walk(root, adj):
    """
    Wrangler traversal from a single source node.

    At each step follows the highest-scoring outgoing edge (greedy).
    Stops when no unvisited outgoing edges remain.

    Returns:
        The single greedy path (list of node indices).
    """
    path = [root]
    visited = {root}

    while True:
        current = path[-1]
        outgoing = [
            (v, sc) for v, sc in adj.get(current, {}).items()
            if v not in visited
        ]
        if not outgoing:
            break
        best_v = max(outgoing, key=lambda x: x[1])[0]
        path.append(best_v)
        visited.add(best_v)

    return path


# ─── Per-Cluster Disentanglement ──────────────────────────────────────────────


def _wrangler_disentangle_cluster(hit_indices, src_arr, dst_arr, scores_arr,
                                   r_coords, x, y, z):
    """
    Run Wrangler on a single CC cluster.

    Steps:
      1. Build directed graph (radially outward) within the cluster.
      2. Find source nodes (in-degree == 0); fall back to innermost-R3d node.
      3. Run Wrangler walk (greedy) from each source.
      4. Assign hits greedily to paths (longest wins on conflict).
      5. Any unassigned hits become a leftover segment.

    Returns:
        List of SegmentInfo objects (one per disentangled path / leftover group).
    """
    hit_set = set(hit_indices.tolist())

    # 3D distance from origin for edge orientation (matches ACORN remove_cycles)
    R3d = r_coords**2 + z**2

    adj = _build_directed_adj(hit_set, src_arr, dst_arr, scores_arr, R3d)

    # Compute in-degree
    in_degree = {n: 0 for n in hit_set}
    for neighbors in adj.values():
        for v in neighbors:
            in_degree[v] += 1

    # Source nodes: in-degree == 0
    sources = [n for n in hit_set if in_degree[n] == 0]
    if not sources:
        # Fallback: hit with smallest 3D distance as single source
        sources = [int(hit_indices[np.argmin(R3d[hit_indices])])]

    # Walk from each source
    all_paths = []
    for src in sources:
        path = _wrangler_walk(src, adj)
        if path:
            all_paths.append(path)

    # Assign hits to paths greedily (longest path wins on conflict)
    all_paths.sort(key=len, reverse=True)
    used = set()
    segments = []
    for path in all_paths:
        clean = [n for n in path if n not in used]
        if clean:
            used.update(clean)
            seg = _build_segment_info(clean, x, y, z, r_coords, cluster_id=-1)
            segments.append(seg)

    # Leftover hits (not reached by any walk) → one fallback segment
    leftover = [n for n in hit_set if n not in used]
    if leftover:
        seg = _build_segment_info(leftover, x, y, z, r_coords, cluster_id=-1)
        segments.append(seg)

    return segments


# ─── Top-Level Entry Point ────────────────────────────────────────────────────


def extract_segments_with_wrangler(graph, score_cut):
    """
    CC clustering followed by Wrangler disentanglement.

    Each CC cluster is disentangled into one or more path segments using
    the Wrangler walk-through algorithm (greedy, no branching). All edges
    already passed score_cut, so no additional score filter is applied.

    Args:
        graph:     PyG Data object with edge_index, edge_scores, hit_r, etc.
        score_cut: Edge score threshold for CC clustering.

    Returns:
        segments:          List of SegmentInfo objects.
        n_wrangler_splits: Number of CC clusters split into >1 segment.
    """
    edge_mask = graph.edge_scores > score_cut
    edges_filtered = graph.edge_index[:, edge_mask]
    scores_filtered = graph.edge_scores[edge_mask].cpu().numpy()
    src_arr = edges_filtered[0].cpu().numpy()
    dst_arr = edges_filtered[1].cpu().numpy()

    num_nodes = graph.hit_x.size(0)
    r = graph.hit_r.cpu().numpy()
    x = graph.hit_x.cpu().numpy()
    y = graph.hit_y.cpu().numpy()
    z = graph.hit_z.cpu().numpy()

    # CC clustering
    edges_clean, _, node_mask = remove_isolated_nodes(edges_filtered, num_nodes=num_nodes)
    n_connected = node_mask.sum().item()
    if n_connected == 0:
        return [], 0

    sparse_mat = to_scipy_sparse_matrix(edges_clean, num_nodes=n_connected)
    n_components, labels_compact = sps.csgraph.connected_components(
        sparse_mat, directed=False, return_labels=True
    )

    labels = torch.ones(num_nodes, dtype=torch.long) * -1
    labels[node_mask] = torch.tensor(labels_compact, dtype=torch.long)
    labels = labels.numpy()

    segments = []
    n_wrangler_splits = 0

    for cluster_id in range(n_components):
        hit_indices = np.where(labels == cluster_id)[0]
        if len(hit_indices) < 1:
            continue

        cluster_segs = _wrangler_disentangle_cluster(
            hit_indices, src_arr, dst_arr, scores_filtered,
            r, x, y, z,
        )

        if len(cluster_segs) > 1:
            n_wrangler_splits += 1

        segments.extend(cluster_segs)

    return segments, n_wrangler_splits
