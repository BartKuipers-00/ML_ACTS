"""
Mini-GNN for per-segment latent embeddings.

Operates on individual detector hit segments: each segment is a small graph
(hits as nodes, fully-connected directed edges) and the GNN produces a
graph-level embedding via global mean pooling.

Training: Supervised Contrastive (SupCon) loss on segment embeddings.
  - Positive pairs  (same majority hit_particle_id): pull together on unit sphere
  - Negative pairs  (different particles):          push apart via softmax denominator

Architecture:
    BatchNorm(4) → [SAGEConv → BatchNorm → ReLU → Dropout] × n_layers
    → global_mean_pool → MLP head(hidden → emb_dim) → L2 normalize
"""

from typing import List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import SAGEConv, global_mean_pool


# Number of node features per hit: [x/1000, y/1000, z/500, r/1000]  (Cartesian + radial)
NODE_FEATURE_DIM = 4


# ─── Model ──────────────────────────────────────────────────────────────────


class SegmentGNN(nn.Module):
    """
    Mini-GNN that maps a hit segment to an L2-normalized embedding vector.

    Input:  Small PyG graph (n_hits nodes, fully-connected edges).
    Output: (emb_dim,) L2-normalized vector per segment (graph-level readout).

    The L2 normalization means: cosine_similarity(a, b) = dot(a, b), which
    makes the embeddings directly comparable via dot product.

    Args:
        node_in_dim: Number of input node features (default 4).
        hidden_dim:  Width of all GNN hidden layers.
        emb_dim:     Dimension of the output embedding space.
        n_layers:    Number of SAGEConv message-passing layers.
        dropout:     Dropout probability after each hidden activation.
        proj_dim:    Width of the projection head hidden layers (default: hidden_dim).
        proj_layers: Number of hidden layers in the projection head (default: 1).
    """

    def __init__(
        self,
        node_in_dim: int,
        hidden_dim: int,
        emb_dim: int,
        n_layers: int,
        dropout: float,
        proj_dim: int,
        proj_layers: int,
    ):
        super().__init__()

        self.input_norm = nn.BatchNorm1d(node_in_dim)

        self.convs = nn.ModuleList()
        self.norms = nn.ModuleList()
        in_dim = node_in_dim
        for _ in range(n_layers):
            self.convs.append(SAGEConv(in_dim, hidden_dim))
            self.norms.append(nn.BatchNorm1d(hidden_dim))
            in_dim = hidden_dim

        self.dropout = nn.Dropout(dropout)
        self.activation = nn.ReLU()

        # Projection head: hidden_dim → [proj_dim → ReLU] × proj_layers → emb_dim
        proj_dim = proj_dim if proj_dim is not None else hidden_dim
        head_layers = []
        head_in = hidden_dim
        for _ in range(proj_layers):
            head_layers.append(nn.Linear(head_in, proj_dim))
            head_layers.append(nn.ReLU())
            head_in = proj_dim
        head_layers.append(nn.Linear(head_in, emb_dim))
        self.head = nn.Sequential(*head_layers)

    def forward(
        self,
        x: torch.Tensor,
        edge_index: torch.Tensor,
        batch: torch.Tensor,
    ) -> torch.Tensor:
        """
        Args:
            x:          (total_nodes, node_in_dim) node feature matrix.
            edge_index: (2, total_edges) edge indices (batch-remapped by Batch).
            batch:      (total_nodes,) graph assignment vector.

        Returns:
            (n_graphs, emb_dim) L2-normalized embedding tensor.
        """
        x = self.input_norm(x)

        for conv, norm in zip(self.convs, self.norms):
            x = conv(x, edge_index)
            x = norm(x)
            x = self.activation(x)
            x = self.dropout(x)

        # Graph-level readout: mean over all nodes in each graph
        x = global_mean_pool(x, batch)   # (n_graphs, hidden_dim)

        # Project to embedding space
        x = self.head(x)                  # (n_graphs, emb_dim)

        # L2 normalize: cos_sim(a, b) = a · b for unit vectors
        return F.normalize(x, p=2, dim=-1)


# ─── Graph Construction ──────────────────────────────────────────────────────


def segment_to_pyg(seg, graph, node_scales=(1000.0, 1000.0, 500.0, 1000.0)) -> Data:
    """
    Build a PyG Data object from a SegmentInfo and its parent event graph.

    Node features (4-dim, per hit):
        [hit_x / node_scales[0], hit_y / node_scales[1],
         hit_z / node_scales[2], hit_r / node_scales[3]]


    Args:
        seg:         SegmentInfo object (from segment_matching.py).
        graph:       PyG Data object for the full event (provides hit coordinates).
        node_scales: Tuple/list of four normalization scales for (x, y, z, r).
                     Default (1000, 1000, 500, 1000).

    Returns:
        PyG Data object with x and edge_index (local node indexing 0..n-1).
    """
    hits = list(seg.hits)
    n = len(hits)

    x_f = graph.hit_x[hits].float() / node_scales[0]
    y_f = graph.hit_y[hits].float() / node_scales[1]
    z_f = graph.hit_z[hits].float() / node_scales[2]
    r_f = graph.hit_r[hits].float() / node_scales[3]

    node_feats = torch.stack([x_f, y_f, z_f, r_f], dim=1)  # (n, 4)

    # Fully-connected directed edges (all i→j with i ≠ j)
    if n > 1:
        idx = torch.arange(n, dtype=torch.long)
        src = idx.repeat_interleave(n)  # [0,0,...,1,1,...,n-1,n-1,...]
        dst = idx.repeat(n)             # [0,1,...,n-1,0,1,...,n-1,...]
        mask = src != dst
        edge_index = torch.stack([src[mask], dst[mask]])  # (2, n*(n-1))
    else:
        edge_index = torch.zeros((2, 0), dtype=torch.long)

    return Data(x=node_feats, edge_index=edge_index)


def build_segment_batch(segments, graph, node_scales=(1000.0, 1000.0, 500.0, 1000.0)):
    """
    Vectorized batched miniGNN input for many segments at once.

    Equivalent to Batch.from_data_list([segment_to_pyg(s, graph) ...]) but built in
    a few tensor ops instead of one PyG Data per segment. Returns (x, edge_index,
    batch) with block-diagonal fully-connected (i != j) edges per segment.
    """
    sizes = torch.tensor([len(s.hits) for s in segments], dtype=torch.long)
    flat = []
    for s in segments:
        flat.extend(s.hits)
    all_hits = torch.as_tensor(flat, dtype=torch.long)

    x = torch.stack([
        graph.hit_x[all_hits].float() / node_scales[0],
        graph.hit_y[all_hits].float() / node_scales[1],
        graph.hit_z[all_hits].float() / node_scales[2],
        graph.hit_r[all_hits].float() / node_scales[3],
    ], dim=1)

    N = sizes.numel()
    batch = torch.repeat_interleave(torch.arange(N), sizes)

    # Per (i,j) slot: recover its block, within-block position p, then a=p//n, b=p%n
    offsets = torch.cumsum(sizes, 0) - sizes
    n_sq = sizes * sizes
    total = int(n_sq.sum())
    if total == 0:
        return x, torch.zeros((2, 0), dtype=torch.long), batch
    elem_block = torch.repeat_interleave(torch.arange(N), n_sq)
    sq_starts = torch.cumsum(n_sq, 0) - n_sq
    p = torch.arange(total) - sq_starts[elem_block]
    n_per = sizes[elem_block]
    off_per = offsets[elem_block]
    a = p // n_per
    b = p % n_per
    mask = a != b
    edge_index = torch.stack([off_per[mask] + a[mask], off_per[mask] + b[mask]])
    return x, edge_index, batch


def get_segment_particle_id(seg, graph) -> int:
    """
    Return the majority hit_particle_id for a segment, excluding noise (pid=0).

    Returns 0 if the segment has no signal hits (all noise or empty).
    """
    particle_ids = graph.hit_particle_id.cpu().numpy()
    pids = particle_ids[np.array(list(seg.hits), dtype=int)]
    signal_pids = pids[pids > 0]
    if len(signal_pids) == 0:
        return 0
    values, counts = np.unique(signal_pids, return_counts=True)
    return int(values[np.argmax(counts)])


def filter_paired_segments(segments, graph, min_hits=2):
    """
    Pairs-only + min-hits filter, applied jointly.

    Keep a segment iff it has >= min_hits hits AND its particle has >= 2 such
    (>= min_hits-hit) segments. The hit cut is applied first so that dropping a
    short arc never orphans its multi-hit partner; particles left with fewer
    than 2 qualifying segments (and all noise, pid=0) are dropped entirely.

    Returns parallel lists (kept_segments, kept_particle_ids).
    """
    candidates = [
        (seg, get_segment_particle_id(seg, graph))
        for seg in segments
        if len(seg.hits) >= min_hits
    ]

    counts = {}
    for _seg, pid in candidates:
        if pid > 0:
            counts[pid] = counts.get(pid, 0) + 1

    kept_segments, kept_pids = [], []
    for seg, pid in candidates:
        if pid > 0 and counts[pid] >= 2:
            kept_segments.append(seg)
            kept_pids.append(pid)

    return kept_segments, kept_pids


# ─── Contrastive Loss ────────────────────────────────────────────────────────


def supcon_loss(
    embeddings: torch.Tensor,
    particle_ids: torch.Tensor,
    temperature: float = 0.1,
) -> torch.Tensor:
    """
    Supervised Contrastive (SupCon / NT-Xent) loss for segment embeddings.

    For each anchor segment, the softmax denominator includes ALL other segments
    in the batch — so every incorrect pairing always contributes gradient,
    with no dead zone.

    For particles with 2 segments (the common looping-track case), this
    reduces to NT-Xent: one positive evaluated against all negatives.

    L(u) = -1/|P(u)| * Σ_{v+∈P(u)} log[ exp(u·v+/τ) / Σ_j≠u exp(u·vj/τ) ]

    Args:
        embeddings:   (N, emb_dim) L2-normalized embedding tensor.
        particle_ids: (N,) long tensor; 0 = noise/invalid (excluded as anchors,
                      but kept in denominator as negatives).
        temperature:  τ — controls cluster tightness. Smaller = tighter clusters.

    Returns:
        Scalar loss tensor.
    """
    N = len(embeddings)
    if N < 2:
        return embeddings.sum() * 0.0

    # Similarity matrix (N, N) scaled by temperature
    sim = (embeddings @ embeddings.T) / temperature

    # Numerical stability: subtract row-wise max before exp
    sim = sim - sim.max(dim=1, keepdim=True)[0].detach()

    # Positive mask: same particle, non-noise, excluding self
    pid = particle_ids
    pos_mask = (pid.unsqueeze(1) == pid.unsqueeze(0)).float()
    pos_mask.fill_diagonal_(0)
    noise_rows = (pid == 0)
    pos_mask[noise_rows] = 0   # noise segments have no valid positives

    # Denominator: all non-self pairs (masked_fill is out-of-place — avoids in-place autograd error)
    self_mask = torch.eye(N, dtype=torch.bool, device=embeddings.device)
    exp_sim = torch.exp(sim).masked_fill(self_mask, 0)
    denom = exp_sim.sum(dim=1, keepdim=True).clamp(min=1e-8)

    log_prob = sim - torch.log(denom)   # (N, N)

    n_pos = pos_mask.sum(dim=1)
    valid = (~noise_rows) & (n_pos > 0)

    if not valid.any():
        return embeddings.sum() * 0.0

    mean_log_prob_pos = (pos_mask * log_prob).sum(dim=1) / n_pos.clamp(min=1)
    return -mean_log_prob_pos[valid].mean()

# ─── Model I/O ───────────────────────────────────────────────────────────────


def load_segment_gnn(model_path: str, device: str = "cpu") -> SegmentGNN:
    """
    Load a trained SegmentGNN from a checkpoint file.

    The checkpoint must be a dict containing 'state_dict' and architecture
    params (hidden_dim, emb_dim, n_layers, dropout, node_scales).
    These are saved automatically by train_mini_GNN_to_match.py.

    Args:
        model_path: Path to the .pt checkpoint.
        device:     Torch device string.

    Returns:
        Loaded model in eval mode on the requested device.
    """
    checkpoint = torch.load(model_path, map_location=device, weights_only=True)

    if not isinstance(checkpoint, dict) or "state_dict" not in checkpoint:
        raise ValueError(
            f"Invalid checkpoint format: {model_path}\n"
            f"Expected a dict with 'state_dict' and architecture params.\n"
            f"Re-train with train_mini_GNN_to_match.py to produce a valid checkpoint."
        )

    model = SegmentGNN(
        node_in_dim=NODE_FEATURE_DIM,
        hidden_dim=checkpoint["hidden_dim"],
        emb_dim=checkpoint["emb_dim"],
        n_layers=checkpoint["n_layers"],
        dropout=checkpoint["dropout"],
        proj_dim=checkpoint.get("proj_dim", None),
        proj_layers=checkpoint.get("proj_layers", 1),
    )
    model.load_state_dict(checkpoint["state_dict"])
    model.node_scales = checkpoint["node_scales"]

    model.to(device)
    model.eval()
    return model
