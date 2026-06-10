#!/usr/bin/env python
"""
Migrate existing graph_constructed_latent graphs to the ACORN-standard
`track_to_edge_map`.

Older builds (build_latent_graphs_fast.py before the fix) stored a custom 2-D
(num_tracks, max_edges) tensor under `track_to_edge_map`, which collides with the
name ACORN's edge_classifier weighting expects (1-D, length = num track_edges).
That mismatch crashes GNN training in construct_weighting with
"IndexError: too many indices for tensor of dimension 1".

This recomputes the correct 1-D map from each graph's stored edge_index +
track_edges via graph_intersection (undirected=True, matching the GNN config) and
overwrites the attribute in place. No model / KNN needed.

Usage:
    python fix_track_to_edge_map.py                 # fix train/val/test under default dir
    python fix_track_to_edge_map.py --input-dir DIR --subsets trainset valset testset
    python fix_track_to_edge_map.py --dry-run       # report only, write nothing
"""
import sys
import os
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "acorn"))

import torch
from tqdm import tqdm

DEFAULT_INPUT = Path(__file__).resolve().parents[1] / "data" / "graph_constructed_latent"


def fast_truth_to_pred(edge_index, track_edges):
    """ACORN-equivalent undirected truth_to_pred via 1-D packed-key lookup.

    For each true edge (track_edges[:, i]) returns its index in edge_index, or -1
    if absent. Matches graph_intersection(..., undirected=True, return_truth_to_pred=True)
    on the matched set, but ~150x faster on CPU (no 2-D torch.unique over ~1M edges).
    """
    N = int(max(edge_index.max(), track_edges.max())) + 1
    pk = edge_index[0].long() * N + edge_index[1].long()
    order = torch.argsort(pk)
    sk = pk[order]

    def lookup(q):
        pos = torch.searchsorted(sk, q).clamp(max=sk.numel() - 1)
        hit = sk[pos] == q
        return torch.where(hit, order[pos], torch.full_like(q, -1)), hit

    fwd = track_edges[0].long() * N + track_edges[1].long()
    rev = track_edges[1].long() * N + track_edges[0].long()
    fi, fh = lookup(fwd)
    ri, rh = lookup(rev)
    return torch.where(fh, fi, torch.where(rh, ri, torch.full_like(fwd, -1)))


def fix_graph(path: Path, dry_run: bool) -> str:
    g = torch.load(path, weights_only=False)
    if not hasattr(g, "track_edges") or not hasattr(g, "edge_index"):
        return "skip-no-truth"

    n_track_edges = g.track_edges.shape[1]
    old = getattr(g, "track_to_edge_map", None)
    already_ok = (
        old is not None and old.dim() == 1 and old.shape[0] == n_track_edges
    )
    if already_ok:
        return "already-1d"

    new_map = fast_truth_to_pred(g.edge_index, g.track_edges)
    assert new_map.shape[0] == n_track_edges, (
        f"{path.name}: map {tuple(new_map.shape)} != track_edges {n_track_edges}"
    )
    if not dry_run:
        g.track_to_edge_map = new_map
        torch.save(g, path)
    return "fixed"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", default=str(DEFAULT_INPUT))
    ap.add_argument("--subsets", nargs="+", default=["trainset", "valset", "testset"])
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    root = Path(args.input_dir)
    totals = {}
    for subset in args.subsets:
        d = root / subset
        if not d.is_dir():
            print(f"[skip] {d} not found")
            continue
        files = sorted(d.glob("*.pyg"))
        counts = {}
        for f in tqdm(files, desc=subset):
            status = fix_graph(f, args.dry_run)
            counts[status] = counts.get(status, 0) + 1
        totals[subset] = counts
        print(f"{subset}: {counts}")

    print("\n=== summary ===")
    for subset, counts in totals.items():
        print(f"  {subset}: {counts}")
    if args.dry_run:
        print("(dry-run — no files written)")


if __name__ == "__main__":
    main()
