#!/usr/bin/env python3
"""
Mini-GNN Segment Data Mining

Mines per-segment PyG graphs from GNN-inferred event graphs (analogous to
mlp_data_mining.py which mines segment *pair* features for the MLP matcher).

For each event, extracts CC segments, converts each segment to a small PyG
Data object (hits as nodes with x/y/z/r features, fully-connected edges),
and records the majority particle ID. Saves one .pyg file per event.

Segment filter (pairs-only + min-hits >=2, applied jointly):
    Keep a segment iff it has >=2 hits AND its particle has >=2 such (>=2-hit)
    segments. Applied jointly so dropping 1-hit arcs never orphans a good
    multi-hit partner. Drops unpaired (singleton) and no-curvature 1-hit
    segments — only matchable pairs are mined.

This is CPU-only (no GPU needed). Run once, then train with
train_mini_GNN_to_match.py which loads the mined data directly onto GPU.

Output: data/track_building/mini_gnn_segments/<split>/ — one .pyg file per event,
each containing {"seg_data_list": [...], "particle_ids": [...]}.

Usage:
    python mini_gnn_data_mining.py
    python mini_gnn_data_mining.py --config path/to/config.yaml
"""

import argparse
import sys
from pathlib import Path
from time import perf_counter

import torch
import yaml
from tqdm import tqdm

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / "acorn"))
sys.path.insert(0, str(PIPELINE_ROOT))

from acorn.utils.loading_utils import load_datafiles_in_dir
from low_pt_custom_utils.segment_matching import (
    extract_segments_from_cc,
    extract_segments_from_ground_truth,
)
from low_pt_custom_utils.wrangler_utils import extract_segments_with_wrangler
from low_pt_custom_utils.mini_gnn_segment_embedding import (
    segment_to_pyg,
    filter_paired_segments,
)


def process_and_save(event_paths, score_cut, use_wrangler, use_ground_truth, node_scales, output_dir, label=""):
    """Process events and save precomputed segments to output_dir."""
    output_dir.mkdir(parents=True, exist_ok=True)
    total_segs = 0
    skipped = 0

    pbar = tqdm(event_paths, desc=f"Precompute {label}", unit="ev")
    for event_path in pbar:
        graph = torch.load(event_path, map_location="cpu", weights_only=False)
        if use_ground_truth:
            segments = extract_segments_from_ground_truth(graph)
        elif use_wrangler:
            segments, _ = extract_segments_with_wrangler(graph, score_cut)
        else:
            segments = extract_segments_from_cc(graph, score_cut)

        # Pairs-only + min-hits >=2 filter (see module docstring)
        kept_segs, particle_ids = filter_paired_segments(segments or [], graph, min_hits=2)
        seg_data_list = [segment_to_pyg(seg, graph, node_scales=node_scales) for seg in kept_segs]

        if len(seg_data_list) < 2:
            skipped += 1
            continue

        # Save with same filename as source event
        out_path = output_dir / Path(event_path).name
        torch.save({"seg_data_list": seg_data_list, "particle_ids": particle_ids}, out_path)
        total_segs += len(seg_data_list)
        pbar.set_postfix(segs=total_segs)

    n_saved = len(event_paths) - skipped
    print(f"  {label}: {n_saved}/{len(event_paths)} events saved "
          f"({skipped} skipped with <2 segments)")
    print(f"  Total segments: {total_segs} (avg {total_segs / max(1, n_saved):.1f}/event)")
    print(f"  Output: {output_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description="Mine segment PyG graphs for mini-GNN training (CPU only)",
    )
    parser.add_argument("--config", type=str, default=None)
    args = parser.parse_args()

    if args.config is None:
        config_file = (
            PIPELINE_ROOT / "acorn_configs" / "track_building_stage_(3)"
            / "mini_gnn_data_mining.yaml"
        )
    else:
        config_file = Path(args.config)

    if not config_file.exists():
        raise FileNotFoundError(f"Config file not found: {config_file}")

    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    input_dir = Path(config["input_dir"])
    if not input_dir.is_absolute():
        input_dir = PIPELINE_ROOT / input_dir

    output_dir = Path(config["precomputed_dir"])
    if not output_dir.is_absolute():
        output_dir = PIPELINE_ROOT / output_dir

    score_cut = config["score_cut"]
    use_wrangler = config["use_wrangler"]
    use_ground_truth = config["use_ground_truth"]
    node_scales = config["node_scales"]
    data_split = config["data_split"]

    print("=" * 65)
    print("MINI-GNN SEGMENT DATA MINING")
    print("=" * 65)
    print(f"Config:      {config_file}")
    print(f"Input:       {input_dir}")
    print(f"Output:      {output_dir}")
    print(f"Score cut:   {score_cut}")
    print(f"Wrangler:    {use_wrangler}")
    print(f"GroundTruth: {use_ground_truth}")
    print()

    t0 = perf_counter()

    splits = [("trainset", data_split[0]), ("valset", data_split[1])]
    if len(data_split) > 2:
        splits.append(("testset", data_split[2]))

    for split, n_events in splits:
        if n_events <= 0:
            continue
        paths = load_datafiles_in_dir(str(input_dir), split, n_events)
        paths.sort()
        print(f"\n{split}: {len(paths)} events")
        process_and_save(paths, score_cut, use_wrangler, use_ground_truth, node_scales, output_dir / split, label=split)

    print(f"\nDone in {perf_counter() - t0:.1f}s")


if __name__ == "__main__":
    main()
