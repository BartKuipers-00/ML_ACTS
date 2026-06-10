#!/usr/bin/env python3
"""
Run mini-GNN segment data mining for a specific chunk of events.

Each job processes a disjoint subset of events from each split and saves
per-event .pyg files to data/track_building/mini_gnn_segments_chunk{N:02d}/<split>/.

After all jobs complete, run:
    python combine_mini_gnn_chunks.py

Usage:
    python run_mini_gnn_mining_chunk.py --chunk 0 --total-chunks 50
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
from low_pt_custom_utils.mini_gnn_segment_embedding import (
    segment_to_pyg,
    filter_paired_segments,
)


def get_chunk_paths(all_paths, chunk, total_chunks):
    """Return the disjoint slice of paths assigned to this chunk."""
    n = len(all_paths)
    events_per_chunk = n // total_chunks
    remainder = n % total_chunks

    if chunk < remainder:
        start = chunk * (events_per_chunk + 1)
        end = start + events_per_chunk + 1
    else:
        start = remainder * (events_per_chunk + 1) + (chunk - remainder) * events_per_chunk
        end = start + events_per_chunk

    return all_paths[start:end]


def mine_paths(paths, score_cut, use_ground_truth, node_scales, output_dir, label=""):
    """Process events and save precomputed segments to output_dir."""
    output_dir.mkdir(parents=True, exist_ok=True)
    total_segs = 0
    skipped = 0

    pbar = tqdm(paths, desc=f"  Mining {label}", unit="ev")
    for event_path in pbar:
        graph = torch.load(event_path, map_location="cpu", weights_only=False)
        if use_ground_truth:
            segments = extract_segments_from_ground_truth(graph)
        else:
            segments = extract_segments_from_cc(graph, score_cut)

        # Pairs-only + min-hits >=2 filter (see mini_gnn_data_mining.py docstring)
        kept_segs, particle_ids = filter_paired_segments(segments or [], graph, min_hits=2)
        seg_data_list = [segment_to_pyg(seg, graph, node_scales=node_scales) for seg in kept_segs]

        if len(seg_data_list) < 2:
            skipped += 1
            continue

        out_path = output_dir / Path(event_path).name
        torch.save({"seg_data_list": seg_data_list, "particle_ids": particle_ids}, out_path)
        total_segs += len(seg_data_list)
        pbar.set_postfix(segs=total_segs)

    n_saved = len(paths) - skipped
    print(f"  {label}: {n_saved}/{len(paths)} events saved "
          f"({skipped} skipped with <2 segments), {total_segs} segments")


def main():
    parser = argparse.ArgumentParser(description="Mine mini-GNN segments for one chunk of events")
    parser.add_argument("--chunk", type=int, required=True, help="Chunk index (0-indexed)")
    parser.add_argument("--total-chunks", type=int, required=True, help="Total number of chunks")
    parser.add_argument("--config", type=str, default=None, help="Path to config YAML")
    args = parser.parse_args()

    if args.config is None:
        config_file = (
            PIPELINE_ROOT / "acorn_configs" / "track_building_stage_(3)"
            / "mini_gnn_data_mining.yaml"
        )
    else:
        config_file = Path(args.config)

    with open(config_file) as f:
        config = yaml.safe_load(f)

    print("=" * 65)
    print(f"MINI-GNN DATA MINING — Chunk {args.chunk}/{args.total_chunks}")
    print("=" * 65)
    print(f"Config: {config_file}")
    print(f"Node:   {__import__('socket').gethostname()}")
    print()

    input_dir = Path(config["input_dir"])
    if not input_dir.is_absolute():
        input_dir = PIPELINE_ROOT / input_dir

    chunk_output_dir = PIPELINE_ROOT / f"data/track_building/mini_gnn_segments_chunk{args.chunk:02d}"
    print(f"Input:  {input_dir}")
    print(f"Output: {chunk_output_dir}")
    print()

    score_cut = config["score_cut"]
    use_ground_truth = config["use_ground_truth"]
    node_scales = config["node_scales"]
    data_split = config["data_split"]
    print(f"GroundTruth: {use_ground_truth}")

    splits = [("trainset", data_split[0]), ("valset", data_split[1])]
    if len(data_split) > 2 and data_split[2] > 0:
        splits.append(("testset", data_split[2]))

    t0 = perf_counter()
    for split_name, n_events in splits:
        if n_events == 0:
            continue

        all_paths = load_datafiles_in_dir(str(input_dir), split_name, n_events)
        all_paths.sort()
        chunk_paths = get_chunk_paths(all_paths, args.chunk, args.total_chunks)

        print(f"{split_name}: {len(chunk_paths)} / {len(all_paths)} events in this chunk")
        if not chunk_paths:
            continue

        mine_paths(chunk_paths, score_cut, use_ground_truth, node_scales,
                   chunk_output_dir / split_name, label=split_name)
        print()

    print(f"Chunk {args.chunk} complete in {perf_counter() - t0:.1f}s.")


if __name__ == "__main__":
    main()
