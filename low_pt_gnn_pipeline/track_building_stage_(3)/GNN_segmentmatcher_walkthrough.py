#!/usr/bin/env python3
"""
GNN Segment Matching Track Builder with Wrangler Disentanglement

Same pipeline as GNN_segment_matching_track_builder.py, but adds an optional
Wrangler (walk-through) step between CC clustering and GNN segment matching.

Algorithm:
    1. CC clustering on GNN edge scores (score > score_cut)
    2. [Optional, use_wrangler=true] Wrangler: convert each CC subgraph to a
       directed graph oriented radially outward, then resolve branching via
       greedy walk from source nodes (in-degree 0). Produces disentangled
       path segments.
    3. GNN embedding + cosine similarity matching of segments
    4. Evaluate and plot

Usage:
    python GNN_segmentmatcher_walkthrough.py testset
    python GNN_segmentmatcher_walkthrough.py testset --use-gt-segments
    python GNN_segmentmatcher_walkthrough.py testset --no-wrangler
    python GNN_segmentmatcher_walkthrough.py testset --score-cut 0.8
    python GNN_segmentmatcher_walkthrough.py testset --skip-build
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
sys.path.insert(0, str(WORKSPACE_ROOT / 'acorn'))
sys.path.insert(0, str(PIPELINE_ROOT))

from acorn.utils.loading_utils import load_datafiles_in_dir

from low_pt_custom_utils.track_evaluation_utils import (
    run_evaluation,
    save_evaluation_results,
    run_plotting,
)
from low_pt_custom_utils.mini_gnn_segment_embedding import load_segment_gnn
from low_pt_custom_utils.segmentGNN import build_tracks_for_event_gnn
from low_pt_custom_utils.segment_evaluation_utils import (
    evaluate_matching_efficiency_for_event,
    make_matching_efficiency_summary,
    plot_matching_efficiency,
)


# ─── Main Algorithm ─────────────────────────────────────────────────────────


def run_gnn_segment_matching(dataset_name, config, model, device,
                              score_cut=None, use_gt_segments=None,
                              use_wrangler=None, matching_eff_results=None, output_name=None):
    """Build tracks for all events in a dataset using GNN-scored segment matching."""
    input_dir = Path(config["input_dir"])
    if not input_dir.is_absolute():
        input_dir = PIPELINE_ROOT / input_dir

    stage_dir = config.get("stage_dir")
    if stage_dir is not None:
        output_dir = Path(stage_dir)
        if not output_dir.is_absolute():
            output_dir = PIPELINE_ROOT / output_dir
    else:
        output_dir = PIPELINE_ROOT / "data" / "track_building"
    dataset_output = output_dir / (output_name or dataset_name)
    dataset_output.mkdir(parents=True, exist_ok=True)

    if score_cut is not None:
        config["score_cut"] = score_cut
    if use_gt_segments is not None:
        config["use_gt_segments"] = use_gt_segments
    if use_wrangler is not None:
        config["use_wrangler"] = use_wrangler

    actual_score_cut = config["score_cut"]
    actual_use_gt    = config["use_gt_segments"]
    actual_wrangler  = config.get("use_wrangler", True)
    gnn_config       = config["gnn"]
    wrangler_cfg     = config.get("wrangler", {})

    data_split = config["data_split"]
    split_map  = {"trainset": data_split[0], "valset": data_split[1], "testset": data_split[2]}
    num_events = split_map[dataset_name]

    input_paths = load_datafiles_in_dir(str(input_dir), dataset_name, num_events)
    input_paths.sort()

    if actual_use_gt:
        segment_mode = "ground truth (hit_segment_id)"
    elif actual_wrangler:
        segment_mode = f"CC + Wrangler (score > {actual_score_cut})"
    else:
        segment_mode = f"CC only (score > {actual_score_cut})"

    geo_on = gnn_config["geo_cut_enabled"]

    print(f"GNN Segment Matching Track Builder (with Wrangler)")
    print(f"  Input:              {input_dir / dataset_name}")
    print(f"  Output:             {dataset_output}")
    print(f"  Events:             {len(input_paths)}")
    print(f"  Segment mode:       {segment_mode}")
    print(f"  Cos sim threshold:  {gnn_config['cos_sim_threshold']}")
    print(f"  Outer R threshold:  {gnn_config['outer_r_threshold']} mm")
    print(f"  Geometric cut:      {'ON' if geo_on else 'OFF'}", end="")
    if geo_on:
        print(f"  (min_3d >= {gnn_config['geo_cut_min_3d_mm']:.0f} mm"
              f"  at r > {gnn_config['geo_cut_r_high_mm']:.0f} mm)", end="")
    print()
    print(f"  Device:             {device}")
    print()

    total_stats = {
        "n_segments": 0,
        "n_wrangler_splits": 0,
        "n_matched_pairs": 0,
        "n_complete_tracks": 0,
        "n_no_match": 0,
        "n_standalone": 0,
        "n_total_tracks": 0,
        "n_assigned_hits": 0,
        "n_unassigned_hits": 0,
    }
    total_time = 0.0

    for event_path in tqdm(input_paths, desc=f"Building tracks for {dataset_name}"):
        t_event = perf_counter()
        graph = torch.load(event_path, map_location="cpu", weights_only=False)

        labels, event_stats, matching_info = build_tracks_for_event_gnn(
            graph, config, model, device, return_matching_info=True
        )

        if matching_eff_results is not None:
            fiducial = config.get("target_segments", {"pt": [0.1, float("inf")], "nhits": [3, float("inf")]})
            event_matching_eff = evaluate_matching_efficiency_for_event(
                graph,
                matching_info["segments"],
                matching_info["matched_tracks"],
                matching_info["unmatched"],
                fiducial,
            )
            matching_eff_results.extend(event_matching_eff)

        graph.hit_track_labels = labels
        graph.time_taken = perf_counter() - t_event
        total_time += graph.time_taken

        for key in total_stats:
            total_stats[key] += event_stats[key]

        event_id = graph.event_id
        if isinstance(event_id, list):
            event_id = event_id[0]
        torch.save(graph, dataset_output / f"event{event_id}.pyg")

    n_events = len(input_paths)
    print(f"\n{'='*70}")
    print(f"GNN SEGMENT MATCHING (WRANGLER) SUMMARY - {dataset_name.upper()}")
    print(f"{'='*70}")

    print(f"\nSegment Statistics (totals across {n_events} events):")
    print(f"  Total segments:        {total_stats['n_segments']:6d}")
    if actual_wrangler and not actual_use_gt:
        print(f"  Wrangler splits:       {total_stats['n_wrangler_splits']:6d}  "
              f"(CC clusters disentangled into >1 segment)")

    print(f"\nMatching Results:")
    print(f"  Matched pairs:         {total_stats['n_matched_pairs']:6d}  (2 segments → 1 track)")
    print(f"  Standalone (exit):     {total_stats['n_complete_tracks']:6d}  (outer_r >= threshold)")
    print(f"  Standalone (no match): {total_stats['n_no_match']:6d}  (no GNN match found)")
    print(f"  Total tracks:          {total_stats['n_total_tracks']:6d}")

    print(f"\nHit Assignment:")
    print(f"  Assigned:              {total_stats['n_assigned_hits']:6d}")
    print(f"  Unassigned:            {total_stats['n_unassigned_hits']:6d}")

    print(f"\nTiming:")
    print(f"  Total time:            {total_time:.2f}s")
    print(f"  Average per event:     {total_time/n_events:.4f}s")
    print(f"  Events per second:     {n_events/total_time:.1f}")

    print(f"\n{'='*70}")
    print(f"Tracks saved to: {dataset_output}")
    print(f"{'='*70}")


# ─── Script Entry Point ─────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="GNN segment embedding track builder with Wrangler disentanglement",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python GNN_segmentmatcher_walkthrough.py testset
  python GNN_segmentmatcher_walkthrough.py testset --use-gt-segments
  python GNN_segmentmatcher_walkthrough.py testset --no-wrangler
  python GNN_segmentmatcher_walkthrough.py testset --score-cut 0.8
  python GNN_segmentmatcher_walkthrough.py testset --skip-build
        """,
    )
    parser.add_argument(
        "dataset",
        type=str,
        choices=["trainset", "valset", "testset"],
        help="Dataset to process",
    )
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to config file (default: acorn_configs/track_building_stage_(3)/gnn_segmentmatcher_walkthrough.yaml)",
    )
    parser.add_argument(
        "--score-cut",
        type=float,
        default=None,
        help="Override CC edge score cut threshold",
    )
    parser.add_argument(
        "--use-gt-segments",
        action="store_true",
        default=None,
        help="Use ground truth segment labels instead of CC clusters",
    )
    parser.add_argument(
        "--no-wrangler",
        action="store_true",
        help="Disable Wrangler: use plain CC without disentanglement",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Skip track building, only re-evaluate existing output",
    )
    parser.add_argument(
        "--device",
        type=str,
        default=None,
        help="Torch device (default: cuda if available, else cpu)",
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default=None,
        help="Override all data paths to subdirectories of this directory",
    )
    parser.add_argument(
        "--ood-config",
        type=str,
        default=None,
        help="Path to ood_experiment.yaml; overrides segment_gnn and evaluation hyperparameters",
    )


    args = parser.parse_args()

    # Load config
    if args.config is None:
        config_file = (
            PIPELINE_ROOT
            / "acorn_configs"
            / "track_building_stage_(3)"
            / "gnn_segmentmatcher_walkthrough.yaml"
        )
    else:
        config_file = Path(args.config)

    if not config_file.exists():
        raise FileNotFoundError(f"Config file not found: {config_file}")

    with open(config_file, "r") as f:
        config = yaml.safe_load(f)
    config["dataset"] = args.dataset

    # Merge OOD experiment config overrides
    if args.ood_config is not None:
        with open(args.ood_config) as f:
            ood = yaml.safe_load(f)
        config["data_split"] = [0, 0, ood["n_events"]]
        config.update(ood.get("evaluation", {}))
        sgnn = ood.get("segment_gnn", {})
        if sgnn:
            config.setdefault("gnn", {}).update(sgnn)
        tb = ood.get("track_building", {})
        for key in ("score_cut", "use_gt_segments", "use_wrangler"):
            if key in tb:
                config[key] = tb[key]

    device = args.device or ("cuda" if torch.cuda.is_available() else "cpu")

    output_dir_suffix = config.get('output_dir') or None
    output_dataset = f"{args.dataset}_{output_dir_suffix}" if output_dir_suffix else args.dataset

    data_dir = Path(args.data_dir) if args.data_dir else None
    if data_dir is not None:
        config["input_dir"] = str(data_dir / "gnn_stage")
        config["stage_dir"] = str(data_dir / "track_building")
        eval_output_dir = data_dir / "track_evaluation" / output_dataset
        plot_output_dir = data_dir / "visuals" / "track_metrics" / output_dataset
    else:
        eval_output_dir = PIPELINE_ROOT / "data" / "track_evaluation" / output_dataset
        plot_output_dir = PIPELINE_ROOT / "data" / "visuals" / "track_metrics" / output_dataset

    # ── Step 1: Track Building ──────────────────────────────────────────────
    matching_eff_results = []
    if not args.skip_build:
        print("=" * 70)
        print("STEP 1: GNN SEGMENT MATCHING WITH WRANGLER")
        print("=" * 70)
        print()

        gnn_config = config["gnn"]
        model_path = gnn_config["model_path"]
        if not Path(model_path).is_absolute():
            model_path = str(PIPELINE_ROOT / model_path)

        if not Path(model_path).exists():
            raise FileNotFoundError(
                f"Trained GNN not found at: {model_path}\n"
                f"Run train_mini_GNN_to_match.py first."
            )

        print(f"Loading mini-GNN from: {model_path}")
        model = load_segment_gnn(model_path, device=device)
        print()

        run_gnn_segment_matching(
            args.dataset,
            config,
            model,
            device,
            score_cut=args.score_cut,
            use_gt_segments=args.use_gt_segments,
            use_wrangler=(False if args.no_wrangler else None),
            matching_eff_results=matching_eff_results,
            output_name=output_dataset,
        )
        print()
    else:
        print("Skipping track building (--skip-build)")
        print()

    # ── Step 2: Evaluation ─────────────────────────────────────────────────
    print("=" * 70)
    print(f"STEP 2: EVALUATING {args.dataset.upper()}")
    print("=" * 70)
    print()

    if data_dir is not None:
        config["input_dir"] = str(data_dir / "track_building")
    else:
        config["input_dir"] = str(PIPELINE_ROOT / "data" / "track_building")
    evaluated_events, summary, summary_text = run_evaluation(output_dataset, config)
    save_evaluation_results(evaluated_events, summary, summary_text, output_dataset, eval_output_dir)

    print()
    print(summary_text)

    # ── Step 3: Plotting ───────────────────────────────────────────────────
    print("=" * 70)
    print("STEP 3: PLOTTING METRICS")
    print("=" * 70)
    print()

    run_plotting(evaluated_events, summary, output_dataset, plot_output_dir, config["plots"])

    # ── Step 3b: Segment Matching Efficiency ───────────────────────────────
    if matching_eff_results:
        print()
        print("=" * 70)
        print("STEP 3b: SEGMENT MATCHING EFFICIENCY")
        print("=" * 70)
        print()
        print(make_matching_efficiency_summary(matching_eff_results, n_events=None))
        print()
        plot_matching_efficiency(
            matching_eff_results,
            plot_output_dir,
            output_dataset,
            config.get("plots", {}).get("segment_evaluation", {}),
        )
    else:
        print("\n(No segment matching efficiency data — run without --skip-build to compute.)")

    # ── Done ───────────────────────────────────────────────────────────────
    print()
    print("=" * 70)
    print("COMPLETE!")
    print("=" * 70)
    print(f"  Tracks:     data/track_building/")
    print(f"  Evaluation: {eval_output_dir.relative_to(PIPELINE_ROOT)}/")
    print(f"  Plots:      {plot_output_dir.relative_to(PIPELINE_ROOT)}/")
    print("=" * 70)


if __name__ == "__main__":
    main()
