#!/usr/bin/env python3
"""
Segment Matching Track Builder — Helix-Based Segment Matching

Builds tracks by fitting helices to GNN-identified segments and matching
segments from the same particle based on helix parameter compatibility
(circle center, radius, pitch).

Algorithm:
    1. Extract segments (CC clustering on GNN edge scores, or ground truth)
    2. Fit circle (Kasa method) and pitch to each segment using ALL hits
    3. Match segments via helix parameter comparison (greedy, highest-score-first)
    4. Attach short segments (1-2 hits) by endpoint proximity
    5. Evaluate and plot

Usage:
    python helix_segmentmatcher_walkthrough.py testset
    python helix_segmentmatcher_walkthrough.py testset --use-gt-segments
    python helix_segmentmatcher_walkthrough.py testset --score-cut 0.8
    python helix_segmentmatcher_walkthrough.py testset --skip-build
"""

import argparse
import sys
from pathlib import Path
from time import process_time, perf_counter

import torch
import yaml
from tqdm import tqdm

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / 'acorn'))
sys.path.insert(0, str(PIPELINE_ROOT))

from acorn.utils.loading_utils import load_datafiles_in_dir

# Import evaluation/plotting from utils module
from low_pt_custom_utils.track_evaluation_utils import (
    run_evaluation,
    save_evaluation_results,
    run_plotting,
)
from low_pt_custom_utils.segment_evaluation_utils import (
    evaluate_matching_efficiency_for_event,
    make_matching_efficiency_summary,
    plot_matching_efficiency,
)

# Import segment matching logic
from low_pt_custom_utils.segment_matching import build_tracks_for_event


# ─── Main Algorithm ─────────────────────────────────────────────────────────

def run_segment_matching(dataset_name, config, score_cut=None, use_gt_segments=None, use_wrangler=None,
                         matching_eff_results=None, data_dir=None, output_name=None):
    """Build tracks for all events in a dataset using helix-based segment matching.

    Args:
        matching_eff_results: If a list is passed, per-event segment matching
            efficiency records are appended to it for downstream plotting.
        data_dir: If provided, override all data paths to subdirectories of this directory.
    """
    if data_dir is not None:
        data_dir = Path(data_dir)
        input_dir = data_dir / config.get('input_dir', 'gnn_stage')
        output_dir = data_dir / config.get('stage_dir', 'track_building')
    else:
        input_dir = Path(config.get('input_dir', 'data/gnn_stage'))
        if not input_dir.is_absolute():
            input_dir = PIPELINE_ROOT / input_dir
        stage_dir = config.get('stage_dir')
        if stage_dir is not None:
            output_dir = Path(stage_dir)
            if not output_dir.is_absolute():
                output_dir = PIPELINE_ROOT / output_dir
        else:
            output_dir = PIPELINE_ROOT / 'data' / 'track_building'
    dataset_output = output_dir / (output_name or dataset_name)
    dataset_output.mkdir(parents=True, exist_ok=True)

    # Override config with CLI args if provided
    if score_cut is not None:
        config['score_cut'] = score_cut
    if use_gt_segments is not None:
        config['use_gt_segments'] = use_gt_segments
    if use_wrangler is not None:
        config['use_wrangler'] = use_wrangler

    actual_score_cut = config.get('score_cut', 0.5)
    actual_use_gt = config.get('use_gt_segments', False)
    actual_use_wrangler = config.get('use_wrangler', False)
    B_field = config.get('B_field', 2.0)

    # Load event files
    data_split = config.get('data_split', [0, 1000, 1000])
    split_map = {'trainset': data_split[0], 'valset': data_split[1], 'testset': data_split[2]}
    num_events = split_map.get(dataset_name, 1000)

    input_paths = load_datafiles_in_dir(str(input_dir), dataset_name, num_events)
    input_paths.sort()

    if actual_use_gt:
        segment_mode = "ground truth (hit_segment_id)"
    elif actual_use_wrangler:
        segment_mode = f"CC + Wrangler (score > {actual_score_cut})"
    else:
        segment_mode = f"CC only (score > {actual_score_cut})"

    print(f"Segment Matching Track Builder (Helix-Based)")
    print(f"  Input:              {input_dir / dataset_name}")
    print(f"  Output:             {dataset_output}")
    print(f"  Events:             {len(input_paths)}")
    print(f"  Segment mode:       {segment_mode}")
    print(f"  B field:            {B_field} T")
    matching = config.get('matching', {})
    print(f"  Max center dist:    {matching.get('max_center_distance', 100.0)} mm")
    print(f"  Sigma center:       {matching.get('sigma_center', 30.0)} mm")
    print(f"  Sigma R:            {matching.get('sigma_R', 0.1)}")
    print()

    # Accumulate statistics
    total_stats = {
        "n_segments": 0,
        "n_wrangler_splits": 0,
        "n_good_fits": 0,
        "n_poor_fits": 0,
        "n_no_fits": 0,
        "n_matched_pairs": 0,
        "n_complete_tracks": 0,
        "n_no_match": 0,
        "n_standalone": 0,
        "n_total_tracks": 0,
        "n_assigned_hits": 0,
        "n_unassigned_hits": 0,
    }
    total_time = 0.0

    for event_idx, event_path in enumerate(tqdm(input_paths, desc=f"Building tracks for {dataset_name}")):
        t_event = perf_counter()
        graph = torch.load(event_path, map_location="cpu", weights_only=False)

        labels, event_stats, matching_info = build_tracks_for_event(
            graph, config, return_matching_info=True
        )

        # Accumulate segment matching efficiency records if requested
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

        # Accumulate statistics
        for key in total_stats:
            total_stats[key] += event_stats[key]

        # Save
        event_id = graph.event_id
        if isinstance(event_id, list):
            event_id = event_id[0]
        torch.save(graph, dataset_output / f"event{event_id}.pyg")

    # Print summary
    n_events = len(input_paths)
    print(f"\n{'='*70}")
    print(f"SEGMENT MATCHING SUMMARY - {dataset_name.upper()}")
    print(f"{'='*70}")

    print(f"\nSegment Statistics (totals across {n_events} events):")
    print(f"  Total segments:        {total_stats['n_segments']:6d}")
    if actual_use_wrangler and not actual_use_gt:
        print(f"  Wrangler splits:       {total_stats['n_wrangler_splits']:6d}  "
              f"(CC clusters disentangled into >1 segment)")
    print(f"  Good fits (3+ hits):   {total_stats['n_good_fits']:6d}")
    print(f"  Poor fits (2 hits):    {total_stats['n_poor_fits']:6d}")
    print(f"  No fits (1 hit):       {total_stats['n_no_fits']:6d}")

    print(f"\nMatching Results:")
    print(f"  Matched pairs:         {total_stats['n_matched_pairs']:6d}  (2 segments → 1 track)")
    print(f"  Standalone (exit):     {total_stats['n_complete_tracks']:6d}  (outer_r ≥ threshold, complete tracks)")
    print(f"  Standalone (no match): {total_stats['n_no_match']:6d}  (loop segment, no partner found)")
    print(f"  Total tracks:          {total_stats['n_total_tracks']:6d}")

    print(f"\nHit Assignment:")
    print(f"  Assigned:              {total_stats['n_assigned_hits']:6d}")
    print(f"  Unassigned:            {total_stats['n_unassigned_hits']:6d}")

    print(f"\nTiming:")
    print(f"  Total time:            {total_time:.2f}s")
    if n_events > 0:
        print(f"  Average per event:     {total_time/n_events:.4f}s")
        print(f"  Events per second:     {n_events/total_time:.1f}" if total_time > 0 else "  Events per second:     N/A")

    print(f"\n{'='*70}")
    print(f"Tracks saved to: {dataset_output}")
    print(f"{'='*70}")


# ─── Script Entry Point ─────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Helix-based segment matching track builder',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python helix_segmentmatcher_walkthrough.py testset
  python helix_segmentmatcher_walkthrough.py testset --use-gt-segments
  python helix_segmentmatcher_walkthrough.py testset --score-cut 0.8
  python helix_segmentmatcher_walkthrough.py testset --skip-build
        """
    )
    parser.add_argument(
        'dataset',
        type=str,
        choices=['trainset', 'valset', 'testset'],
        help='Dataset to process'
    )
    parser.add_argument(
        '--config',
        type=str,
        default=None,
        help='Path to config file (default: acorn_configs/track_building_stage_(3)/helix_segmentmatcher_walkthrough.yaml)'
    )
    parser.add_argument(
        '--score-cut',
        type=float,
        default=None,
        help='Override score cut threshold for CC clustering'
    )
    parser.add_argument(
        '--use-gt-segments',
        action='store_true',
        default=None,
        help='Use ground truth segment labels instead of CC clusters'
    )
    parser.add_argument(
        '--no-wrangler',
        action='store_true',
        help='Disable Wrangler: use plain CC without walk-through disentanglement'
    )
    parser.add_argument(
        '--skip-build',
        action='store_true',
        help='Skip track building, only re-evaluate'
    )
    parser.add_argument(
        '--data-dir',
        type=str,
        default=None,
        help='Override all data paths (gnn_stage, track_building, track_evaluation, visuals) to subdirs of this directory'
    )
    parser.add_argument(
        '--ood-config',
        type=str,
        default=None,
        help='Path to ood_experiment.yaml; overrides track_building and evaluation hyperparameters'
    )


    args = parser.parse_args()

    # Load config
    if args.config is None:
        config_file = PIPELINE_ROOT / 'acorn_configs' / 'track_building_stage_(3)' / 'helix_segmentmatcher_walkthrough.yaml'
    else:
        config_file = Path(args.config)

    if not config_file.exists():
        raise FileNotFoundError(f"Config file not found: {config_file}")

    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    config['dataset'] = args.dataset

    # Merge OOD experiment config overrides
    if args.ood_config is not None:
        with open(args.ood_config) as f:
            ood = yaml.safe_load(f)
        config['data_split'] = [0, 0, ood['n_events']]
        config.update(ood.get('track_building', {}))
        config.update(ood.get('evaluation', {}))

    output_dir_suffix = config.get('output_dir') or None
    output_dataset = f"{args.dataset}_{output_dir_suffix}" if output_dir_suffix else args.dataset

    # data_set_path (config) acts as a default for --data-dir: all stage dirs
    # (input_dir, stage_dir, track_evaluation, visuals) live under it.
    data_root = args.data_dir or config.get('data_set_path')
    data_dir = Path(data_root) if data_root else None
    if data_dir is not None and not data_dir.is_absolute():
        data_dir = PIPELINE_ROOT / data_dir
    if data_dir is not None:
        eval_output_dir = data_dir / 'track_evaluation' / output_dataset
        plot_output_dir = data_dir / 'visuals' / 'track_metrics' / output_dataset
    else:
        eval_output_dir = PIPELINE_ROOT / 'data' / 'track_evaluation' / output_dataset
        plot_output_dir = PIPELINE_ROOT / 'data' / 'visuals' / 'track_metrics' / output_dataset

    # ── Step 1: Track Building ──────────────────────────────────────────
    matching_eff_results = []
    if not args.skip_build:
        print("=" * 70)
        print("STEP 1: HELIX-BASED SEGMENT MATCHING")
        print("=" * 70)
        print()

        run_segment_matching(
            args.dataset, config,
            score_cut=args.score_cut,
            use_gt_segments=args.use_gt_segments,
            use_wrangler=(False if args.no_wrangler else None),
            matching_eff_results=matching_eff_results,
            data_dir=data_dir,
            output_name=output_dataset,
        )
        print()
    else:
        print("Skipping track building (--skip-build)")
        print()

    # ── Step 2: Evaluation ──────────────────────────────────────────────
    print("=" * 70)
    print(f"STEP 2: EVALUATING {args.dataset.upper()}")
    print("=" * 70)
    print()

    # pre_build_tracks_dir (only with --skip-build): read prebuilt tracks from this
    # dataset dir (relative to data_set_path), decoupled from output_dir naming.
    prebuilt = config.get('pre_build_tracks_dir') if args.skip_build else None
    if prebuilt:
        prebuilt = Path(prebuilt)
        if not prebuilt.is_absolute():
            prebuilt = (data_dir / prebuilt) if data_dir is not None else (PIPELINE_ROOT / prebuilt)
        config['input_dir'] = str(prebuilt.parent)
        read_dataset = prebuilt.name
    else:
        if data_dir is not None:
            config['input_dir'] = str(data_dir / config.get('stage_dir', 'track_building'))
        else:
            config['input_dir'] = str(PIPELINE_ROOT / 'data' / 'track_building')
        read_dataset = output_dataset
    evaluated_events, summary, summary_text = run_evaluation(read_dataset, config)
    save_evaluation_results(evaluated_events, summary, summary_text, output_dataset, eval_output_dir)

    print()
    print(summary_text)

    # ── Step 3: Plotting ────────────────────────────────────────────────
    print("=" * 70)
    print("STEP 3: PLOTTING METRICS")
    print("=" * 70)
    print()

    run_plotting(evaluated_events, summary, output_dataset, plot_output_dir, config.get('plots', {}))

    # ── Step 3b: Segment Matching Efficiency ────────────────────────────
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
            config.get('plots', {}).get('segment_evaluation', {}),
        )
    else:
        print("\n(No segment matching efficiency data — run without --skip-build to compute.)")

    # ── Done ────────────────────────────────────────────────────────────
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
