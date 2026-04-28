#!/usr/bin/env python3
"""
Segment-level evaluation: purity and completeness of CC/Wrangler segments
before the matching stage.

For each true segment (particle, segment_id) passing fiducial cuts, the
proposed CC/Wrangler segment with the most shared hits is found.
Completeness and purity are then reported as a function of pT.

  completeness = shared_hits / n_true_hits
  purity       = shared_hits / n_proposed_hits

No matching threshold: if no proposed segment shares any hits, both = 0.

Usage:
    python evaluate_segments.py testset
    python evaluate_segments.py testset --use-gt-segments   # GT as both true AND proposed
    python evaluate_segments.py testset --no-wrangler       # CC only (skip wrangler)
    python evaluate_segments.py testset --score-cut 0.6
"""

import argparse
import sys
from pathlib import Path

import torch
import yaml
from tqdm import tqdm

SCRIPT_DIR    = Path(__file__).resolve().parent
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
from low_pt_custom_utils.segment_evaluation_utils import (
    evaluate_segments_for_event,
    make_segment_evaluation_summary,
    plot_segment_evaluation,
)


def run_segment_evaluation(dataset_name: str, config: dict, args) -> None:
    """Evaluate segment purity/completeness over a full dataset."""

    input_dir = Path(config["input_dir"])
    if not input_dir.is_absolute():
        input_dir = PIPELINE_ROOT / input_dir

    # Override config from CLI
    score_cut       = args.score_cut if args.score_cut is not None else config.get("score_cut", 0.4)
    use_gt_segments = args.use_gt_segments or config.get("use_gt_segments", False)
    use_wrangler    = (not args.no_wrangler) and config.get("use_wrangler", True) and not use_gt_segments

    fiducial = config.get("target_segments", {
        "pt":    [0.1, float("inf")],
        "nhits": [3,   float("inf")],
    })
    plot_cfg = config.get("plots", {}).get("segment_evaluation", {})

    data_split = config["data_split"]
    split_map  = {"trainset": data_split[0], "valset": data_split[1], "testset": data_split[2]}
    num_events = split_map[dataset_name]

    input_paths = load_datafiles_in_dir(str(input_dir), dataset_name, num_events)
    input_paths.sort()

    if use_gt_segments:
        mode_str = "ground truth segments (hit_segment_id)"
    elif use_wrangler:
        mode_str = f"CC + Wrangler (score > {score_cut})"
    else:
        mode_str = f"CC only (score > {score_cut})"

    print(f"Segment Evaluation")
    print(f"  Dataset:         {dataset_name}  ({len(input_paths)} events)")
    print(f"  Input:           {input_dir / dataset_name}")
    print(f"  Proposed segs:   {mode_str}")
    print(f"  Fiducial pt:     {fiducial['pt']}")
    print(f"  Fiducial nhits:  {fiducial['nhits']}")
    print()

    all_results = []

    # ── Diagnostic: check first graph for required attributes ───────────────
    first_graph = torch.load(input_paths[0], map_location="cpu", weights_only=False)
    keys = list(first_graph.keys()) if hasattr(first_graph, "keys") else dir(first_graph)
    has_segment_id = hasattr(first_graph, "hit_segment_id")
    print(f"First graph attributes: {[k for k in first_graph.keys()]}")
    if not has_segment_id:
        print()
        print("WARNING: hit_segment_id not found in graphs.")
        print("  True segment evaluation requires hit_segment_id on each node.")
        print("  It is listed in convert_csv_to_pyg_sets.yaml hit_features, but")
        print("  may have been dropped by a downstream stage (graph construction")
        print("  or GNN inference). Re-check those stages' feature propagation.")
        print()

    for event_path in tqdm(input_paths, desc="Evaluating segments"):
        graph = torch.load(event_path, map_location="cpu", weights_only=False)

        # Extract proposed segments
        if use_gt_segments:
            proposed = extract_segments_from_ground_truth(graph)
        elif use_wrangler:
            proposed, _ = extract_segments_with_wrangler(graph, score_cut)
        else:
            proposed = extract_segments_from_cc(graph, score_cut)

        # Evaluate against true segments
        event_results = evaluate_segments_for_event(graph, proposed, fiducial)
        all_results.extend(event_results)

    # ── Summary ─────────────────────────────────────────────────────────────
    print()
    print("=" * 60)
    print("SEGMENT EVALUATION SUMMARY")
    print("=" * 60)
    print(make_segment_evaluation_summary(all_results, len(input_paths)))
    print("=" * 60)

    # ── Plot ─────────────────────────────────────────────────────────────────
    plot_out = PIPELINE_ROOT / "data" / "visuals" / "segment_evaluation" / dataset_name
    plot_segment_evaluation(all_results, plot_out, dataset_name, plot_cfg)

    print()
    print(f"Plots saved to: {plot_out.relative_to(PIPELINE_ROOT)}/")


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate CC/Wrangler segment purity and completeness vs pT",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python evaluate_segments.py testset
  python evaluate_segments.py testset --no-wrangler
  python evaluate_segments.py testset --use-gt-segments
  python evaluate_segments.py testset --score-cut 0.6
        """,
    )
    parser.add_argument(
        "dataset",
        choices=["trainset", "valset", "testset"],
        help="Dataset to evaluate",
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Config YAML (default: acorn_configs/track_building_stage_(3)/gnn_segmentmatcher_walkthrough.yaml)",
    )
    parser.add_argument("--score-cut",       type=float, default=None)
    parser.add_argument("--use-gt-segments", action="store_true", default=False)
    parser.add_argument("--no-wrangler",     action="store_true", default=False)

    args = parser.parse_args()

    config_path = (
        Path(args.config)
        if args.config
        else PIPELINE_ROOT / "acorn_configs" / "track_building_stage_(3)" / "segment_evaluation.yaml"
    )
    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")

    with open(config_path) as f:
        config = yaml.safe_load(f)

    run_segment_evaluation(args.dataset, config, args)


if __name__ == "__main__":
    main()
