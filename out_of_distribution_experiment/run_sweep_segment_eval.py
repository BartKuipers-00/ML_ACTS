#!/usr/bin/env python3
"""
OOD Segment Evaluation — single dataset worker

uses ood_experiment.yaml

Reads ood_experiment.yaml, takes one data directory (containing gnn_stage/testset/),
runs CC + optional walkthrough, then helix or mini-GNN segment matching, and evaluates.

condor_submit segment_eval.sub

Usage:
    python run_sweep_segment_eval.py --data-dir data_50
    python run_sweep_segment_eval.py --data-dir /abs/path/data_50 --ood-config ood_experiment.yaml
    python run_sweep_segment_eval.py --data-dir data_50 --skip-build
"""

import argparse
import sys
from pathlib import Path

import torch
import yaml

EXPERIMENT_DIR = Path(__file__).resolve().parent
WORKSPACE_ROOT = EXPERIMENT_DIR.parent
PIPELINE_ROOT = WORKSPACE_ROOT / "low_pt_gnn_pipeline"
TRACK_BUILD_DIR = PIPELINE_ROOT / "track_building_stage_(3)"

sys.path.insert(0, str(WORKSPACE_ROOT / "acorn"))
sys.path.insert(0, str(PIPELINE_ROOT))
sys.path.insert(0, str(TRACK_BUILD_DIR))

from low_pt_custom_utils.track_evaluation_utils import (
    run_evaluation,
    save_evaluation_results,
    run_plotting,
)


def _build_helix_config(ood: dict, data_dir: Path, track_build_dir: Path) -> dict:
    base = PIPELINE_ROOT / "acorn_configs" / "track_building_stage_(3)" / "helix_segmentmatcher_walkthrough.yaml"
    with open(base) as f:
        config = yaml.safe_load(f)
    config["data_split"] = [0, 0, ood["n_events"]]
    config["input_dir"] = str(data_dir / "gnn_stage")
    config["stage_dir"] = str(track_build_dir)
    config.update(ood.get("track_building", {}))
    config.update(ood.get("evaluation", {}))
    return config


def _build_gnn_config(ood: dict, data_dir: Path) -> dict:
    base = (
        PIPELINE_ROOT
        / "acorn_configs"
        / "track_building_stage_(3)"
        / "gnn_segmentmatcher_walkthrough.yaml"
    )
    with open(base) as f:
        config = yaml.safe_load(f)
    config["data_split"] = [0, 0, ood["n_events"]]
    config["input_dir"] = str(data_dir / "gnn_stage")
    config["stage_dir"] = str(data_dir / "track_building")
    config.update(ood.get("evaluation", {}))
    sgnn = ood.get("segment_gnn", {})
    if sgnn:
        config.setdefault("gnn", {}).update(sgnn)
    tb = ood.get("track_building", {})
    for key in ("score_cut", "use_gt_segments", "use_wrangler"):
        if key in tb:
            config[key] = tb[key]
    return config


def main():
    parser = argparse.ArgumentParser(description="OOD segment eval for one data directory")
    parser.add_argument(
        "--data-dir",
        type=str,
        required=True,
        help="Path to data directory (e.g. data_50 or /abs/path/data_50)",
    )
    parser.add_argument(
        "--ood-config",
        type=str,
        default=str(EXPERIMENT_DIR / "ood_experiment.yaml"),
        help="Path to ood_experiment.yaml",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Skip track building; only re-run evaluation on existing output",
    )
    parser.add_argument(
        "--device",
        type=str,
        default=None,
        help="Torch device for mini-GNN (default: cuda if available, else cpu)",
    )
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    if not data_dir.is_absolute():
        data_dir = EXPERIMENT_DIR / data_dir
    data_dir = data_dir.resolve()

    with open(args.ood_config) as f:
        ood = yaml.safe_load(f)

    use_segment_gnn = ood.get("use_segment_gnn", False)
    dataset_name = "testset"
    matcher_suffix = "segmentgnn" if use_segment_gnn else "helix"

    gnn_stage_dir = data_dir / "gnn_stage" / dataset_name
    if not gnn_stage_dir.exists():
        print(f"ERROR: gnn_stage/{dataset_name} not found in {data_dir}")
        sys.exit(1)

    matcher_name = "mini-GNN segment embedder" if use_segment_gnn else "helix segment matching"
    print("=" * 70)
    print(f"OOD Segment Eval: {data_dir.name}")
    print("=" * 70)
    print(f"  Config:     {args.ood_config}")
    print(f"  Matcher:    {matcher_name}")
    print(f"  Data dir:   {data_dir}")
    print(f"  Skip build: {args.skip_build}")
    print()

    track_build_dir = data_dir / f"track_building_{matcher_suffix}"
    eval_output_dir = data_dir / f"track_evaluation_{matcher_suffix}" / dataset_name
    plot_output_dir = data_dir / f"visuals_{matcher_suffix}" / "track_metrics" / dataset_name

    if use_segment_gnn:
        from GNN_segmentmatcher_walkthrough import run_gnn_segment_matching
        config = _build_gnn_config(ood, data_dir)
        config["stage_dir"] = str(track_build_dir)

        if not args.skip_build:
            device = args.device or ("cuda" if torch.cuda.is_available() else "cpu")
            from low_pt_custom_utils.mini_gnn_segment_embedding import load_segment_gnn
            model_path = ood["segment_gnn"]["model_path"]
            if not Path(model_path).is_absolute():
                model_path = str(PIPELINE_ROOT / model_path)
            if not Path(model_path).exists():
                raise FileNotFoundError(f"Mini-GNN checkpoint not found: {model_path}")
            print(f"Loading mini-GNN from: {model_path}")
            model = load_segment_gnn(model_path, device=device)
            print()
            run_gnn_segment_matching(dataset_name, config, model, device)
    else:
        from helix_segmentmatcher_walkthrough import run_segment_matching
        config = _build_helix_config(ood, data_dir, track_build_dir)

        if not args.skip_build:
            run_segment_matching(dataset_name, config)

    # Evaluation reads from track_building output
    config["input_dir"] = str(track_build_dir)
    evaluated_events, summary, summary_text = run_evaluation(dataset_name, config)
    save_evaluation_results(evaluated_events, summary, summary_text, dataset_name, eval_output_dir)
    run_plotting(evaluated_events, summary, dataset_name, plot_output_dir, config.get("plots", {}))

    print()
    print(summary_text)
    print("=" * 70)
    print(f"Done: {data_dir.name}")
    print(f"  Evaluation: {eval_output_dir}")
    print("=" * 70)


if __name__ == "__main__":
    main()
