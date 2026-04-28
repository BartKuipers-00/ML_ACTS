#!/usr/bin/env python3
"""
Full GNN Chain — Step 2: Reconstruction & Evaluation

Run this after simulate.py, in the acorn conda environment:

    conda activate acorn
    cd /data/alice/bkuipers/low_pt_gnn_pipeline/full_chain_GNN
    python full_chain_gnn.py [--config ../acorn_configs/full_chain_gnn.yaml]

Runs stages 1–5 of the GNN tracking pipeline on CSV data produced by simulate.py:
  1. CSV → PyG          →  feature store graphs
  2. Latent graphs      →  KNN edges via metric learning embedding
  3. GNN inference      →  edge scores (InteractionGNN)
  4. Track building     →  hit labels (GNN segment matching)
  5. Evaluation         →  efficiency / fake rate / clone rate + plots

Analogous to acts/Examples/Scripts/Python/perfect_spacepoints_multigen.py
but using the GNN tracking chain instead of the Kalman Filter.
All pretrained model paths are specified in the YAML config.
"""

import argparse
import sys
from pathlib import Path
import yaml

SCRIPT_DIR    = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent

sys.path.insert(0, str(PIPELINE_ROOT.parent / "acorn"))
sys.path.insert(0, str(PIPELINE_ROOT))

from low_pt_custom_utils.pipeline_chain_utils import (
    stage_csv_to_pyg,
    stage_build_graphs,
    stage_gnn_inference,
    stage_build_tracks,
    stage_evaluate,
)


def _resolve(path_str):
    p = Path(path_str)
    return p if p.is_absolute() else PIPELINE_ROOT / p


def _banner(title):
    print(f"\n{'─' * 70}")
    print(f"  {title}")
    print(f"{'─' * 70}\n")


def main():
    parser = argparse.ArgumentParser(
        description="GNN reconstruction & evaluation step of the full chain",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Run simulate.py first (ACTS env), then this script (acorn env).",
    )
    parser.add_argument(
        "--config",
        default=str(PIPELINE_ROOT / "acorn_configs" / "full_chain_gnn.yaml"),
        help="Path to full_chain_gnn.yaml (default: ../acorn_configs/full_chain_gnn.yaml)",
    )
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = PIPELINE_ROOT / config_path

    with open(config_path) as f:
        config = yaml.safe_load(f)

    run_dir       = _resolve(config.get("run_dir", "data/full_chain_run"))
    dataset       = config.get("dataset", "testset")
    num_events    = config["simulation"]["num_events"]
    debug_files   = config.get("debug_files", False)
    models        = config["models"]
    latent_ckpt   = _resolve(models["latent_checkpoint"])
    gnn_ckpt      = _resolve(models["gnn_checkpoint"])
    mini_gnn_ckpt = _resolve(models["mini_gnn_checkpoint"])

    csv_dir = run_dir / "csv"
    if not csv_dir.exists():
        print(f"ERROR: CSV data not found at {csv_dir}")
        print("Run simulate.py first (in the ACTS/LCG environment).")
        raise SystemExit(1)

    print("=" * 70)
    print("FULL GNN CHAIN — STEP 2: RECONSTRUCTION & EVALUATION")
    print("=" * 70)
    print(f"  Run dir:     {run_dir}")
    print(f"  Dataset:     {dataset}  ({num_events} events)")
    print(f"  Latent:      {latent_ckpt.name}")
    print(f"  GNN:         {gnn_ckpt.name}")
    print(f"  Mini-GNN:    {mini_gnn_ckpt.name}")
    print(f"  Debug files: {debug_files}")
    print("=" * 70)

    _banner("STAGE 1: CSV → PyG")
    stage_csv_to_pyg(run_dir, num_events, dataset)

    _banner("STAGE 2: LATENT GRAPH CONSTRUCTION")
    stage_build_graphs(run_dir, latent_ckpt, config["graph_construction"], dataset, debug_files=debug_files)

    _banner("STAGE 3: GNN EDGE CLASSIFICATION")
    stage_gnn_inference(run_dir, gnn_ckpt, num_events, dataset, config.get("gnn_inference"), debug_files=debug_files)

    _banner("STAGE 4: GNN SEGMENT MATCHING")
    stage_build_tracks(run_dir, mini_gnn_ckpt, num_events, dataset, config["track_building"], debug_files=debug_files)

    _banner("STAGE 5: EVALUATION")
    summary = stage_evaluate(run_dir, dataset, config["evaluation"], debug_files=debug_files)

    print("\n" + "=" * 70)
    print("FULL CHAIN COMPLETE")
    print("=" * 70)
    eff   = summary.get("efficiency",  "N/A")
    fake  = summary.get("fake_rate",   "N/A")
    clone = summary.get("clone_rate",  "N/A")
    if isinstance(eff,   float): eff   = f"{eff:.1%}"
    if isinstance(fake,  float): fake  = f"{fake:.1%}"
    if isinstance(clone, float): clone = f"{clone:.1%}"
    print(f"  Efficiency:  {eff}")
    print(f"  Fake rate:   {fake}")
    print(f"  Clone rate:  {clone}")
    print(f"  Results:     {run_dir / 'track_evaluation' / dataset}/")
    print(f"  Plots:       {run_dir / 'visuals' / 'track_metrics' / dataset}/")
    print("=" * 70)


if __name__ == "__main__":
    main()
