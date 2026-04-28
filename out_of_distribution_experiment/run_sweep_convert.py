#!/usr/bin/env python3
"""
Convert one multiplicity-sweep dataset (data_N/) from CSV to PyG graphs.

Only produces a testset (no train/val split needed for evaluation sweeps).
Called by jobs/run_multiplicity_sweep_job.sh after simulation and cleaning.

Usage:
    python run_sweep_convert.py --data-dir /path/to/data_30 --n-events 500
"""

import sys
import argparse
import yaml
from pathlib import Path

EXPERIMENT_DIR = Path(__file__).resolve().parent
WORKSPACE_ROOT = EXPERIMENT_DIR.parent
PIPELINE_ROOT = WORKSPACE_ROOT / "low_pt_gnn_pipeline"
LATENT_DIR = PIPELINE_ROOT / "latent_stage_(1)"

sys.path.insert(0, str(WORKSPACE_ROOT / "acorn"))
sys.path.insert(0, str(PIPELINE_ROOT))
sys.path.insert(0, str(LATENT_DIR))

from acts_custom_low_pt_reader import ActsCustomLowPTReader


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=str, required=True,
                        help="Absolute path to the dataset directory (e.g. /path/to/data_30)")
    parser.add_argument("--n-events", type=int, default=500,
                        help="Number of events (all go into testset)")
    args = parser.parse_args()

    data_dir = Path(args.data_dir).resolve()

    config_path = PIPELINE_ROOT / "acorn_configs" / "latent_stage_(1)" / "convert_csv_to_pyg_sets.yaml"
    with open(config_path) as f:
        config = yaml.safe_load(f)

    config["input_dir"] = str(data_dir / "csv")
    config["stage_dir"] = str(data_dir / "feature_store")
    config["detector_path"] = str(data_dir / "csv" / "detectors.csv")
    config["data_split"] = [0, 0, args.n_events]
    config["input_sets"] = ["testset"]

    print("="*70)
    print("CSV → PyG Conversion (multiplicity sweep)")
    print("="*70)
    print(f"input_dir  : {config['input_dir']}")
    print(f"stage_dir  : {config['stage_dir']}")
    print(f"n_events   : {args.n_events} (testset only)")
    print()

    ActsCustomLowPTReader.infer(config)

    print()
    print("="*70)
    print("Conversion complete.")
    print("="*70)


if __name__ == "__main__":
    main()
