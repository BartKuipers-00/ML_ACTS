#!/usr/bin/env python3
"""
Run GNN edge-score inference for one multiplicity-sweep dataset.

Usage:
    python run_sweep_infer_gnn.py --data-dir /path/to/data_30 \
        --ood-config /path/to/ood_experiment.yaml
"""

import sys
import argparse
import yaml
import tempfile
import os
from pathlib import Path

EXPERIMENT_DIR = Path(__file__).resolve().parent
WORKSPACE_ROOT = EXPERIMENT_DIR.parent
PIPELINE_ROOT = WORKSPACE_ROOT / "low_pt_gnn_pipeline"

sys.path.insert(0, str(WORKSPACE_ROOT / "acorn"))
sys.path.insert(0, str(PIPELINE_ROOT))

from acorn.core.infer_stage import infer as acorn_infer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=str, required=True)
    parser.add_argument("--ood-config", type=str,
                        default=str(EXPERIMENT_DIR / "ood_experiment.yaml"))
    args = parser.parse_args()

    data_dir = Path(args.data_dir).resolve()

    with open(args.ood_config) as f:
        ood = yaml.safe_load(f)

    gnn = ood["gnn_inference"]
    n_events = ood["n_events"]

    # Load base pipeline config then override with ood params
    base_config_path = PIPELINE_ROOT / "acorn_configs" / "gnn_stage_(2)" / "gnn_infer.yaml"
    with open(base_config_path) as f:
        config = yaml.safe_load(f)

    config["input_dir"] = str(data_dir / "graph_constructed_latent")
    config["stage_dir"] = str(data_dir / "gnn_stage")
    config["data_split"] = [0, 0, n_events]
    config["checkpoint"] = gnn["checkpoint"]
    config["accelerator"] = gnn["accelerator"]
    config["node_features"] = gnn["node_features"]
    config["node_scales"] = gnn["node_scales"]

    checkpoint_path = config["checkpoint"]
    if not Path(checkpoint_path).is_absolute():
        checkpoint_path = str(PIPELINE_ROOT / checkpoint_path)
    if not Path(checkpoint_path).exists():
        raise FileNotFoundError(f"GNN checkpoint not found: {checkpoint_path}")

    print("="*70)
    print("GNN Inference (OOD sweep)")
    print("="*70)
    print(f"ood_config : {args.ood_config}")
    print(f"input_dir  : {config['input_dir']}")
    print(f"stage_dir  : {config['stage_dir']}")
    print(f"checkpoint : {Path(checkpoint_path).name}")
    print(f"n_events   : {n_events}")
    print()

    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False, prefix="sweep_gnn_infer_")
    try:
        yaml.dump(config, tmp)
        tmp.close()
        acorn_infer(tmp.name, checkpoint=checkpoint_path)
    finally:
        os.unlink(tmp.name)

    print()
    print("="*70)
    print("GNN inference complete.")
    print("="*70)


if __name__ == "__main__":
    main()
