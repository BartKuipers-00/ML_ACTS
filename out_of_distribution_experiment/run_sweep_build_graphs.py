#!/usr/bin/env python3
"""
Build latent-space KNN graphs for one multiplicity-sweep dataset.

Usage:
    python run_sweep_build_graphs.py --data-dir /path/to/data_30 \
        --ood-config /path/to/ood_experiment.yaml
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

from build_latent_graphs_fast import load_config, load_model, run_graph_construction


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=str, required=True)
    parser.add_argument("--ood-config", type=str,
                        default=str(EXPERIMENT_DIR / "ood_experiment.yaml"))
    args = parser.parse_args()

    data_dir = Path(args.data_dir).resolve()

    with open(args.ood_config) as f:
        ood = yaml.safe_load(f)

    lgc = ood["latent_graph_construction"]
    n_events = ood["n_events"]

    # Load base pipeline config (resolves absolute paths) then override
    config = load_config()
    config["input_dir"] = str(data_dir / "feature_store")
    config["output_dir"] = str(data_dir / "graph_constructed_latent")
    config["datasets"] = {"testset": n_events}
    config["k_max"] = lgc["k_max"]
    config["r_max"] = lgc["r_max"]
    config["r_max_geometric"] = lgc["r_max_geometric"]
    config["segmented"] = lgc["segmented"]
    config["one_in_one_out"] = lgc["one_in_one_out"]
    config["dr_same_layer_cut"] = lgc["dr_same_layer_cut"]
    config["device"] = lgc["device"]

    checkpoint_name = lgc["checkpoint"]
    if not checkpoint_name.endswith(".ckpt"):
        checkpoint_name += ".ckpt"
    checkpoint_path = PIPELINE_ROOT / "saved_models" / checkpoint_name
    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Latent model checkpoint not found: {checkpoint_path}")

    print("="*70)
    print("Latent Graph Construction (OOD sweep)")
    print("="*70)
    print(f"ood_config : {args.ood_config}")
    print(f"input_dir  : {config['input_dir']}")
    print(f"output_dir : {config['output_dir']}")
    print(f"checkpoint : {checkpoint_path.name}")
    print(f"n_events   : {n_events}")
    print()

    model, hparams = load_model(checkpoint_path)
    run_graph_construction(model, hparams, config)

    print()
    print("="*70)
    print("Graph construction complete.")
    print("="*70)


if __name__ == "__main__":
    main()
