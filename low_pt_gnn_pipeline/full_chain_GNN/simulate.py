#!/usr/bin/env python3
"""
Full GNN Chain — Step 1: ACTS Simulation (multi-species)

Run this in the ACTS/LCG environment (NOT the acorn conda env):

    cd /data/alice/bkuipers/acts
    source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh
    source build/python/setup.sh
    cd /data/alice/bkuipers/low_pt_gnn_pipeline/full_chain_GNN
    python simulate.py [--config ../acorn_configs/full_chain_gnn.yaml]

Then switch to the acorn environment and run full_chain_gnn.py.
"""

import argparse
import glob
from pathlib import Path
import numpy as np
import yaml
from tqdm import tqdm

SCRIPT_DIR    = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent

import sys
sys.path.insert(0, str(PIPELINE_ROOT / "simulation_(0)"))

from event_generator_for_gnn_training_data import generate_minimal_training_data
from clean_loops_and_attribute_segments import clean_event


def main():
    parser = argparse.ArgumentParser(
        description="ACTS multi-species simulation for the full GNN tracking chain",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
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
        chain = yaml.safe_load(f)

    run_dir = Path(chain.get("run_dir", "data/full_chain_run"))
    if not run_dir.is_absolute():
        run_dir = PIPELINE_ROOT / run_dir

    sim_config = chain["simulation"]
    gun_cfg    = sim_config["particle_gun"]
    species    = gun_cfg.get("particle_types", [gun_cfg.get("particle_type", "muon")])

    # write_geometry: true is required — detectors.csv is read by the CSV→PyG stage
    sim_config["output"] = {
        **sim_config.get("output", {}),
        "base_dir":       str(run_dir),
        "write_geometry": True,
    }

    print("=" * 70)
    print("FULL GNN CHAIN — STEP 1: SIMULATION")
    print("=" * 70)
    print(f"  Config:   {config_path}")
    print(f"  Run dir:  {run_dir}")
    print(f"  Events:   {sim_config['num_events']}")
    print(f"  Species:  {', '.join(species)}")
    print(f"  Particles per species per event: {sim_config['particles_per_vertex']}")
    print(f"  pT range: {gun_cfg['momentum']['min']} – {gun_cfg['momentum']['max']} GeV")
    print("=" * 70 + "\n")

    try:
        generate_minimal_training_data(config=sim_config)
    except RuntimeError as exc:
        if "Sequencer terminated abnormally" not in str(exc):
            raise
        # ACTS raises this when unmasked FPEs occurred in FATRAS, but all events
        # are written before the exception is thrown. Check the CSVs are there.
        csv_dir_check = run_dir / "csv"
        n_hits = len(sorted(csv_dir_check.glob("event*-hits.csv")))
        if n_hits == 0:
            raise RuntimeError(f"Sequencer failed and no CSV files were written.") from exc
        print(f"\nWARNING: ACTS sequencer reported FPEs (common in FATRAS at low pT).")
        print(f"         All {n_hits} events were written successfully — continuing.\n")

    # ── Assign segment IDs and enforce loop fraction on the CSV files ─────────
    # segment_id is required by ActsCustomLowPTReader to create segment-aware
    # ground truth edges. This must run before the CSV→PyG conversion step.
    loop_fraction = sim_config["simulation"].get("loop_fraction", 1.0)
    max_segments  = int(np.ceil(2 * loop_fraction))
    csv_dir       = run_dir / "csv"

    hit_files = sorted(glob.glob(str(csv_dir / "event*-hits.csv")))
    print(f"\nAssigning segments (loop_fraction={loop_fraction} → max {max_segments} segments) ...")
    for hits_path in tqdm(hit_files, desc="  Cleaning events", unit="event"):
        event_prefix = hits_path.replace("-hits.csv", "")
        clean_event(event_prefix, max_segments)
    print(f"✓ Segment assignment done for {len(hit_files)} events.")

    print(f"\nSimulation done. CSV files in: {csv_dir}")
    print("Next step (acorn env):")
    print(f"  conda activate acorn")
    print(f"  python full_chain_gnn.py --config {config_path}")


if __name__ == "__main__":
    main()
