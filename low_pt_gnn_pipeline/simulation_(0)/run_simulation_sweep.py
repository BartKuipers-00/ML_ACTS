#!/usr/bin/env python3
"""
Run ACTS simulation for one multiplicity-sweep dataset.

Loads the base config, overrides particles_per_vertex and output directory,
then runs the simulation. Called by jobs/run_multiplicity_sweep_job.sh.

Usage:
    python run_simulation_sweep.py --particles-per-vertex 30
"""

import sys
import argparse
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent

sys.path.insert(0, str(SCRIPT_DIR))
from event_generator_for_gnn_training_data import load_simulation_config, generate_minimal_training_data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--particles-per-vertex", type=int, required=True,
                        help="Number of particles per vertex (per species)")
    args = parser.parse_args()

    ppv = args.particles_per_vertex
    config = load_simulation_config()

    config["particles_per_vertex"] = ppv
    config["num_events"] = 500
    config["output"]["base_dir"] = str(PIPELINE_ROOT / f"data_{ppv}")

    print(f"particles_per_vertex = {ppv}")
    print(f"num_events           = {config['num_events']}")
    print(f"output               = {config['output']['base_dir']}")

    generate_minimal_training_data(config=config)


if __name__ == "__main__":
    main()
