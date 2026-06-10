#!/usr/bin/env python3
"""
Generate grid search configurations and HTCondor batch submission for Latent Stage

Usage:
    cd /data/alice/bkuipers/low_pt_gnn_pipeline/latent_stage_(1)
    python generate_grid_search.py

This creates:
- Individual config files for each grid point
- HTCondor submission file to run all configs in parallel
"""

import yaml
import itertools
from pathlib import Path
import copy


def expand_grid(grid_params, zipped_params=None):
    """Convert grid parameters dict into list of all combinations
    
    Args:
        grid_params: Dict of parameters with their possible values
        zipped_params: List of parameter names that should be zipped together (not Cartesian product)
    """
    # Separate zipped and regular parameters
    if zipped_params:
        # Extract zipped parameters
        zipped_values = {k: grid_params.pop(k) for k in zipped_params if k in grid_params}
        
        # Verify all zipped params have same length
        lengths = [len(v) if isinstance(v, list) else 1 for v in zipped_values.values()]
        if len(set(lengths)) > 1:
            raise ValueError(f"Zipped parameters must have same length: {zipped_values}")
        
        # Create tuples for zipped params
        if zipped_values:
            zipped_combos = list(zip(*[v if isinstance(v, list) else [v] for v in zipped_values.values()]))
        else:
            zipped_combos = [()]
    else:
        zipped_values = {}
        zipped_combos = [()]
    
    # Get regular parameter names and their possible values
    param_names = list(grid_params.keys())
    param_values = [grid_params[k] if isinstance(grid_params[k], list) else [grid_params[k]]
                    for k in param_names]

    # Generate all combinations
    regular_combos = list(itertools.product(*param_values)) if param_names else [()]

    # Combine regular and zipped
    configs = []
    for regular in regular_combos:
        for zipped in zipped_combos:
            config = {name: value for name, value in zip(param_names, regular)}
            config.update({name: value for name, value in zip(zipped_values.keys(), zipped)})
            configs.append(config)

    return configs


def create_run_name(config):
    """Generate a descriptive run name from config"""
    # Key parameters for naming
    h = config['emb_hidden']
    nl = config['nb_layer']
    dim = config['emb_dim']
    lr = config['lr']
    margin = config['margin']
    r = config['r_train']
    act = config['activation']
    ppb = config.get('points_per_batch', 5000)
    knn = config.get('knn', 50)

    return f"h{h}_nl{nl}_dim{dim}_lr{lr}_m{margin}_r{r}_ppb{ppb}_knn{knn}_act-{act}"


def main():
    script_dir = Path(__file__).parent
    pipeline_root = script_dir.parent
    grid_file = pipeline_root / 'acorn_configs' / 'latent_stage_(1)' / 'latent_grid.yml'
    output_dir = pipeline_root / 'acorn_configs' / 'latent_stage_(1)' / 'grid_configs'
    output_dir.mkdir(exist_ok=True)

    print("="*70)
    print("LATENT SPACE (METRIC LEARNING) GRID SEARCH GENERATOR")
    print("="*70)
    print()

    # Load grid specification
    with open(grid_file) as f:
        grid_spec = yaml.safe_load(f)

    grid_params = grid_spec['grid_params']
    fixed_params = grid_spec['fixed_params']
    zipped_params = grid_spec.get('zipped_params', None)

    # Generate all combinations
    print("Generating grid combinations...")
    grid_configs = expand_grid(grid_params, zipped_params)

    print(f"Grid search space: {len(grid_configs)} configurations")
    print()

    # Create individual config files
    for idx, grid_config in enumerate(grid_configs):
        # Merge fixed and grid parameters
        full_config = copy.deepcopy(fixed_params)
        full_config.update(grid_config)

        # Generate run name
        run_name = create_run_name(grid_config)
        full_config['run_name'] = run_name

        # Save config
        config_file = output_dir / f"config_{idx:04d}_{run_name}.yaml"
        with open(config_file, 'w') as f:
            yaml.dump(full_config, f, default_flow_style=False, sort_keys=False)

        if idx < 3:  # Show first few
            print(f"Config {idx}: {run_name}")

    print(f"... ({len(grid_configs) - 3} more)")
    print()
    print(f"Configs saved to: {output_dir}/")
    print()

    # Create HTCondor submission file
    jobs_dir = Path('/data/alice/bkuipers/jobs')
    jobs_dir.mkdir(exist_ok=True)
    (jobs_dir / 'output').mkdir(exist_ok=True)  # Ensure output directory exists
    sub_file = jobs_dir / 'train_latent_gridsearch.sub'

    submission = f"""#!/usr/bin/env condor_submit
# Latent Space (Metric Learning) Grid Search - {len(grid_configs)} parallel CPU jobs

executable              = run_latent_gridsearch.sh
log                     = /data/alice/bkuipers/jobs/output/latent_grid_$(ClusterId).log
output                  = /data/alice/bkuipers/jobs/output/latent_grid_$(ClusterId).$(Process).out
error                   = /data/alice/bkuipers/jobs/output/latent_grid_$(ClusterId).$(Process).err

# Pass config index
arguments               = $(Process)

# Use Alma Linux 9 environment
+UseOS                  = "el9"

# Short jobs - training (up to 4 hours)
+JobCategory            = "short"

# CPU resource requirements
request_cpus            = 4
request_memory          = 32G

# File transfer
should_transfer_files   = NO
initialdir              = /data/alice/bkuipers/low_pt_gnn_pipeline

# Environment
environment             = "CONDA_ENV=acorn"

# Queue one job per config
queue {len(grid_configs)}
"""

    with open(sub_file, 'w') as f:
        f.write(submission)

    print(f"HTCondor submission created: {sub_file}")
    print()

    # Create runner script
    runner_file = jobs_dir / 'run_latent_gridsearch.sh'

    runner = f"""#!/bin/bash
# Run Latent Space (Metric Learning) training with specific grid search config

set -e

CONFIG_IDX=${{1:-0}}

echo "=========================================="
echo "Latent Space Grid Search - Config $CONFIG_IDX"
echo "=========================================="
echo "Started: $(date)"
echo "Node: $(hostname -f)"
echo "CPUs: $(nproc)"
echo ""

# Activate conda
source /data/alice/bkuipers/miniconda3/etc/profile.d/conda.sh
conda activate acorn

# Find config file
CONFIG_FILE=$(ls acorn_configs/latent_stage_\(1\)/grid_configs/config_$(printf "%04d" $CONFIG_IDX)_*.yaml)

echo "Config: $CONFIG_FILE"
echo ""

# Run training
cd /data/alice/bkuipers/low_pt_gnn_pipeline/latent_stage_\(1\)
python train_latent_cluster_learning.py --config "$CONFIG_FILE"

echo ""
echo "=========================================="
echo "Training complete: $(date)"
echo "=========================================="
"""

    with open(runner_file, 'w') as f:
        f.write(runner)

    runner_file.chmod(0o755)

    print(f"Runner script created: {runner_file}")
    print()

    print("="*70)
    print("SETUP COMPLETE!")
    print("="*70)
    print()
    print("To launch grid search:")
    print(f"  cd {jobs_dir}")
    print(f"  condor_submit train_latent_gridsearch.sub")
    print()
    print(f"This will launch {len(grid_configs)} CPU jobs in parallel!")
    print("All results will be logged to W&B under project: Low_pt_latent_map_gridsearch")
    print()


if __name__ == "__main__":
    main()
