#!/usr/bin/env python3
"""
Convert ACTS CSV files to PyG format using ActsReader

This script takes the raw ACTS simulation CSV files and converts them to
PyG graph format with train/val/test splits. The resulting graphs contain:
- Hit coordinates (r, φ, z, (t)) 
- Particle IDs
- Ground truth (sequential) edges 

Usage:
    python convert_csv_to_pyg_sets.py
    
Configuration:
    acorn_configs/convert_csv_to_pyg_sets.yaml
"""

import sys
import argparse
from pathlib import Path
import yaml

# Add acorn to path
SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / 'acorn'))
sys.path.insert(0, str(PIPELINE_ROOT))

# Use custom reader for low-pT data with time-based trajectory ordering
from acts_custom_low_pt_reader import ActsCustomLowPTReader as ActsReader


def main():
    """Convert ACTS CSV files to PyG graphs using ActsReader"""
    
    print("="*80)
    print("CONVERT ACTS CSV TO PYG FORMAT")
    print("="*80)
    print()
    
    # Load configuration (default: muon config; override with --config for other samples)
    ap = argparse.ArgumentParser(description="Convert ACTS CSV files to PyG graphs")
    ap.add_argument('--config', default=None,
                    help="Path to convert config YAML (default: convert_csv_to_pyg_sets.yaml)")
    cli_args, _ = ap.parse_known_args()
    config_path = (Path(cli_args.config).resolve() if cli_args.config
                   else PIPELINE_ROOT / 'acorn_configs' / 'latent_stage_(1)' / 'convert_csv_to_pyg_sets.yaml')
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Make paths absolute relative to pipeline root directory
    for path_key in ['input_dir', 'stage_dir', 'detector_path']:
        if path_key in config and not Path(config[path_key]).is_absolute():
            config[path_key] = str(PIPELINE_ROOT / config[path_key])

    print(f"Configuration: {config_path.name}")
    print(f"Input:  {config['input_dir']}")
    print(f"Output: {config['stage_dir']}")
    print(f"Data split: {config['data_split']}")
    print()
    
    # Create reader and convert
    print("Starting conversion...")
    print("Using ActsCustomLowPTReader for time-based trajectory ordering")
    reader = ActsReader.infer(config)  # ActsReader is now ActsCustomLowPTReader
    
    print()
    print("="*80)
    print("✓ CSV TO PYG CONVERSION COMPLETE!")
    print("="*80)
    print(f"\nPyG graphs saved to: {config['stage_dir']}")
    print(f"  - trainset: {config['data_split'][0]} events")
    print(f"  - valset:   {config['data_split'][1]} events")
    print(f"  - testset:  {config['data_split'][2]} events")
    print()
    print("Next step:")
    print("  python train_latent_cluster_learning.py")
    print()


if __name__ == "__main__":
    main()
