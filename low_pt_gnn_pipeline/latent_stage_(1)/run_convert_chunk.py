#!/usr/bin/env python3
"""
Run CSV to PyG conversion for a specific event range 

Usage:
    python run_convert_chunk.py --chunk 0 --total-chunks 10
"""

import sys
import argparse
import glob
from pathlib import Path
import yaml

# Add acorn to path
SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / 'acorn'))
sys.path.insert(0, str(PIPELINE_ROOT))

from acts_custom_low_pt_reader import ActsCustomLowPTReader
from low_pt_custom_utils.region_label_utils import VALID_DETECTORS


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chunk', type=int, required=True, help='Chunk number (0-indexed)')
    parser.add_argument('--total-chunks', type=int, required=True, help='Total number of chunks')
    parser.add_argument('--detector', type=str, default=None, choices=list(VALID_DETECTORS),
                        help="Detector geometry selecting the region map (generic=FATRAS, "
                             "odd=Geant4/ODD). Prompts if omitted in an interactive shell.")
    parser.add_argument('--config', type=str, default='acorn_configs/latent_stage_(1)/convert_csv_to_pyg_sets.yaml')
    parser.add_argument('--input-dir', type=str, default=None,
                        help='Override input_dir (CSV source). detector_path is set to '
                             '<input-dir>/detectors.csv. Default: value from config.')
    parser.add_argument('--stage-dir', type=str, default=None,
                        help='Override output stage_dir. Default: data/feature_store_chunk{NN}/.')
    args = parser.parse_args()
    
    print("="*70)
    print(f"CSV to PyG Conversion - Chunk {args.chunk}/{args.total_chunks}")
    print("="*70)
    
    # Load base config
    config_path = PIPELINE_ROOT / args.config
    with open(config_path) as f:
        config = yaml.safe_load(f)

    # Detector selects the region map (resolved in ActsCustomLowPTReader)
    config['detector'] = args.detector

    # Optional directory overrides (used by the geant clean+convert flow)
    if args.input_dir is not None:
        config['input_dir'] = args.input_dir
        config['detector_path'] = str(Path(args.input_dir) / "detectors.csv")

    # Calculate event range for this chunk based on actual CSV files present
    input_dir = config['input_dir']
    if not Path(input_dir).is_absolute():
        input_dir = str(PIPELINE_ROOT / input_dir)
    total_events = len(glob.glob(str(Path(input_dir) / "event*-hits.csv")))
    if total_events == 0:
        print(f"ERROR: No hit files found in {input_dir}")
        sys.exit(1)
    print(f"Total events (from CSV files): {total_events}")
    events_per_chunk = total_events // args.total_chunks
    remainder = total_events % args.total_chunks
    
    # Distribute remainder across first chunks
    if args.chunk < remainder:
        start_event = args.chunk * (events_per_chunk + 1)
        end_event = start_event + events_per_chunk + 1
    else:
        start_event = remainder * (events_per_chunk + 1) + (args.chunk - remainder) * events_per_chunk
        end_event = start_event + events_per_chunk
    
    chunk_size = end_event - start_event

    stage_dir = args.stage_dir if args.stage_dir is not None \
        else f"data/feature_store_chunk{args.chunk:02d}/"

    print(f"Input:  {input_dir}")
    print(f"Processing events: {start_event} to {end_event-1} ({chunk_size} events)")
    print(f"Output: {stage_dir}")
    print()

    # Modify config for this chunk
    config['data_split'] = [chunk_size, 0, 0]  # All in trainset
    config['input_sets'] = ['trainset']
    config['stage_dir'] = stage_dir
    
    # Override the dataset to only include our event range
    config['event_range'] = [start_event, end_event]
    
    print(f"Workers: {config.get('max_workers', 1)}")
    print()
    
    # Run conversion
    reader = ActsCustomLowPTReader.infer(config)
    
    print()
    print("="*70)
    print(f"✓ Chunk {args.chunk} complete!")
    print("="*70)


if __name__ == "__main__":
    main()
