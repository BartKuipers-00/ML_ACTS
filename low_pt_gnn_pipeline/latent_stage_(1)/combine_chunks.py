#!/usr/bin/env python3
"""
Combine chunked feature_store data and split into train/val/test sets.
Preserves original event numbers (no renaming).

Usage:
    python combine_chunks.py --data-split TRAIN VAL TEST

Parses chunk directories to determine available events.
Requires explicit train/val/test split sizes.
"""

import argparse
import shutil
from pathlib import Path
import re
from tqdm import tqdm


def extract_event_number(filepath):
    """Extract event number from filename like 'event000003045.pyg'"""
    match = re.search(r'event(\d+)', filepath.name)
    if match:
        return int(match.group(1))
    return -1


def get_event_files(pyg_file):
    """Get all files for an event: -graph.pyg, -truth.csv, -particles.csv"""
    parent_dir = pyg_file.parent
    # Extract event prefix (e.g., 'event000134996' from 'event000134996-graph.pyg')
    # pyg_file.name is like 'event000134996-graph.pyg'
    event_prefix = re.match(r'(event\d+)', pyg_file.name).group(1)

    files_to_move = []
    # Always move the .pyg file (named event######-graph.pyg)
    files_to_move.append(pyg_file)

    # Move associated CSV files if they exist
    truth_csv = parent_dir / f'{event_prefix}-truth.csv'
    if truth_csv.exists():
        files_to_move.append(truth_csv)

    particles_csv = parent_dir / f'{event_prefix}-particles.csv'
    if particles_csv.exists():
        files_to_move.append(particles_csv)

    return files_to_move


def main():
    parser = argparse.ArgumentParser(description='Combine chunked feature_store data and split into train/val/test sets')
    parser.add_argument('--chunk-dir', type=str, default='../data/feature_store_chunk*',
                        help='Pattern for chunk directories (default: ../data/feature_store_chunk*)')
    parser.add_argument('--output', type=str, default='../data/feature_store',
                        help='Output directory for combined data (default: ../data/feature_store)')
    parser.add_argument('--data-split', type=int, nargs=3, required=True,
                        help='Train/val/test split sizes (e.g., 148000 1000 1000)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print plan without moving files')
    args = parser.parse_args()

    # Find all chunk directories
    chunk_dirs = sorted(Path('.').glob(args.chunk_dir))
    print(f"Found {len(chunk_dirs)} chunk directories")

    # Collect all event files from all chunks
    all_events = []
    for chunk_dir in chunk_dirs:
        trainset_dir = chunk_dir / 'trainset'
        if trainset_dir.exists():
            # Get .pyg files (the actual graph data)
            pyg_files = list(trainset_dir.glob('event*.pyg'))
            all_events.extend(pyg_files)
            print(f"  {chunk_dir.name}: {len(pyg_files)} events")

    print(f"\nTotal events collected: {len(all_events)}")

    if len(all_events) == 0:
        print("ERROR: No events found! Check that chunk directories exist and contain .pyg files.")
        print(f"  Searched pattern: {args.chunk_dir}")
        print(f"  Current directory: {Path.cwd()}")
        return

    # Sort by event number
    all_events.sort(key=extract_event_number)

    # Verify no duplicates
    event_numbers = [extract_event_number(f) for f in all_events]
    if len(event_numbers) != len(set(event_numbers)):
        print("WARNING: Duplicate event numbers detected!")
        duplicates = [n for n in event_numbers if event_numbers.count(n) > 1]
        print(f"Duplicates: {set(duplicates)}")

    print(f"Event range: {event_numbers[0]} to {event_numbers[-1]}")

    # Use user-specified split
    train_size, val_size, test_size = args.data_split
    total_requested = train_size + val_size + test_size

    print(f"\nData split: [{train_size}, {val_size}, {test_size}]")
    print(f"Output: {args.output}")
    print()

    if len(all_events) < total_requested:
        print(f"ERROR: Requested split {train_size}+{val_size}+{test_size}={total_requested} events, "
              f"but only {len(all_events)} events are available.")
        return

    train_events = all_events[:train_size]
    val_events = all_events[train_size:train_size + val_size]
    test_events = all_events[train_size + val_size:train_size + val_size + test_size]

    print(f"\nSplit plan:")
    print(f"  Train: {len(train_events)} events (events {extract_event_number(train_events[0])} - {extract_event_number(train_events[-1])})")
    if val_events:
        print(f"  Val:   {len(val_events)} events (events {extract_event_number(val_events[0])} - {extract_event_number(val_events[-1])})")
    if test_events:
        print(f"  Test:  {len(test_events)} events (events {extract_event_number(test_events[0])} - {extract_event_number(test_events[-1])})")

    if args.dry_run:
        print("\nDry run - no files moved")
        return

    # Create output directories
    output_dir = Path(args.output)
    train_dir = output_dir / 'trainset'
    val_dir = output_dir / 'valset'
    test_dir = output_dir / 'testset'

    train_dir.mkdir(parents=True, exist_ok=True)
    val_dir.mkdir(parents=True, exist_ok=True)
    test_dir.mkdir(parents=True, exist_ok=True)

    # Move files (preserving event numbers!)
    # Includes .pyg + CSV files (truth, particles)
    print("\nMoving files (.pyg + CSV)...")

    for pyg_file in tqdm(train_events, desc="Train"):
        for file_to_move in get_event_files(pyg_file):
            shutil.move(str(file_to_move), train_dir / file_to_move.name)

    for pyg_file in tqdm(val_events, desc="Val"):
        for file_to_move in get_event_files(pyg_file):
            shutil.move(str(file_to_move), val_dir / file_to_move.name)

    for pyg_file in tqdm(test_events, desc="Test"):
        for file_to_move in get_event_files(pyg_file):
            shutil.move(str(file_to_move), test_dir / file_to_move.name)

    print(f"\n✓ Combined data saved to {output_dir}")
    print(f"  Train: {len(list(train_dir.glob('*.pyg')))} files")
    print(f"  Val:   {len(list(val_dir.glob('*.pyg')))} files")
    print(f"  Test:  {len(list(test_dir.glob('*.pyg')))} files")


     # remove empty chunk directories
    for chunk_dir in chunk_dirs:
        shutil.rmtree(chunk_dir)



if __name__ == "__main__":
    main()
