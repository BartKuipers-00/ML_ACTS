#!/usr/bin/env python3
"""
Combine chunked mini-GNN mining output into per-split directories.

Scans data/track_building/mini_gnn_segments_chunk*/<split>/ for per-event .pyg
files and moves them into data/track_building/mini_gnn_segments/<split>/.

The streaming training script (train_mini_GNN_streaming.py) loads these
individual files on demand via a DataLoader — no consolidation into a single
.pt needed.

Usage:
    python combine_mini_gnn_chunks.py
    python combine_mini_gnn_chunks.py --dry-run
    python combine_mini_gnn_chunks.py --keep-chunks   # don't delete chunk dirs
"""

import argparse
import shutil
import sys
from pathlib import Path

from tqdm import tqdm

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent


def main():
    parser = argparse.ArgumentParser(description="Combine chunked mini-GNN mining output")
    parser.add_argument("--dry-run", action="store_true", help="Print plan without moving or deleting")
    parser.add_argument("--keep-chunks", action="store_true", help="Keep chunk directories after combining")
    args = parser.parse_args()

    track_building_dir = PIPELINE_ROOT / "data" / "track_building"
    chunk_dirs = sorted(track_building_dir.glob("mini_gnn_segments_chunk*"))

    if not chunk_dirs:
        print("ERROR: No chunk directories found in data/track_building/mini_gnn_segments_chunk*/")
        sys.exit(1)

    output_dir = track_building_dir / "mini_gnn_segments"
    print(f"Found {len(chunk_dirs)} chunk directories")
    print(f"Output: {output_dir}")
    print()

    for split_name in ["trainset", "valset", "testset"]:
        all_files = []
        for chunk_dir in chunk_dirs:
            split_dir = chunk_dir / split_name
            if split_dir.exists():
                all_files.extend(split_dir.glob("*.pyg"))

        if not all_files:
            print(f"{split_name}: no files found — skipping")
            continue

        all_files.sort(key=lambda f: f.name)
        print(f"{split_name}: {len(all_files)} event files")

        if args.dry_run:
            print(f"  [dry-run] Would move into {output_dir / split_name}/")
            print()
            continue

        dest_dir = output_dir / split_name
        dest_dir.mkdir(parents=True, exist_ok=True)

        for src in tqdm(all_files, desc=f"  Moving {split_name}", unit="file"):
            dest = dest_dir / src.name
            if dest.exists():
                # Prefer the file with the lower chunk number (already sorted)
                continue
            shutil.move(str(src), dest)

        print(f"  Done → {dest_dir}/")
        print()

    if not args.dry_run and not args.keep_chunks:
        print("Removing chunk directories...")
        for d in chunk_dirs:
            shutil.rmtree(d)
            print(f"  Removed {d.name}")

    print("Done!")


if __name__ == "__main__":
    main()
