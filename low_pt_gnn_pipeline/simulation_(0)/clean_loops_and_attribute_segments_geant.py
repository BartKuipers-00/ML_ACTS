#!/usr/bin/env python3
"""
Clean Geant4 simulation CSV data: assign radial segments and enforce loop fraction.

Same logic as clean_loops_and_attribute_segments.py but adapted for Geant4 output,
which encodes particle identity across five columns:
  particle_id_pv, particle_id_sv, particle_id_part, particle_id_gen, particle_id_subpart

These are combined into a single 64-bit integer using the ACTS barcode encoding:
  particle_id = (pv << 44) | (sv << 32) | (part << 22) | (gen << 16) | subpart

The combined column is added to hits.csv so that downstream readers that expect a
single 'particle_id' column will work correctly.

Usage:
    python clean_loops_and_attribute_segments_geant.py --data-dir /path/to/data_geant
    python clean_loops_and_attribute_segments_geant.py --data-dir /path/to/data_geant --loop-fraction 1.0
"""

import sys
import argparse
import glob
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent

PARTICLE_ID_COLS = [
    "particle_id_pv", "particle_id_sv", "particle_id_part",
    "particle_id_gen", "particle_id_subpart",
]


def encode_particle_id(df):
    """
    Combine five ACTS barcode columns into a single 64-bit particle_id integer.

    Encoding: (pv << 44) | (sv << 32) | (part << 22) | (gen << 16) | subpart
    A combined value of 0 (all components zero) represents noise/background.
    """
    pv   = df["particle_id_pv"].to_numpy(dtype=np.int64)
    sv   = df["particle_id_sv"].to_numpy(dtype=np.int64)
    part = df["particle_id_part"].to_numpy(dtype=np.int64)
    gen  = df["particle_id_gen"].to_numpy(dtype=np.int64)
    sub  = df["particle_id_subpart"].to_numpy(dtype=np.int64)
    pid = pd.Series(
        (pv << 44) | (sv << 32) | (part << 22) | (gen << 16) | sub,
        index=df.index,
    )
    return pid


def assign_segments(hits_df):
    """
    Assign radial segment IDs to each hit based on radial momentum sign flips.

    Segment 1 = first outgoing arc, 2 = first incoming arc, 3 = second outgoing, etc.
    Hits with particle_id == 0 (noise) get segment_id = 0.

    Parameters
    ----------
    hits_df : pd.DataFrame
        Hits dataframe. Must have a 'particle_id' column (combined from Geant components),
        plus tx, ty, tpx, tpy, tt.

    Returns
    -------
    pd.DataFrame
        Same dataframe with 'segment_id' column added.
    """
    hits_df["segment_id"] = 0

    signal = hits_df[hits_df["particle_id"] != 0].copy()
    if len(signal) == 0:
        return hits_df

    signal = signal.sort_values(["particle_id", "tt"])

    r = np.sqrt(signal["tx"].values**2 + signal["ty"].values**2)
    r = np.where(r == 0, 1e-10, r)
    p_r = (signal["tpx"].values * signal["tx"].values
           + signal["tpy"].values * signal["ty"].values) / r
    radial_sign = np.sign(p_r)

    segment_ids = np.ones(len(signal), dtype=np.int64)
    particle_ids = signal["particle_id"].values

    particle_boundaries = np.where(np.diff(particle_ids) != 0)[0] + 1
    particle_starts = np.concatenate([[0], particle_boundaries])
    particle_ends = np.concatenate([particle_boundaries, [len(signal)]])

    for start, end in zip(particle_starts, particle_ends):
        current_segment = 1
        last_nonzero_sign = radial_sign[start]

        for k in range(start + 1, end):
            d = radial_sign[k]
            if d != 0 and last_nonzero_sign != 0 and d != last_nonzero_sign:
                current_segment += 1
            if d != 0:
                last_nonzero_sign = d
            segment_ids[k] = current_segment

    hits_df.loc[signal.index, "segment_id"] = segment_ids
    return hits_df


def clean_event(event_prefix, max_segments):
    """
    Clean a single Geant4 event: combine particle ID columns, assign segments,
    remove excess hits, and update related CSV files.

    Modifies in-place:
      - hits.csv:                   Adds 'particle_id' (combined), 'simhit_id',
                                    'segment_id'; removes excess hits.
      - measurement-simhit-map.csv: Removes entries for deleted hits.
      - measurements.csv:           Removes measurements whose simhit was deleted.
      - particles_simulated.csv:    Removes particles that lost all hits.

    Parameters
    ----------
    event_prefix : str
        Path prefix for event files (e.g. '.../event000000000').
    max_segments : int
        Maximum allowed segment ID.

    Returns
    -------
    tuple (int, int)
        (hits_removed, particles_removed)
    """
    hits_path = f"{event_prefix}-hits.csv"
    hits = pd.read_csv(hits_path)

    # Combine particle ID components into a single integer (if not already done)
    if "particle_id" not in hits.columns:
        hits["particle_id"] = encode_particle_id(hits)

    # Record original SimHitContainer positions BEFORE any removal
    if "simhit_id" not in hits.columns:
        hits["simhit_id"] = hits.index

    # Assign segments
    hits = assign_segments(hits)

    # Remove hits beyond allowed segments (keep noise hits with segment_id=0)
    mask_keep = (hits["segment_id"] <= max_segments) | (hits["segment_id"] == 0)
    # Remove hits with |z| > 599 mm
    mask_keep = mask_keep & (hits["tz"].abs() <= 599)
    hits_removed = (~mask_keep).sum()

    removed_simhit_ids = set(hits.loc[~mask_keep, "simhit_id"])

    hits = hits[mask_keep].reset_index(drop=True)
    hits.to_csv(hits_path, index=False)

    # Update measurement-simhit-map and measurements.csv
    map_path = f"{event_prefix}-measurement-simhit-map.csv"
    meas_path = f"{event_prefix}-measurements.csv"
    if Path(map_path).exists() and len(removed_simhit_ids) > 0:
        simhit_map = pd.read_csv(map_path)
        mask_map_keep = ~simhit_map["hit_id"].isin(removed_simhit_ids)
        removed_meas_ids = set(simhit_map.loc[~mask_map_keep, "measurement_id"])
        simhit_map = simhit_map[mask_map_keep]
        simhit_map.to_csv(map_path, index=False)

        if Path(meas_path).exists() and len(removed_meas_ids) > 0:
            meas = pd.read_csv(meas_path)
            meas = meas[~meas["measurement_id"].isin(removed_meas_ids)]
            meas.to_csv(meas_path, index=False)

    # Clean particles_simulated.csv (no particles_initial in Geant output)
    surviving_pids = set(hits["particle_id"].unique()) - {0}
    particles_removed = 0
    for suffix in ["particles_initial", "particles_simulated"]:
        p_path = f"{event_prefix}-{suffix}.csv"
        if not Path(p_path).exists():
            continue
        particles = pd.read_csv(p_path)
        # Build combined particle_id for the particles file
        if "particle_id" not in particles.columns:
            particles["particle_id"] = encode_particle_id(particles)
        n_before = len(particles)
        particles = particles[particles["particle_id"].isin(surviving_pids)]
        if suffix in ("particles_initial", "particles_simulated"):
            particles_removed += n_before - len(particles)
        if len(particles) < n_before:
            particles.to_csv(p_path, index=False)

    return hits_removed, particles_removed


def main():
    parser = argparse.ArgumentParser(
        description="Assign radial segments and enforce loop fraction on Geant4 simulation CSVs"
    )
    parser.add_argument(
        "--loop-fraction", "-f", type=float, default=None,
        help="Loop fraction (e.g. 1.0 = one full loop = 2 segments). "
             "Defaults to value in event_generator_simulation.yaml."
    )
    parser.add_argument(
        "--data-dir", type=str, default=None,
        help="Data directory containing a csv/ subdirectory. "
             "Defaults to base_dir from config."
    )
    parser.add_argument(
        "--chunk", type=int, default=None,
        help="Chunk index (0-indexed) for parallel processing."
    )
    parser.add_argument(
        "--total-chunks", type=int, default=None,
        help="Total number of chunks for parallel processing."
    )
    args = parser.parse_args()

    import yaml
    config_path = PIPELINE_ROOT / "acorn_configs" / "simulation_(0)" / "event_generator_simulation.yaml"
    with open(config_path, 'r') as f_cfg:
        config = yaml.safe_load(f_cfg)

    if args.loop_fraction is None:
        config_loop_fraction = config.get('simulation', {}).get('loop_fraction', None)
        if config_loop_fraction is not None:
            args.loop_fraction = float(config_loop_fraction)
        else:
            print("ERROR: --loop-fraction is required (not found in config either).")
            sys.exit(1)

    max_segments = int(np.ceil(2 * args.loop_fraction))

    if args.data_dir is not None:
        csv_dir = Path(args.data_dir) / "csv"
    else:
        base_dir = config['output']['base_dir']
        if not Path(base_dir).is_absolute():
            base_dir = str(PIPELINE_ROOT / base_dir)
        csv_dir = Path(base_dir) / "csv"

    if not csv_dir.exists():
        print(f"ERROR: CSV directory not found: {csv_dir}")
        sys.exit(1)

    all_hit_files = sorted(glob.glob(str(csv_dir / "event*-hits.csv")))
    if not all_hit_files:
        print(f"No hit files found in {csv_dir}")
        sys.exit(1)

    if args.chunk is not None and args.total_chunks is not None:
        n = len(all_hit_files)
        events_per_chunk = n // args.total_chunks
        remainder = n % args.total_chunks
        if args.chunk < remainder:
            start = args.chunk * (events_per_chunk + 1)
            end = start + events_per_chunk + 1
        else:
            start = remainder * (events_per_chunk + 1) + (args.chunk - remainder) * events_per_chunk
            end = start + events_per_chunk
        hit_files = all_hit_files[start:end]
    else:
        hit_files = all_hit_files

    print("=" * 80)
    print("CLEAN LOOPS AND ASSIGN SEGMENTS (Geant4)")
    print("=" * 80)
    print(f"CSV directory:    {csv_dir}")
    print(f"Events found:     {len(all_hit_files)}")
    if args.chunk is not None:
        print(f"Chunk:            {args.chunk}/{args.total_chunks} ({len(hit_files)} events)")
    print(f"Loop fraction:    {args.loop_fraction}")
    print(f"Max segments:     {max_segments}")
    print(f"  (segment 1=outgoing, 2=incoming, 3=outgoing, ...)")
    print()

    total_hits_removed = 0
    total_particles_removed = 0
    events_affected = 0

    for hits_path in tqdm(hit_files, desc="Cleaning events", unit="event"):
        event_prefix = hits_path.replace("-hits.csv", "")

        hits_removed, particles_removed = clean_event(event_prefix, max_segments)

        if hits_removed > 0 or particles_removed > 0:
            events_affected += 1
        total_hits_removed += hits_removed
        total_particles_removed += particles_removed

    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Events processed:    {len(hit_files)} / {len(all_hit_files)}")
    print(f"Events affected:     {events_affected}")
    print(f"Total hits removed:  {total_hits_removed}")
    print(f"Total particles removed (all hits cut): {total_particles_removed}")
    print()


if __name__ == "__main__":
    main()
