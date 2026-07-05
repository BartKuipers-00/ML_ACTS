#!/usr/bin/env python3
"""
Convert ColliderML (HF: CERN/ColliderML-Release-1) parquet -> ACTS-style per-event CSV.

ColliderML stores one collision event per parquet ROW, with per-hit / per-particle
quantities as list columns. This emits the raw ACTS CSV files that the low-pT
pipeline's reader (ActsCustomLowPTReader, geometry_from_csv path) consumes:

  event{N:09d}-hits.csv                  particle_id, tt, segment_id, simhit_id, tx, ty, tz
  event{N:09d}-measurements.csv          measurement_id, volume_id, layer_id, module_id,
                                         global_x, global_y, global_z
  event{N:09d}-measurement-simhit-map.csv   measurement_id, hit_id
  event{N:09d}-particles_simulated.csv   particle_id, particle_type, vx,vy,vz, px,py,pz, m, q

ColliderML tracker hits ARE the digitized measurements, so simhit_id == hit_id ==
measurement_id == row index (identity map). These are high-pT ODD tracks that do not
loop, so every hit gets segment_id = 1 (one segment per particle) and clean_loops is
skipped entirely.

Run with a pyarrow-capable python (e.g. the LCG_108 view), NOT the acorn env:
    source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh
    python colliderml_to_acts_csv.py
"""

import argparse
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent


def _event_row(tbl, i, cols):
    """Extract row i as {col: np.ndarray} straight from Arrow (no to_pydict blow-up).

    List cells -> per-event numpy via .values.to_numpy(); scalar cells (event_id)
    via .as_py(). Avoids materializing all ~25M list elements as Python objects.
    """
    out = {}
    for c in cols:
        cell = tbl.column(c)[i]
        vals = getattr(cell, "values", None)
        out[c] = vals.to_numpy(zero_copy_only=False) if vals is not None else cell.as_py()
    return out


def _read_events(path):
    """Yield (event_id, {col: np.ndarray}) per parquet row, memory-frugally."""
    tbl = pq.read_table(path)
    cols = [c for c in tbl.column_names if c != "event_id"]
    for i in range(tbl.num_rows):
        yield tbl.column("event_id")[i].as_py(), _event_row(tbl, i, cols)


def convert(hits_pq, particles_pq, out_dir, detector_src, limit=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # particles indexed by event_id for the join
    particles_by_event = {ev: d for ev, d in _read_events(particles_pq)}

    written = 0
    for ev, h in _read_events(hits_pq):
        if limit is not None and written >= limit:
            break
        nhits = len(h["x"])
        idx = np.arange(nhits, dtype=np.int64)  # simhit_id == hit_id == measurement_id

        # --- hits.csv (truth simhits: only the fields the reader maps by simhit_id) ---
        pd.DataFrame({
            "particle_id": h["particle_id"],
            "tt": h["time"],
            "segment_id": np.ones(nhits, dtype=np.int64),   # non-looping -> single segment
            "simhit_id": idx,
            "tx": h["true_x"], "ty": h["true_y"], "tz": h["true_z"],
        }).to_csv(out_dir / f"event{ev:09d}-hits.csv", index=False)

        # --- measurements.csv (global positions taken straight from CSV) ---
        pd.DataFrame({
            "measurement_id": idx,
            "volume_id": h["volume_id"],
            "layer_id": h["layer_id"],
            "module_id": h["surface_id"],
            "global_x": h["x"], "global_y": h["y"], "global_z": h["z"],
        }).to_csv(out_dir / f"event{ev:09d}-measurements.csv", index=False)

        # --- measurement-simhit-map.csv (identity: one measurement per sim hit) ---
        pd.DataFrame({"measurement_id": idx, "hit_id": idx}).to_csv(
            out_dir / f"event{ev:09d}-measurement-simhit-map.csv", index=False)

        # --- particles_simulated.csv ---
        p = particles_by_event.get(ev)
        if p is None:
            raise RuntimeError(f"event {ev}: present in hits but missing from particles parquet")
        pd.DataFrame({
            "particle_id": p["particle_id"],
            "particle_type": p["pdg_id"],
            "vx": p["vx"], "vy": p["vy"], "vz": p["vz"],
            "px": p["px"], "py": p["py"], "pz": p["pz"],
            "m": p["mass"],
            "q": p["charge"],
        }).to_csv(out_dir / f"event{ev:09d}-particles_simulated.csv", index=False)

        written += 1

    # detectors.csv is not used for geometry in the geometry_from_csv path, but the
    # reader's __init__ still loads it -> drop in a valid ODD one.
    shutil.copy(detector_src, out_dir / "detectors.csv")
    print(f"Wrote {written} events to {out_dir}  (+ detectors.csv)")
    return written


def main():
    ap = argparse.ArgumentParser(description="ColliderML parquet -> ACTS CSV")
    ap.add_argument("--sample-dir", default=str(PIPELINE_ROOT / "data_coliderML"),
                    help="dir containing ttbar_pu0_{tracker_hits,particles}_*.parquet")
    ap.add_argument("--prefix", default="ttbar_pu0", help="channel_pileup prefix of the parquet files")
    ap.add_argument("--suffix", default="_100ev", help="filename suffix before .parquet")
    ap.add_argument("--out-dir", default=None, help="output csv dir (default: <sample-dir>/csv)")
    ap.add_argument("--detector-src",
                    default=str(PIPELINE_ROOT / "data_60" / "csv" / "detectors.csv"),
                    help="a valid ODD detectors.csv to copy (loaded but unused for geometry)")
    ap.add_argument("--limit", type=int, default=None, help="max events to convert")
    args = ap.parse_args()

    sd = Path(args.sample_dir)
    hits_pq = sd / f"{args.prefix}_tracker_hits{args.suffix}.parquet"
    particles_pq = sd / f"{args.prefix}_particles{args.suffix}.parquet"
    out_dir = args.out_dir or (sd / "csv")
    for f in (hits_pq, particles_pq):
        if not f.exists():
            raise FileNotFoundError(f)
    convert(hits_pq, particles_pq, out_dir, args.detector_src, limit=args.limit)


if __name__ == "__main__":
    main()
