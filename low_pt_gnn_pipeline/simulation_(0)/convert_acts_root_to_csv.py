#!/usr/bin/env python3
"""
Convert ACTS ROOT output (RootSimHitWriter / RootMeasurementWriter / RootParticleWriter)
into the per-event ACORN CSV format consumed by the low-pT pipeline.

Source ROOT files (all events concatenated, one row/entry per hit/measurement):
    hits.root             tree 'hits'          -> event*-hits.csv
    measurements.root     tree 'measurements'  -> event*-measurements.csv
                                               -> event*-measurement-simhit-map.csv
    particles_simulation  tree 'particles'     -> event*-particles_simulated.csv

This was written for the ALICE 3 (IRIS) detector sample, whose custom geometry has no
detectors.csv. We therefore copy geometry (volume/layer/module) and global positions
straight from the ROOT branches into the CSV, and the custom reader reads them via the
`geometry_from_csv` flag (no local->global recompute needed).

particle_id is the 64-bit ACTS barcode, packed from its five components:
    (vertex_primary<<52)|(vertex_secondary<<40)|(particle<<24)|(generation<<16)|sub_particle

Usage:
    python convert_acts_root_to_csv.py \
        --root-dir /data/alice/pbutti/.../<run> \
        --out-dir  ../data_alice3_10/csv \
        --first 98000 --n 10
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import uproot

U = np.uint64


def pack_barcode(vp, vs, pa, ge, su):
    """Pack the five ACTS barcode sub-fields into the 64-bit particle_id."""
    return (
        (vp.astype(U) << U(52))
        | (vs.astype(U) << U(40))
        | (pa.astype(U) << U(24))
        | (ge.astype(U) << U(16))
        | su.astype(U)
    )


def _first(vecs):
    """First contributor of each measurement (jagged vector -> flat array)."""
    return np.array([v[0] if len(v) > 0 else 0 for v in vecs], dtype=np.uint64)


def write_event(parts_df, hits_df, meas_df, ev, out_idx, out_dir):
    """Build and write the four CSV files for one event (already-sliced DataFrames)."""
    nh = len(hits_df)
    hits_df = hits_df.copy()
    hits_df["simhit_id"] = np.arange(nh, dtype=np.int64)  # SimHitContainer row position
    hits_df["segment_id"] = 1  # single time-ordered segment (no loop splitting for this sample)

    # Keep only particles that actually produced a hit (drops ~170k->~few-k secondaries
    # with no hits; the rest are dead weight for training/evaluation).
    hit_pids = set(np.unique(hits_df["particle_id"].to_numpy()))
    parts_df = parts_df[parts_df["particle_id"].isin(hit_pids)].reset_index(drop=True)

    nm = len(meas_df)
    meas_df = meas_df.copy()
    meas_df["measurement_id"] = np.arange(nm, dtype=np.int64)

    # ---- link each measurement to its sim hit by (particle_id, exact truth pos) --
    link = meas_df.merge(
        hits_df[["particle_id", "tx", "ty", "tz", "simhit_id", "geometry_id"]],
        left_on=["_pid", "_tx", "_ty", "_tz"],
        right_on=["particle_id", "tx", "ty", "tz"], how="left",
    )
    link = link.drop_duplicates("measurement_id", keep="first")  # resolve rare ambiguous matches
    matched = link[link["simhit_id"].notna()].copy()
    matched["simhit_id"] = matched["simhit_id"].astype(np.int64)
    matched["geometry_id"] = matched["geometry_id"].astype(np.uint64)  # un-float (NaN-promoted by join)

    meas_out = matched[["measurement_id", "geometry_id", "volume_id", "layer_id",
                        "module_id", "local0", "local1", "global_x", "global_y", "global_z"]]
    simhit_map = matched[["measurement_id", "simhit_id"]].rename(columns={"simhit_id": "hit_id"})

    # ---- write -------------------------------------------------------------
    tag = f"event{out_idx:09d}"  # ACORN expects events numbered from 0
    hits_df.drop(columns=["event_id", "volume_id", "layer_id"]).to_csv(out_dir / f"{tag}-hits.csv", index=False)
    meas_out.to_csv(out_dir / f"{tag}-measurements.csv", index=False)
    simhit_map.to_csv(out_dir / f"{tag}-measurement-simhit-map.csv", index=False)
    parts_df.to_csv(out_dir / f"{tag}-particles_simulated.csv", index=False)
    # cells unused (use_cell_information=false) -- header-only placeholder
    (out_dir / f"{tag}-cells.csv").write_text(
        "geometry_id,measurement_id,channel0,channel1,timestamp,value\n")
    return nh, len(meas_out)


def write_minimal_detectors(hits_t, meas_t, out_dir):
    """A tiny detectors.csv (one row per volume/layer) just to satisfy reader init.

    The custom reader reads geometry+global positions from the CSV columns directly
    (geometry_from_csv), so the translation/rotation values here are never used.
    """
    cols = ["volume_id", "layer_id"]
    hv = hits_t.arrays(cols, library="np")
    pairs = pd.DataFrame({"volume_id": hv["volume_id"], "layer_id": hv["layer_id"]}).drop_duplicates()
    pairs["module_id"] = 0
    pairs = pairs.sort_values(["volume_id", "layer_id"]).reset_index(drop=True)
    pairs["geometry_id"] = np.arange(1, len(pairs) + 1)
    for c, val in [("boundary_id", 0), ("extra_id", 0),
                   ("cx", 0.0), ("cy", 0.0), ("cz", 0.0),
                   ("rot_xu", 1.0), ("rot_xv", 0.0), ("rot_xw", 0.0),
                   ("rot_yu", 0.0), ("rot_yv", 1.0), ("rot_yw", 0.0),
                   ("rot_zu", 0.0), ("rot_zv", 0.0), ("rot_zw", 1.0),
                   ("module_t", -1.0), ("pitch_u", -1.0), ("pitch_v", -1.0)]:
        pairs[c] = val
    order = ["geometry_id", "volume_id", "boundary_id", "layer_id", "module_id", "extra_id",
             "cx", "cy", "cz", "rot_xu", "rot_xv", "rot_xw", "rot_yu", "rot_yv", "rot_yw",
             "rot_zu", "rot_zv", "rot_zw", "module_t", "pitch_u", "pitch_v"]
    pairs[order].to_csv(out_dir / "detectors.csv", index=False)


def entry_bounds(tree, branch, lo, hi):
    """Entry [start, stop) covering event ids [lo, hi).

    Events are written sequentially by ACTS, so the id branch is sorted and we can
    read only the relevant baskets instead of the whole (26M-row) tree.
    """
    ev = tree[branch].array(library="np")
    assert np.all(np.diff(ev) >= 0), f"{branch} not sorted; entry-range read unsafe"
    return int(np.searchsorted(ev, lo, "left")), int(np.searchsorted(ev, hi, "left"))


def read_hits(hits_t, lo, hi):
    s, e = entry_bounds(hits_t, "event_id", lo, hi)
    a = hits_t.arrays(
        ["event_id", "geometry_id", "volume_id", "layer_id",
         "tx", "ty", "tz", "tt", "tpx", "tpy", "tpz", "te",
         "deltapx", "deltapy", "deltapz", "deltae", "index",
         "barcode_vertex_primary", "barcode_vertex_secondary", "barcode_particle",
         "barcode_generation", "barcode_sub_particle"],
        entry_start=s, entry_stop=e, library="np",
    )
    df = pd.DataFrame({
        "event_id": a["event_id"],
        "particle_id": pack_barcode(a["barcode_vertex_primary"], a["barcode_vertex_secondary"],
                                    a["barcode_particle"], a["barcode_generation"], a["barcode_sub_particle"]),
        "geometry_id": a["geometry_id"],
        "tx": a["tx"], "ty": a["ty"], "tz": a["tz"], "tt": a["tt"],
        "tpx": a["tpx"], "tpy": a["tpy"], "tpz": a["tpz"], "te": a["te"],
        "deltapx": a["deltapx"], "deltapy": a["deltapy"],
        "deltapz": a["deltapz"], "deltae": a["deltae"],
        "index": a["index"], "volume_id": a["volume_id"], "layer_id": a["layer_id"],
    })
    return df


def read_meas(meas_t, lo, hi):
    s, e = entry_bounds(meas_t, "event_nr", lo, hi)
    a = meas_t.arrays(
        ["event_nr", "volume_id", "layer_id", "surface_id",
         "rec_loc0", "rec_loc1", "rec_gx", "rec_gy", "rec_gz",
         "true_x", "true_y", "true_z",
         "particles_vertex_primary", "particles_vertex_secondary",
         "particles_particle", "particles_generation", "particles_sub_particle"],
        entry_start=s, entry_stop=e, library="np",
    )
    df = pd.DataFrame({
        "event_id": a["event_nr"],
        "volume_id": a["volume_id"], "layer_id": a["layer_id"], "module_id": a["surface_id"],
        "local0": a["rec_loc0"], "local1": a["rec_loc1"],
        "global_x": a["rec_gx"], "global_y": a["rec_gy"], "global_z": a["rec_gz"],
        "_pid": pack_barcode(_first(a["particles_vertex_primary"]), _first(a["particles_vertex_secondary"]),
                             _first(a["particles_particle"]), _first(a["particles_generation"]),
                             _first(a["particles_sub_particle"])),
        "_tx": a["true_x"], "_ty": a["true_y"], "_tz": a["true_z"],
    })
    return df


def read_particles(parts_t, lo, hi):
    s, e = entry_bounds(parts_t, "event_id", lo, hi)
    a = parts_t.arrays(
        ["event_id", "particle_type", "process", "vx", "vy", "vz", "vt",
         "px", "py", "pz", "m", "q",
         "vertex_primary", "vertex_secondary", "particle", "generation", "sub_particle"],
        entry_start=s, entry_stop=e, library="np",
    )
    rows = {}
    for i, ev in enumerate(a["event_id"]):
        rows[int(ev)] = pd.DataFrame({
            "particle_id": pack_barcode(a["vertex_primary"][i], a["vertex_secondary"][i],
                                        a["particle"][i], a["generation"][i], a["sub_particle"][i]),
            "particle_type": a["particle_type"][i], "process": a["process"][i],
            "vx": a["vx"][i], "vy": a["vy"][i], "vz": a["vz"][i], "vt": a["vt"][i],
            "px": a["px"][i], "py": a["py"][i], "pz": a["pz"][i],
            "m": a["m"][i], "q": a["q"][i],
            # barcode sub-fields -> lets the reader flag primary vs secondary
            "particle_id_sv": a["vertex_secondary"][i],
            "particle_id_gen": a["generation"][i],
            "particle_id_subpart": a["sub_particle"][i],
        })
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root-dir", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--first", type=int, default=98000, help="first event id")
    p.add_argument("--n", type=int, default=10, help="number of events")
    args = p.parse_args()

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    lo, hi = args.first, args.first + args.n

    hits_t = uproot.open(f"{args.root_dir}/hits.root")["hits"]
    meas_t = uproot.open(f"{args.root_dir}/measurements.root")["measurements"]
    parts_t = uproot.open(f"{args.root_dir}/particles_simulation.root")["particles"]

    print(f"Writing detectors.csv (minimal) -> {out_dir}")
    write_minimal_detectors(hits_t, meas_t, out_dir)

    print(f"Reading events [{lo}, {hi}) ...")
    hits_all = read_hits(hits_t, lo, hi)
    meas_all = read_meas(meas_t, lo, hi)
    parts_by_ev = read_particles(parts_t, lo, hi)
    hg = dict(tuple(hits_all.groupby("event_id")))
    mg = dict(tuple(meas_all.groupby("event_id")))

    for out_idx, ev in enumerate(range(lo, hi)):
        if ev not in hg or ev not in parts_by_ev:
            print(f"  event {ev}: missing, skipped")
            continue
        nh, nm = write_event(parts_by_ev[ev],
                             hg[ev].reset_index(drop=True),
                             mg.get(ev, hits_all.iloc[:0]).reset_index(drop=True),
                             ev, out_idx, out_dir)
        print(f"  event {ev} -> {out_idx:09d}: {nh} hits, {nm} measurements")

    print(f"\nDone. {args.n} events written to {out_dir}")


if __name__ == "__main__":
    main()
