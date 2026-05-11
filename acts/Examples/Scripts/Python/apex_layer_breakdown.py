#!/usr/bin/env python3
"""Per-track apex-layer arc-direction breakdown for Fatras simhits.

Categorises each truth particle by:
  - apex layer = largest-r barrel layer touched (via volume_id in {8,13,17}),
    where r is computed from the simhits' (tx,ty),
  - on the apex layer: how many hits have outward radial momentum
    (tpx*tx + tpy*ty > 0, i.e. trajectory still going outward),
    how many have inward (incoming arc).

For tracks below the V17L4 loop threshold (≈ 0.306 GeV pT), ideal physics is
1 outgoing + 1 incoming = clean loop. For tracks above the threshold, the
trajectory's apex sits outside the detector world volume and only the
outgoing arc is captured.

Usage:
  python3 apex_layer_breakdown.py /path/to/output_dir
or:
  python3 apex_layer_breakdown.py off=/path/A on=/path/B  (compare two runs)

Reads <dir>/hits.root and <dir>/particles_simulation.root.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import uproot


BARREL_VOLUMES = {8, 13, 17}  # Pixel, SStrip, LStrip in the Generic detector


def _load_hits_with_pt(run_dir: Path):
    """Return numpy arrays for hits + a (event,pid)->pt dict."""
    f = uproot.open(run_dir / "hits.root")
    htree = f[next(t for t in f.keys() if "hit" in t.lower())]
    hits = htree.arrays(
        [
            "event_id",
            "particle_id",
            "volume_id",
            "layer_id",
            "sensitive_id",
            "tx",
            "ty",
            "tz",
            "tpx",
            "tpy",
            "tpz",
        ],
        library="np",
    )

    p = uproot.open(run_dir / "particles_simulation.root")
    ptree = p[next(t for t in p.keys() if "particle" in t.lower())]
    parr = ptree.arrays(["event_id", "particle_id", "pt"], library="np")
    pt_map = {
        (int(e), int(pid)): float(pt)
        for e, pid, pt in zip(parr["event_id"], parr["particle_id"], parr["pt"])
    }
    return hits, pt_map


def _per_particle_apex_layer(hits):
    """For each (event, particle), return the apex barrel layer (vol, lay)
    and the list of hit indices on that layer."""
    # group hit indices by (event, particle)
    by_pid = defaultdict(list)
    for i, (ev, pid, vol) in enumerate(
        zip(hits["event_id"], hits["particle_id"], hits["volume_id"])
    ):
        if int(vol) in BARREL_VOLUMES:
            by_pid[(int(ev), int(pid))].append(i)

    apex = {}  # (event, pid) -> ((vol, lay), [hit_idx...])
    for key, idxs in by_pid.items():
        idxs_arr = np.asarray(idxs)
        vols = hits["volume_id"][idxs_arr]
        lays = hits["layer_id"][idxs_arr]
        tx = hits["tx"][idxs_arr]
        ty = hits["ty"][idxs_arr]
        r = np.sqrt(tx * tx + ty * ty)

        # for each (vol, lay) touched, mean r:
        layer_r = {}
        layer_idxs = defaultdict(list)
        for j, (v, l, ri) in enumerate(zip(vols, lays, r)):
            layer_idxs[(int(v), int(l))].append(idxs_arr[j])
            layer_r.setdefault((int(v), int(l)), []).append(float(ri))
        layer_meanr = {k: float(np.mean(v)) for k, v in layer_r.items()}
        apex_layer = max(layer_meanr, key=layer_meanr.get)
        apex[key] = (apex_layer, layer_idxs[apex_layer])
    return apex


def _classify(hits, apex):
    """Per-(event, pid): (n_out, n_in) on the apex layer."""
    classes = {}
    for key, (layer, idxs) in apex.items():
        if not idxs:
            classes[key] = (layer, 0, 0)
            continue
        idxs_arr = np.asarray(idxs)
        tx = hits["tx"][idxs_arr]
        ty = hits["ty"][idxs_arr]
        tpx = hits["tpx"][idxs_arr]
        tpy = hits["tpy"][idxs_arr]
        # radial momentum sign: positive = outgoing
        pr = tpx * tx + tpy * ty
        n_out = int((pr > 0).sum())
        n_in = int((pr < 0).sum())
        classes[key] = (layer, n_out, n_in)
    return classes


def _summary(classes, pt_map, label: str):
    """Group by apex layer + pT bin and report fractions."""
    print(f"\n=== {label} ===")
    print(f"Total particles touching ≥1 barrel layer: {len(classes)}")

    # Per-apex-layer aggregate
    by_layer = defaultdict(lambda: {"clean": 0, "miss_in": 0, "miss_out": 0,
                                     "over": 0, "total": 0})
    for (ev, pid), (layer, n_out, n_in) in classes.items():
        bucket = by_layer[layer]
        bucket["total"] += 1
        if n_out >= 2 or n_in >= 2:
            bucket["over"] += 1
        elif n_out == 1 and n_in == 1:
            bucket["clean"] += 1
        elif n_out == 1 and n_in == 0:
            bucket["miss_in"] += 1
        elif n_out == 0 and n_in == 1:
            bucket["miss_out"] += 1
        else:
            bucket["over"] += 1  # zero or weird

    print(f"\n{'apex layer':<14} {'total':>6} {'clean%':>8} "
          f"{'miss-in%':>9} {'miss-out%':>10} {'over%':>7}")
    for layer in sorted(by_layer):
        b = by_layer[layer]
        t = b["total"]
        if t == 0:
            continue
        print(f"  v{layer[0]:2d}l{layer[1]:1d}      {t:>6} "
              f"{100*b['clean']/t:>7.1f}% {100*b['miss_in']/t:>8.1f}% "
              f"{100*b['miss_out']/t:>9.1f}% {100*b['over']/t:>6.1f}%")

    # pT-binned breakdown across all apex layers
    bin_edges = np.arange(0.10, 0.51, 0.02)  # 20-MeV bins
    bins = defaultdict(lambda: {"clean": 0, "miss_in": 0, "miss_out": 0,
                                 "over": 0, "total": 0})
    for (ev, pid), (layer, n_out, n_in) in classes.items():
        pt = pt_map.get((ev, pid))
        if pt is None:
            continue
        bin_idx = int(np.clip(np.searchsorted(bin_edges, pt) - 1,
                              0, len(bin_edges) - 2))
        b = bins[bin_idx]
        b["total"] += 1
        if n_out >= 2 or n_in >= 2:
            b["over"] += 1
        elif n_out == 1 and n_in == 1:
            b["clean"] += 1
        elif n_out == 1 and n_in == 0:
            b["miss_in"] += 1
        elif n_out == 0 and n_in == 1:
            b["miss_out"] += 1
        else:
            b["over"] += 1
    print(f"\n{'pT range [GeV]':<18} {'total':>6} {'clean%':>8} "
          f"{'miss-in%':>9} {'miss-out%':>10} {'over%':>7}")
    for i in range(len(bin_edges) - 1):
        b = bins.get(i)
        if b is None or b["total"] == 0:
            continue
        lo, hi = bin_edges[i], bin_edges[i + 1]
        t = b["total"]
        print(f"  [{lo:.3f}, {hi:.3f})    {t:>6} "
              f"{100*b['clean']/t:>7.1f}% {100*b['miss_in']/t:>8.1f}% "
              f"{100*b['miss_out']/t:>9.1f}% {100*b['over']/t:>6.1f}%")


def analyse(run_dir: Path, label: str):
    hits, pt_map = _load_hits_with_pt(run_dir)
    apex = _per_particle_apex_layer(hits)
    classes = _classify(hits, apex)
    _summary(classes, pt_map, label)
    return classes


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("paths", nargs="+",
                        help="One or more run dirs. Use 'label=path' to label.")
    args = parser.parse_args()

    specs = []
    for p in args.paths:
        if "=" in p:
            label, path = p.split("=", 1)
        else:
            label, path = Path(p).name, p
        specs.append((label, Path(path)))

    for label, path in specs:
        analyse(path, label)


if __name__ == "__main__":
    main()
