#!/usr/bin/env python3
"""
Plot tracking efficiency vs particles_per_vertex across the OOD sweep datasets.

Overall efficiency comes from summary_testset.json.
Per-species efficiency is computed from particles_testset.csv
(particle_type, is_reconstructable, is_reconstructed columns).

Also saves ood_efficiency_data.json next to this script — a minimal snapshot
of all data used for the plot, so the plot can be recreated without the full
pipeline graphs/evaluation files.

Usage:
    python plot_ood_efficiency.py            # collect from pipeline + plot
    python plot_ood_efficiency.py --from-cache  # plot from saved JSON only
"""

import json
import sys
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

EXPERIMENT_DIR = Path(__file__).resolve().parent
WORKSPACE_ROOT = EXPERIMENT_DIR.parent
PIPELINE_ROOT = WORKSPACE_ROOT / "low_pt_gnn_pipeline"

with open(EXPERIMENT_DIR / "ood_experiment.yaml") as f:
    OOD = yaml.safe_load(f)

max_ppv = OOD["max_particles_per_vertex"]
PPV_VALUES = list(range(10, max_ppv + 1, 10))
MATCHER_SUFFIX = "segmentgnn" if OOD.get("use_segment_gnn", False) else "helix"
N_SPECIES = 4  # muon, pion, electron, proton

SPECIES_COLORS = {
    'Muon':     '#1f77b4',
    'Pion':     '#ff7f0e',
    'Electron': '#2ca02c',
    'Proton':   '#d62728',
}
SPECIES_MARKERS = {'Muon': 'o', 'Pion': 's', 'Electron': '^', 'Proton': 'D'}

CACHE_PATH = EXPERIMENT_DIR / "ood_efficiency_data.json"

# ── Load data ─────────────────────────────────────────────────────────────────

if "--from-cache" in sys.argv:
    with open(CACHE_PATH) as f:
        cache = json.load(f)
    x_all = cache["overall"]["x"]
    eff_all = cache["overall"]["eff"]
    species_data = cache["per_species"]
    print(f"Loaded from cache: {CACHE_PATH}")
else:
    x_all, eff_all = [], []
    species_data = {}   # species -> {x: [], eff: []}
    missing = []

    for ppv in PPV_VALUES:
        eval_dir = PIPELINE_ROOT / f"data_{ppv}" / f"track_evaluation_{MATCHER_SUFFIX}" / "testset"
        summary_path = eval_dir / "summary_testset.json"
        particles_path = eval_dir / "particles_testset.csv"

        if not summary_path.exists():
            missing.append(ppv)
            continue

        x = ppv * N_SPECIES  # total particles per event

        with open(summary_path) as f:
            summary = json.load(f)
        x_all.append(x)
        eff_all.append(summary["efficiency"] * 100)

        if particles_path.exists():
            df = pd.read_csv(particles_path)
            df = df[df["is_reconstructable"]]
            df = df.drop_duplicates(subset=["event_id", "particle_id"])
            for species, grp in df.groupby("particle_type"):
                eff = grp["is_reconstructed"].sum() / len(grp) * 100
                species_data.setdefault(species, {"x": [], "eff": []})
                species_data[species]["x"].append(x)
                species_data[species]["eff"].append(eff)

    if missing:
        print(f"Missing datasets (not yet evaluated): ppv = {missing}")
    if not x_all:
        print("No results found. Run the eval sweep first.")
        raise SystemExit(1)

    cache = {"overall": {"x": x_all, "eff": eff_all}, "per_species": species_data}
    with open(CACHE_PATH, "w") as f:
        json.dump(cache, f, indent=2)
    print(f"Saved plot data: {CACHE_PATH}")

fig, ax = plt.subplots(figsize=(8, 5))

# Overall efficiency
ax.plot(x_all, eff_all, marker="o", linewidth=2, markersize=6,
        color="black", linestyle="--", label="Overall")

# Per-species lines
for species, d in sorted(species_data.items()):
    ax.plot(d["x"], d["eff"],
            marker=SPECIES_MARKERS.get(species, "o"),
            color=SPECIES_COLORS.get(species, "gray"),
            linewidth=2, markersize=6, label=species)

ax.set_xlabel("Number of particles per event", fontsize=13)
ax.set_ylabel("Tracking efficiency [%]", fontsize=13)
ax.set_title("OOD Sweep: Tracking Efficiency vs Multiplicity", fontsize=14)
ax.set_xticks(x_all)
ax.set_ylim(0, 105)
ax.grid(True, alpha=0.4)
ax.legend(fontsize=11)

out_path = EXPERIMENT_DIR / "ood_efficiency_vs_multiplicity.pdf"
fig.savefig(out_path, bbox_inches="tight")
print(f"Saved: {out_path}")
