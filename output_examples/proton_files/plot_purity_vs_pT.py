#!/usr/bin/env python3
"""
Compute purity vs pT from per-track tuples in `tracksummary_ckf.root` using uproot.

This reads `t_pT` and `matchingFraction` from the `tracksummary` tree and
computes a binned mean (like a TProfile) with configurable binning.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

try:
    import uproot
except Exception as e:
    raise SystemExit("Missing dependency 'uproot'. Install with 'pip install uproot' or use conda-forge.") from e

# Resolve files relative to this script
script_dir = Path(__file__).resolve().parent
tracksummary_file = script_dir / "tracksummary_ckf.root"
if not tracksummary_file.exists():
    candidate = script_dir / "output" / "tracksummary_ckf.root"
    if candidate.exists():
        tracksummary_file = candidate
    else:
        raise SystemExit(f"tracksummary_ckf.root not found (tried {tracksummary_file} and {candidate})")

# Configurable parameters
pt_max = 2.0   # from your config momentum.max
n_bins = 30
out_png = script_dir / "purity_vs_pT_from_tracksummary.png"

print(f"Opening {tracksummary_file}")
with uproot.open(str(tracksummary_file)) as f:
    # Tree name used by RootTrackSummaryWriter is 'tracksummary'
    if "tracksummary" not in f:
        # show keys to help debugging
        print("Keys in file:", f.keys())
        raise SystemExit("'tracksummary' TTree not found in file")
    tree = f["tracksummary"]

    # Branch names: t_pT (truth pT) and nMajorityHits / nMeasurements
    # Read as numpy arrays
    arr_pT = tree["t_pT"].array(library="np")
    arr_nMajorityHits = tree["nMajorityHits"].array(library="np")
    arr_nMeasurements = tree["nMeasurements"].array(library="np")

    # The tracksummary branches are stored as per-event vectors (arrays of tracks per event).
    # Flatten them into flat 1D arrays. Use awkward if available for robust flattening.
    try:
        import awkward as ak
        pT = ak.to_numpy(ak.flatten(arr_pT))
        nMajorityHits = ak.to_numpy(ak.flatten(arr_nMajorityHits))
        nMeasurements = ak.to_numpy(ak.flatten(arr_nMeasurements))
    except Exception:
        # fallback: concatenate object-dtype arrays
        try:
            pT = np.concatenate([np.asarray(x) for x in arr_pT if len(x) > 0])
            nMajorityHits = np.concatenate([np.asarray(x) for x in arr_nMajorityHits if len(x) > 0])
            nMeasurements = np.concatenate([np.asarray(x) for x in arr_nMeasurements if len(x) > 0])
        except Exception as e:
            raise SystemExit("Failed to flatten per-event arrays; install 'awkward' or ensure branches are flat") from e

    # Now pT, nMajorityHits and nMeasurements are flat numpy arrays
    pT = np.asarray(pT, dtype=float)
    nMajorityHits = np.asarray(nMajorityHits, dtype=np.uint64)
    nMeasurements = np.asarray(nMeasurements, dtype=np.uint64)

    # Filter invalid entries: require measurements > 0 and majority hit sentinel not present
    UINT32_MAX = np.iinfo(np.uint32).max
    valid_mask = (nMeasurements > 0) & (nMajorityHits < UINT32_MAX) & (nMajorityHits <= nMeasurements)
    pT = pT[valid_mask]
    purity = (nMajorityHits.astype(float) / nMeasurements.astype(float))[valid_mask]

if pT.size == 0:
    raise SystemExit("No valid (pT, purity) entries found in tracksummary tree")

# Bin the data and compute mean & stderr per bin (like TProfile)
bins = np.linspace(0.0, pt_max, n_bins + 1)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
bin_idx = np.digitize(pT, bins) - 1
means = np.full(n_bins, np.nan)
stderr = np.full(n_bins, np.nan)
counts = np.zeros(n_bins, dtype=int)

for i in range(n_bins):
    sel = bin_idx == i
    cnt = np.count_nonzero(sel)
    counts[i] = cnt
    if cnt > 0:
        vals = purity[sel]
        means[i] = np.mean(vals)
        stderr[i] = np.std(vals, ddof=1) / np.sqrt(cnt) if cnt > 1 else 0.0

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(bin_centers, means, yerr=stderr, xerr=(bins[1]-bins[0])/2.0, fmt='o', ms=4, color='tab:blue', ecolor='tab:gray', capsize=3)
ax.set_xlim(0.0, pt_max)
ax.set_xlabel("pT [GeV/c]")
ax.set_ylabel("Purity (matching fraction)")
ax.set_title("Purity vs pT (from tracksummary)")
ax.grid(True)
plt.tight_layout()
fig.savefig(str(out_png))
print(f"Wrote {out_png}")