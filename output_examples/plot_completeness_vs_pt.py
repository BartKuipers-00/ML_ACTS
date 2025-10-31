#!/usr/bin/env python3
import pathlib
import sys
try:
    import uproot
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except Exception as e:
    print('Import error:', e)
    print('Make sure uproot, numpy, matplotlib are installed in this Python environment.')
    sys.exit(2)

p = pathlib.Path('performance_finding_ckf_matchingdetails.root')
if not p.exists():
    print('File not found:', p.resolve())
    sys.exit(1)

outdir = pathlib.Path('.')
print('Opening', p)
with uproot.open(str(p)) as f:
    if 'matchingdetails' not in f:
        print('matchingdetails tree not found in', p)
        sys.exit(1)
    tree = f['matchingdetails']
    print('Branches:', list(tree.keys()))
    # Read relevant branches
    want = ['track_index', 'track_pt', 'nMatchedHitsOnTrack', 'nParticleTruthHits', 'completeness', 'matched']
    want = [b for b in want if b in tree.keys()]
    arr = tree.arrays(want, library='np')

# Convert to arrays
track_idx = arr['track_index']
track_pt = arr['track_pt']
n_matched = arr['nMatchedHitsOnTrack'].astype(float)
n_particle = arr['nParticleTruthHits'].astype(float)
# pref existing completeness branch if present
if 'completeness' in arr:
    completeness_branch = arr['completeness'].astype(float)
else:
    completeness_branch = np.where(n_particle>0, n_matched / n_particle, np.nan)

matched_flag = arr['matched'] if 'matched' in arr else np.zeros_like(track_idx, dtype=bool)

# We'll compute one completeness per track. Prefer row with matched==True (majority),
# else take the row with maximum n_matched for that track.

from collections import defaultdict
rows_by_track = defaultdict(list)
for i, tid in enumerate(track_idx):
    rows_by_track[int(tid)].append(i)

track_list = []
for tid, indices in rows_by_track.items():
    # find majority if any
    maj = None
    for i in indices:
        if matched_flag[i]:
            maj = i
            break
    if maj is None:
        # choose index with largest n_matched (tie-breaker: largest completeness)
        best = max(indices, key=lambda j: (n_matched[j], completeness_branch[j] if not np.isnan(completeness_branch[j]) else -1))
        maj = best
    i = maj
    # only include tracks where particle truth hits > 0
    if n_particle[i] <= 0:
        comp = np.nan
    else:
        comp = n_matched[i] / n_particle[i]
    track_list.append((tid, track_pt[i], comp, n_matched[i], n_particle[i]))

# Convert to arrays and filter finite completeness
track_list = sorted(track_list, key=lambda x: x[1])
if len(track_list) == 0:
    print('No tracks found in matchingdetails')
    sys.exit(0)

track_arr = np.array(track_list, dtype=[('track', 'i8'), ('pt', 'f4'), ('completeness', 'f4'), ('n_matched', 'f4'), ('n_particle', 'f4')])
valid_mask = np.isfinite(track_arr['completeness'])
pt = track_arr['pt'][valid_mask]
comp = track_arr['completeness'][valid_mask]

print(f'Total tracks: {len(track_arr)}, with valid completeness: {len(pt)}')

# Scatter plot
plt.figure(figsize=(7,5))
plt.scatter(pt, comp, s=6, alpha=0.5)
plt.xscale('log')
plt.ylim(-0.05, 1.05)
plt.xlabel('track_pt')
plt.ylabel('completeness')
plt.title('Per-track completeness vs pT')
plt.grid(True, linestyle=':', alpha=0.4)
scatter_png = outdir / 'completeness_vs_pt_scatter.png'
plt.tight_layout()
plt.savefig(scatter_png)
plt.close()

# Binned median with 16/84 percentiles
# Define log-spaced bins between min and max pt (exclude pt<=0)
pt_pos = pt[pt>0]
comp_pos = comp[pt>0]
if pt_pos.size > 0:
    nbins = 10
    bins = np.logspace(np.log10(pt_pos.min()), np.log10(pt_pos.max()), nbins+1)
    bin_centers = np.sqrt(bins[:-1] * bins[1:])
    med = []
    p16 = []
    p84 = []
    counts = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        sel = (pt_pos >= lo) & (pt_pos < hi)
        selc = comp_pos[sel]
        counts.append(selc.size)
        if selc.size > 0:
            med.append(np.nanmedian(selc))
            p16.append(np.nanpercentile(selc, 16))
            p84.append(np.nanpercentile(selc, 84))
        else:
            med.append(np.nan)
            p16.append(np.nan)
            p84.append(np.nan)
    med = np.array(med)
    p16 = np.array(p16)
    p84 = np.array(p84)

    plt.figure(figsize=(7,5))
    plt.errorbar(bin_centers, med, yerr=[med - p16, p84 - med], fmt='o-', capsize=4)
    plt.xscale('log')
    plt.ylim(-0.05, 1.05)
    plt.xlabel('track_pt')
    plt.ylabel('median completeness (16/84 pct)')
    plt.title('Binned completeness vs pT')
    plt.grid(True, linestyle=':', alpha=0.4)
    binned_png = outdir / 'completeness_vs_pt_binned.png'
    plt.tight_layout()
    plt.savefig(binned_png)
    plt.close()
else:
    print('No positive pt values to bin; skipping binned plot')

# Save CSV summary
import csv
csv_file = outdir / 'completeness_per_track.csv'
with open(csv_file, 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['track_index', 'track_pt', 'completeness', 'n_matched_on_track', 'n_particle_truth_hits'])
    for row in track_arr:
        w.writerow([int(row['track']), float(row['pt']), float(row['completeness']) if np.isfinite(row['completeness']) else '', int(row['n_matched']), int(row['n_particle'])])

print('Wrote:', scatter_png, binned_png if pt_pos.size>0 else '(no binned)', csv_file)
print('Done.')
