"""
create_track_table <out_prefix> [--hits HITS_ROOT]

Reads `hits.root` (or the provided file) in the current directory, assumes a
single-particle run, and produces a per-hit table and a short summary of how
many unique layers were traversed per detector volume. If available, the
script will also infer the primary particle transverse momentum (Pt) and
pseudorapidity (eta) using momentum stored in `hits.root`.

Output files (written to current working directory):
  <out_prefix>.txt   - short textual summary (this script intentionally
                      writes only the TXT file so it can be inspected easily)

Usage examples:
  python create_track_table first_table
  python create_track_table first_table --hits /path/to/hits.root

"""
import argparse
import sys
from pathlib import Path

# lazy imports (install if missing)

import uproot
import numpy as np


# Known endcap disk z-positions (absolute values, mm) - removed (not used)

# Optional explicit mapping from volume_id -> detector region. If a
# volume_id appears here it will be used (1:1 mapping). Extend this dict
# when you know the geometry volume numbering. Example entries below were
# inferred from sample runs in this workspace.
VOLUME_MAP = {
    # barrel volumes
    8:  {'type': 'Pixel',      'region': 'barrel'},
    13: {'type': 'ShortStrip', 'region': 'barrel'},
    17: {'type': 'LongStrip',  'region': 'barrel'},
    # endcap volumes (examples observed in forward run)
    7:  {'type': 'Pixel',      'region': 'endcap'},
    12: {'type': 'Strip',      'region': 'endcap'},
    # additional endcap volumes observed on the opposite side / other runs
    9:  {'type': 'Pixel',      'region': 'endcap'},
    14: {'type': 'Strip',      'region': 'endcap'},
}


# Classification is driven by `VOLUME_MAP`.


def main():
    p = argparse.ArgumentParser(description='Create track table from hits.root')
    p.add_argument('out_prefix', help='output prefix (files written as <prefix>.txt)')
    p.add_argument('--hits', '-i', help='path to hits.root (default: ./hits.root)', default='hits.root')
    args = p.parse_args()

    out_prefix = Path(args.out_prefix)
    hits_path = Path(args.hits)
    if not hits_path.exists():
        print(f"Could not find hits file: {hits_path}")
        sys.exit(1)

    # open file and read branches
    try:
        tree = uproot.open(str(hits_path))
    except Exception as e:
        print('Error opening ROOT file:', e)
        sys.exit(1)

    # find 'hits' tree (common name) - accept single tree
    tree_names = list(tree.keys())
    if len(tree_names) == 0:
        print('No objects found in file')
        sys.exit(1)

    # pick tree that starts with 'hits' if available
    chosen = None
    for k in tree_names:
        if k.lower().startswith('hits'):
            chosen = k
            break
    if chosen is None:
        chosen = tree_names[0]

    t = tree[chosen]

    # required branches
    required = ['volume_id', 'layer_id', 'sensitive_id', 'tx', 'ty', 'tz']
    missing = [r for r in required if r not in t.keys()]
    if missing:
        print('Missing required branches in', chosen, 'missing:', missing)
        sys.exit(1)

    vol = np.asarray(t['volume_id'].array())
    lay = np.asarray(t['layer_id'].array())
    mod = np.asarray(t['sensitive_id'].array())
    tx = np.asarray(t['tx'].array())
    ty = np.asarray(t['ty'].array())
    tz = np.asarray(t['tz'].array())

    r = np.sqrt(tx**2 + ty**2)

    # detect whether per-hit momentum exists in the file; don't read arrays
    # until (and unless) we actually need them to avoid copying large arrays.
    has_momentum_in_hits = all(k in t.keys() for k in ('tpx', 'tpy', 'tpz'))
    hit_tpx = hit_tpy = hit_tpz = None


    rows = []
    for i in range(len(vol)):
        rows.append({
            'orig_index': int(i),
            'volume_id': int(vol[i]),
            'layer_id': int(lay[i]),
            'module_id': int(mod[i]),
            'x_mm': float(tx[i]),
            'y_mm': float(ty[i]),
            'z_mm': float(tz[i]),
            'r_mm': float(r[i]),
        })

    # Sort by decreasing r (outer -> inner)
    rows = sorted(rows, key=lambda r_: r_['r_mm'], reverse=True)

    # Unique layer counts per volume (only count entries where module_id != 0)
    uniques = set()
    for row in rows:
        if row['module_id'] != 0:
            uniques.add((row['volume_id'], row['layer_id']))
    total_unique_layers = len(uniques)

    # Map volume_id to inferred detector type
    vol_info = {}
    from collections import defaultdict
    grouped = defaultdict(list)
    for row in rows:
        grouped[row['volume_id']].append(row)

    # Build per-volume info. If a volume_id is unknown in VOLUME_MAP, warn
    # but continue — classify it as 'Unknown' instead of aborting.
    for vid, sub in grouped.items():
        mask = [r for r in sub if r['module_id'] != 0]
        n_layers = len(set([r['layer_id'] for r in mask]))
        vid_int = int(vid)
        if vid_int in VOLUME_MAP:
            vm = VOLUME_MAP[vid_int]
            is_barrel = (vm['region'] == 'barrel')
            is_endcap = (vm['region'] == 'endcap')
            barrel_type = vm['type'] if is_barrel else None
            endcap_type = vm['type'] if is_endcap else None
            detname = vm['type'] if is_barrel else (vm['type'] + 'Endcap')
        else:
            print(f"Warning: volume_id {vid_int} not present in VOLUME_MAP; classifying as 'Unknown'")
            is_barrel = False
            is_endcap = False
            barrel_type = None
            endcap_type = None
            detname = 'Unknown'

        vol_info[int(vid)] = {
            'detname': detname,
            'n_layers': n_layers,
            'is_barrel': is_barrel,
            'is_endcap': is_endcap,
            'barrel_type': barrel_type,
            'endcap_type': endcap_type,
        }

    # Infer primary particle (single-particle run assumed) and compute Pt/eta from hit momenta if available
    particle_summary = None
    measurement_rows = [r for r in rows if r['module_id'] != 0]
    if len(measurement_rows) > 0:
        # determine primary particle id if particle_id branch exists
        primary_pid = None
        if 'particle_id' in t.keys():
            pid_arr = np.asarray(t['particle_id'].array())
            meas_indices = [r['orig_index'] for r in measurement_rows]
            sel_pids = pid_arr[meas_indices]
            unique_pids, counts = np.unique(sel_pids, return_counts=True)
            if unique_pids.size > 0:
                primary_pid = int(unique_pids[np.argmax(counts)])

        # compute average momentum from hits if momentum branches are present
        if has_momentum_in_hits:
            # lazy-read momentum arrays now that we know which indices we will use
            hit_tpx = np.asarray(t['tpx'].array())
            hit_tpy = np.asarray(t['tpy'].array())
            hit_tpz = np.asarray(t['tpz'].array())
            meas_indices = [r['orig_index'] for r in measurement_rows]
            sel_px = hit_tpx[meas_indices]
            sel_py = hit_tpy[meas_indices]
            sel_pz = hit_tpz[meas_indices]
            avg_px = float(np.mean(sel_px))
            avg_py = float(np.mean(sel_py))
            avg_pz = float(np.mean(sel_pz))
            pt_val = float(np.hypot(avg_px, avg_py))
            p_val = float(np.sqrt(avg_px**2 + avg_py**2 + avg_pz**2))
            # robust pseudorapidity calculation: avoid division by zero or log(0).
            # If p == |pz| (i.e. pt == 0) then eta should be signed infinity when
            # pz != 0 (particle going exactly along +z/-z) and NaN only if the
            # momentum is zero/invalid. Use a small epsilon for numerical safety.
            eps = 1e-12
            if p_val > abs(avg_pz) + eps:
                eta_val = 0.5 * np.log((p_val + avg_pz) / (p_val - avg_pz))
            else:
                if pt_val == 0.0 and avg_pz != 0.0:
                    eta_val = float('inf') if avg_pz > 0.0 else float('-inf')
                else:
                    eta_val = float('nan')
            particle_summary = {
                'particle_id': primary_pid if primary_pid is not None else 0,
                'n_hits': len(measurement_rows),
                'pt': pt_val,
                'eta': eta_val,
            }
        else:
            particle_summary = {
                'particle_id': primary_pid if primary_pid is not None else 0,
                'n_hits': len(measurement_rows),
                'pt': None,
                'eta': None,
            }

    # Build summary text

    per_vol_lines = []
    seen_vids = []
    display_rows = [r for r in rows if r['module_id'] != 0]
    for rrow in display_rows:
        vid = int(rrow['volume_id'])
        if vid not in seen_vids:
            seen_vids.append(vid)

   
    for vid in seen_vids:
        info = vol_info[int(vid)]
        det = info['detname']
        nlay = info['n_layers']
        parts = [f"volume {vid}:"]
        if info.get('is_barrel') and info.get('is_endcap'):
            parts.append(f"{info['barrel_type']} (barrel & {info['endcap_type']})")
        elif info.get('is_barrel'):
            parts.append(f"{info['barrel_type']} Barrel")
        elif info.get('is_endcap'):
            parts.append(f"{info['endcap_type']} Endcap")
        else:
            parts.append(det)

        per_vol_lines.append(', '.join(parts) + f", unique layers hit={nlay}")

    summary_lines = []
    if particle_summary is not None:
        summary_lines.append(f"Primary particle id={particle_summary['particle_id']} (contributing hits={particle_summary['n_hits']})")
        # Only display Pt/eta if we actually computed momentum values.
        if particle_summary.get('pt') is not None:
            # eta may be inf/-inf/nan; format will render these sensibly.
            summary_lines.append(f"  Pt = {particle_summary['pt']:.3f} , eta = {particle_summary['eta']:.3f}")
    summary_lines.append(f"Total unique layers entered: {total_unique_layers}")
    

    # Prepare TXT output only (summary + table)
    txt_path = Path(f"{out_prefix}.txt")
    with open(txt_path, 'w') as f:
        for line in summary_lines:
            f.write(line + "\n")
        f.write('\nPer-volume summary:\n')
        for l in per_vol_lines:
            f.write('  ' + l + '\n')

        f.write('\nPer-measurement table (measurements only, outer -> inner):\n')
        hdr = ('idx', 'orig_idx', 'vol', 'lay', 'mod', 'r_mm', 'x_mm', 'y_mm', 'z_mm')
        widths = [4, 8, 5, 5, 7, 9, 9, 9, 9]
        # header
        header_line = ' '.join(h.center(w) for h, w in zip(hdr, widths))
        f.write(header_line + '\n')
        f.write('-' * sum(widths) + '\n')
        # write only rows with module_id != 0
        display_rows = [r for r in rows if r['module_id'] != 0]
        for i, rrow in enumerate(display_rows):
            parts = [
                str(i).rjust(widths[0]),
                str(rrow['orig_index']).rjust(widths[1]),
                str(rrow['volume_id']).rjust(widths[2]),
                str(rrow['layer_id']).rjust(widths[3]),
                str(rrow['module_id']).rjust(widths[4]),
                f"{rrow['r_mm']:.1f}".rjust(widths[5]),
                f"{rrow['x_mm']:.1f}".rjust(widths[6]),
                f"{rrow['y_mm']:.1f}".rjust(widths[7]),
                f"{rrow['z_mm']:.1f}".rjust(widths[8]),
            ]
            f.write(' '.join(parts) + '\n')

    # Print short confirmation
    print('\n'.join(summary_lines))
    print('\nPer-volume summary:')
    for l in per_vol_lines:
        print(' ', l)
    print('\nWrote:', txt_path.resolve())


if __name__ == '__main__':
    main()
