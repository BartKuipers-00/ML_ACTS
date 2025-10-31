"""
create_track_table <out_prefix> [--hits HITS_ROOT]

Reads `hits.root` (or the provided file) in the current directory, assumes a
single-particle run, and produces a per-hit table and a short summary of how
many unique layers were traversed per detector volume. If available, the
script will also infer the primary particle transverse momentum (Pt) and
pseudorapidity (eta) using momentum stored in `hits.root` or falling back to
`particles.root` / `particles_simulation.root`.

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
try:
    import uproot
    import numpy as np
except Exception as e:
    print("Missing Python packages or import error:", e)
    print("Please ensure uproot and numpy are installed.")
    sys.exit(2)

LONG_STRIP_RADII = np.array([820., 1020.])

# Known endcap disk z-positions (absolute values, mm) - from your detector note
PIXEL_ENDCAP_Z = np.array([600., 700., 820., 960., 1100., 1300., 1500.])
STRIP_ENDCAP_Z = np.array([1220., 1500., 1800., 2150., 2550., 2950.])

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


# Note: classification is driven solely by VOLUME_MAP. Heuristics removed.


def main():
    p = argparse.ArgumentParser(description='Create track table from hits.root')
    p.add_argument('out_prefix', help='output prefix (files written as <prefix>.txt)')
    p.add_argument('--hits', '-i', help='path to hits.root (default: ./hits.root)', default='hits.root')
    p.add_argument('--volume-map-only', action='store_true', help='Only use explicit VOLUME_MAP to classify volumes; do not apply geometric heuristics')
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

    # optional per-hit particle id and momentum in hits.root
    pid_hit_arr = np.asarray(t['particle_id'].array()) if 'particle_id' in t.keys() else None
    has_momentum_in_hits = all(k in t.keys() for k in ('tpx', 'tpy', 'tpz'))
    if has_momentum_in_hits:
        hit_tpx = np.asarray(t['tpx'].array())
        hit_tpy = np.asarray(t['tpy'].array())
        hit_tpz = np.asarray(t['tpz'].array())
    else:
        hit_tpx = hit_tpy = hit_tpz = None

    # Fallback: try to read a particles file nearby to map particle_id -> (px,py,pz)
    particles_map = None
    for alt in ('particles.root', 'particles_simulation.root', 'particles.root'):
        altp = Path(alt)
        if not altp.exists():
            continue
        try:
            pr = uproot.open(str(altp))
            pname = None
            for k in pr.keys():
                if k.lower().startswith('particles'):
                    pname = k
                    break
            if pname is None:
                continue
            ptree = pr[pname]
            if not all(k in ptree.keys() for k in ('particle_id', 'px', 'py', 'pz')):
                continue
            p_pid = np.asarray(ptree['particle_id'].array())
            p_px = np.asarray(ptree['px'].array())
            p_py = np.asarray(ptree['py'].array())
            p_pz = np.asarray(ptree['pz'].array())
            particles_map = {}
            for i, pid in enumerate(p_pid):
                if int(pid) not in particles_map:
                    particles_map[int(pid)] = (float(p_px[i]), float(p_py[i]), float(p_pz[i]))
            break
        except Exception:
            particles_map = None
            continue

    # Build simple rows list (avoid pandas to keep script light)
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

    # Enforce strict volume-id mapping: require that every encountered volume_id
    # exists in VOLUME_MAP. Do not fall back to geometric heuristics. If any
    # unknown volume ids are present, print them and exit with non-zero code.
    encountered_vids = set(grouped.keys())
    unknown_vids = sorted([v for v in encountered_vids if v not in VOLUME_MAP])
    if unknown_vids:
        print('Error: encountered volume_id(s) not present in VOLUME_MAP:', unknown_vids)
        print('Please extend VOLUME_MAP in create_track_table.py to include these volume ids.')
        sys.exit(2)

    for vid, sub in grouped.items():
        mask = [r for r in sub if r['module_id'] != 0]
        if mask:
            median_r = float(np.median([r['r_mm'] for r in mask]))
            median_abs_z = float(np.median([abs(r['z_mm']) for r in mask]))
        else:
            median_r = float(np.median([r['r_mm'] for r in sub]))
            median_abs_z = float(np.median([abs(r['z_mm']) for r in sub]))

        # Classification is driven solely by VOLUME_MAP. We already enforced
        # that all encountered volume_ids are present in VOLUME_MAP earlier.
        vid_int = int(vid)
        n_layers = len(set([r['layer_id'] for r in mask]))
        vm = VOLUME_MAP[vid_int]
        is_barrel = (vm['region'] == 'barrel')
        is_endcap = (vm['region'] == 'endcap')
        barrel_type = vm['type'] if is_barrel else None
        endcap_type = vm['type'] if is_endcap else None
        detname = vm['type'] if is_barrel else (vm['type'] + 'Endcap')
        # optionally compute nearest disk z for endcaps for informative output
        nearest_z = None
        if is_endcap:
            if vm['type'].lower().startswith('pixel'):
                nearest_z = float(PIXEL_ENDCAP_Z[np.argmin(np.abs(PIXEL_ENDCAP_Z - median_abs_z))])
            else:
                nearest_z = float(STRIP_ENDCAP_Z[np.argmin(np.abs(STRIP_ENDCAP_Z - median_abs_z))])

        vol_info[int(vid)] = {
            'detname': detname,
            'median_r': median_r,
            'median_abs_z': median_abs_z,
            'n_layers': n_layers,
            'is_barrel': is_barrel,
            'is_endcap': is_endcap,
            'barrel_type': barrel_type,
            'endcap_type': endcap_type,
            'disk_z': nearest_z,
        }

    # If possible, infer the primary particle (by majority of measurement hits)
    particle_summary = None
    measurement_orig_indices = [r['orig_index'] for r in rows if r['module_id'] != 0]
    if len(measurement_orig_indices) > 0 and pid_hit_arr is not None:
        measurement_mask_orig = np.zeros(len(vol), dtype=bool)
        measurement_mask_orig[measurement_orig_indices] = True
        sel_pids = pid_hit_arr[measurement_mask_orig]
        if sel_pids.size > 0:
            unique_pids, counts = np.unique(sel_pids, return_counts=True)
            primary_idx = np.argmax(counts)
            primary_pid = int(unique_pids[primary_idx])
            # prefer per-hit momentum if present
            if hit_tpx is not None:
                sel_px = hit_tpx[measurement_mask_orig]
                sel_py = hit_tpy[measurement_mask_orig]
                sel_pz = hit_tpz[measurement_mask_orig]
                mask_pid = (sel_pids == primary_pid)
                avg_px = float(np.mean(sel_px[mask_pid]))
                avg_py = float(np.mean(sel_py[mask_pid]))
                avg_pz = float(np.mean(sel_pz[mask_pid]))
            elif particles_map is not None and primary_pid in particles_map:
                avg_px, avg_py, avg_pz = particles_map[primary_pid]
            else:
                avg_px = avg_py = avg_pz = None

            if avg_px is not None:
                pt_val = float(np.hypot(avg_px, avg_py))
                p_val = float(np.sqrt(avg_px**2 + avg_py**2 + avg_pz**2))
                if p_val > abs(avg_pz):
                    eta_val = 0.5 * np.log((p_val + avg_pz) / (p_val - avg_pz))
                else:
                    eta_val = float('nan')
                particle_summary = {
                    'particle_id': primary_pid,
                    'n_hits': int(counts[primary_idx]),
                    'pt': pt_val,
                    'eta': eta_val,
                }

    # Build summary text
    # Build per-volume lines in the same order as the per-measurement table
    per_vol_lines = []
    seen_vids = []
    display_rows = [r for r in rows if r['module_id'] != 0]
    for rrow in display_rows:
        vid = int(rrow['volume_id'])
        if vid not in seen_vids:
            seen_vids.append(vid)

    # totals are intentionally omitted; per-volume summary follows measurement order
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
