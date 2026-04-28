#!/usr/bin/env python3
"""
Parameter tuning for helix-based segment matching.

Segments and helix fits are precomputed once per event; the sweep only
re-runs match_segments + segments_to_track_labels for each combination.

Usage:
    python tune_segment_matching.py                         # default: 20 val events
    python tune_segment_matching.py --n-events 20
    python tune_segment_matching.py --n-events 5 --score-cut 0.7

The parameter grid is defined in the PARAM_GRID dict near the top of this script.
Edit it to focus the sweep on the parameters you care about.
"""

import argparse
import copy
import sys
from pathlib import Path
from itertools import product

import pandas as pd
import torch
import yaml
from tqdm import tqdm

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / 'acorn'))
sys.path.insert(0, str(PIPELINE_ROOT))

from acorn.utils.loading_utils import load_datafiles_in_dir
from low_pt_custom_utils.segment_matching import (
    extract_segments_from_ground_truth,
    extract_segments_from_cc,
    fit_helices_to_segments,
    match_segments,
    segments_to_track_labels,
)
from low_pt_custom_utils.track_evaluation_utils import safe_evaluate_labelled_graph


# Parameter grid is defined in the config YAML under the `param_grid` key.


# ─── Helpers ─────────────────────────────────────────────────────────────────

def set_nested(d, dotted_key, value):
    """Set d['a']['b'] from 'a.b'."""
    keys = dotted_key.split('.')
    for k in keys[:-1]:
        d = d.setdefault(k, {})
    d[keys[-1]] = value


def precompute_segments(events, base_config):
    """
    Extract segments and fit helices once per event.
    Returns a list of (fitted_segments, graph) tuples.
    """
    use_gt = base_config.get("use_gt_segments", False)
    score_cut = base_config.get("score_cut", 0.5)
    B_field = base_config.get("B_field", 2.0)
    outlier_rejection = base_config.get("outlier_rejection", False)

    cached = []
    for graph in tqdm(events, desc="Precomputing segments + helix fits"):
        if use_gt:
            segments = extract_segments_from_ground_truth(graph)
        else:
            segments = extract_segments_from_cc(graph, score_cut)
        segments = fit_helices_to_segments(segments, graph, B_field=B_field, outlier_rejection=outlier_rejection)
        cached.append((segments, graph))
    return cached


def evaluate_cached(cached_events, matching_config, eval_config):
    """
    Run only match_segments + segments_to_track_labels on precomputed segments.
    Returns (efficiency, fake_rate, clone_rate).
    """
    sel_conf = eval_config.get('target_tracks', {})
    matching_fraction = eval_config.get('matching_fraction', 0.5)
    matching_style = eval_config.get('matching_style', 'ATLAS')
    min_track_length = eval_config.get('min_track_length', 1)

    dfs = []
    for segments, graph in cached_events:
        matched_tracks, unmatched = match_segments(segments, matching_config)
        num_nodes = graph.hit_x.size(0)
        hit_t = graph.hit_t.cpu().numpy() if hasattr(graph, 'hit_t') else None
        labels = segments_to_track_labels(matched_tracks, unmatched, num_nodes, hit_t=hit_t)

        g = copy.copy(graph)
        g.hit_track_labels = labels

        df = safe_evaluate_labelled_graph(
            g,
            sel_conf=sel_conf,
            matching_fraction=matching_fraction,
            matching_style=matching_style,
            min_track_length=min_track_length,
        )
        dfs.append(df)

    ev = pd.concat(dfs, ignore_index=True)

    particles = ev[ev["is_reconstructable"]]
    reco = particles[particles["is_reconstructed"] & particles["is_matchable"]]
    tracks = ev[ev["is_matchable"]]
    matched = tracks[tracks["is_matched"]]

    n_p   = len(particles.drop_duplicates(subset=["event_id", "particle_id"]))
    n_r   = len(reco.drop_duplicates(subset=["event_id", "particle_id"]))
    n_t   = len(tracks.drop_duplicates(subset=["event_id", "track_id"]))
    n_m   = len(matched.drop_duplicates(subset=["event_id", "track_id"]))
    n_dup = len(reco) - n_r

    eff   = n_r / n_p if n_p > 0 else 0.0
    fake  = 1 - (n_m / n_t) if n_t > 0 else 0.0
    clone = n_dup / n_r if n_r > 0 else 0.0

    return eff, fake, clone


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Sweep segment matching parameters')
    parser.add_argument('--n-events', type=int, default=20,
                        help='Number of validation events to use (default: 20)')
    parser.add_argument('--input-dir', type=str, default='data_300_mixed/gnn_stage',
                        help='Input data directory')
    parser.add_argument('--score-cut', type=float, default=None,
                        help='Override score cut (default: from config)')
    parser.add_argument('--config', type=str, default=None,
                        help='Path to config YAML (default: acorn_configs/...)')
    parser.add_argument('--top', type=int, default=20,
                        help='Show top N results sorted by efficiency (default: 20)')
    args = parser.parse_args()

    # ── Load config ──────────────────────────────────────────────────────────
    if args.config is None:
        config_file = PIPELINE_ROOT / 'acorn_configs' / 'track_building_stage_(3)' / 'tune_grid_helix.yaml'
    else:
        config_file = Path(args.config)

    with open(config_file) as f:
        base_config = yaml.safe_load(f)

    if args.score_cut is not None:
        base_config['score_cut'] = args.score_cut

    # ── Load graphs once ──────────────────────────────────────────────────────
    input_dir = Path(args.input_dir)
    if not input_dir.is_absolute():
        input_dir = PIPELINE_ROOT / input_dir

    paths = sorted(load_datafiles_in_dir(str(input_dir), 'valset', args.n_events))[:args.n_events]
    if not paths:
        print(f"ERROR: no valset events found in {input_dir}/valset/")
        sys.exit(1)

    print(f"Loading {len(paths)} validation events ...")
    events = [torch.load(p, map_location='cpu', weights_only=False) for p in tqdm(paths)]

    # ── Precompute segments + helix fits (done once) ──────────────────────────
    cached_events = precompute_segments(events, base_config)

    # ── Sweep ─────────────────────────────────────────────────────────────────
    param_grid = base_config.get('param_grid', {})
    if not param_grid:
        raise ValueError("No param_grid found in config — define it under the 'param_grid' key.")
    param_keys = list(param_grid.keys())
    param_values = list(param_grid.values())
    combos = list(product(*param_values))
    print(f"Done. Running parameter sweep over {len(combos)} combinations...\n")

    results = []
    for combo in tqdm(combos, desc='Sweeping'):
        cfg = copy.deepcopy(base_config)
        for key, val in zip(param_keys, combo):
            set_nested(cfg, key, val)

        # Auto-couple weight_pitch = 1 - weight_center - weight_R (floored at 0)
        w_c = cfg['matching'].get('weight_center', 0.5)
        w_r = cfg['matching'].get('weight_R', 0.0)
        cfg['matching']['weight_pitch'] = round(max(0.0, 1.0 - w_c - w_r), 4)

        eff, fake, clone = evaluate_cached(cached_events, cfg['matching'], cfg)
        results.append({
            **dict(zip(param_keys, combo)),
            'matching.weight_pitch': cfg['matching']['weight_pitch'],
            'eff': eff, 'fake': fake, 'clone': clone,
        })

    # ── Report ────────────────────────────────────────────────────────────────
    df = pd.DataFrame(results).sort_values('eff', ascending=False)

    short = {k: k.replace('matching.', '') for k in param_keys + ['matching.weight_pitch']}
    df_show = df.rename(columns=short).head(args.top)
    df_show['eff']   = df_show['eff'].map('{:.3f}'.format)
    df_show['fake']  = df_show['fake'].map('{:.3f}'.format)
    df_show['clone'] = df_show['clone'].map('{:.3f}'.format)

    print(f"\nTop {args.top} results (sorted by efficiency):")
    print(df_show.to_string(index=False))

    best = df.iloc[0]
    print(f"\n{'='*60}")
    print(f"Best: eff={best['eff']:.3f}  fake={best['fake']:.3f}  clone={best['clone']:.3f}")
    print(f"  matching:")
    for k in param_keys:
        print(f"    {k.replace('matching.', '')}: {best[k]}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
