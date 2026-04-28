#!/usr/bin/env python3
"""
Find optimal score_cut for CC clustering at the segment level.

For each score_cut in a range, runs CC clustering on the already-inferred
validation set and matches CC clusters to ground-truth segments (unique
particle_id × segment_id pairs). Reports efficiency (recall) and precision
at the segment level.

Usage:
    python find_cc_threshold.py --num_events 200 --range 0.6 1.0 --stepsize 0.02
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sps
import torch
from torch_geometric.utils import remove_isolated_nodes, to_scipy_sparse_matrix
from tqdm import tqdm

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / 'acorn'))
sys.path.insert(0, str(PIPELINE_ROOT))

from acorn.utils.loading_utils import load_datafiles_in_dir


# ─── CC Clustering ───────────────────────────────────────────────────────────

def run_cc(graph, score_cut: float):
    """
    Run Connected Components on edges with score > score_cut.

    Returns:
        labels: np.ndarray of shape (num_nodes,) with cluster id per node,
                -1 for isolated nodes.
    """
    scores = graph.edge_scores if hasattr(graph, 'edge_scores') else graph.scores
    edge_mask = scores > score_cut
    edges = graph.edge_index[:, edge_mask]
    num_nodes = graph.hit_x.size(0)

    edges_clean, _, node_mask = remove_isolated_nodes(edges, num_nodes=num_nodes)
    n_connected = node_mask.sum().item()

    labels = np.full(num_nodes, -1, dtype=np.int64)

    if n_connected == 0:
        return labels

    sparse_edges = to_scipy_sparse_matrix(edges_clean, num_nodes=n_connected)
    n_components, candidate_labels = sps.csgraph.connected_components(
        sparse_edges, directed=False, return_labels=True
    )

    labels[node_mask.numpy()] = candidate_labels
    return labels


# ─── Segment-Level Evaluation ────────────────────────────────────────────────

def evaluate_event(graph, score_cut: float, match_fraction: float = 0.5):
    """
    For one event, run CC at score_cut and evaluate against GT segments.

    GT segment = unique (particle_id, segment_id) pair with particle_id > 0.

    A GT segment is *found* (matched) if some CC cluster contains at least
    `match_fraction` of that segment's hits.

    A CC cluster is *pure* if at least `match_fraction` of its hits belong
    to a single GT segment (i.e. it does not span multiple segments heavily).

    Returns:
        n_gt_segs  : number of GT segments in this event
        n_found    : GT segments matched by at least one CC cluster
        n_cc_segs  : number of CC clusters with >= 2 hits (single-hit clusters ignored)
        n_pure     : pure CC clusters
    """
    particle_ids = graph.hit_particle_id.cpu().numpy().astype(np.int64)
    segment_ids  = graph.hit_segment_id.cpu().numpy().astype(np.int64)

    # ── Build GT segment map using tuple keys to avoid int64 overflow ──
    # (ACTS particle IDs are large 64-bit values; pid * 10000 overflows)
    signal_mask    = particle_ids > 0
    signal_indices = np.where(signal_mask)[0]

    pair_to_seg: dict[tuple, int] = {}
    gt_segments: dict[int, list] = {}
    hit_to_gt = np.full(len(particle_ids), -1, dtype=np.int64)

    for idx in signal_indices:
        key = (int(particle_ids[idx]), int(segment_ids[idx]))
        if key not in pair_to_seg:
            seg_idx = len(pair_to_seg)
            pair_to_seg[key] = seg_idx
            gt_segments[seg_idx] = []
        seg_idx = pair_to_seg[key]
        gt_segments[seg_idx].append(int(idx))
        hit_to_gt[idx] = seg_idx  # always non-negative for signal hits

    n_gt_segs = len(gt_segments)

    # ── Run CC ──
    cc_labels = run_cc(graph, score_cut)

    # ── Build CC cluster map: cluster_id → hit indices (size >= 2 only) ──
    unique_cc, counts = np.unique(cc_labels[cc_labels >= 0], return_counts=True)
    cc_clusters: dict[int, np.ndarray] = {}
    for cid, cnt in zip(unique_cc, counts):
        if cnt >= 2:
            cc_clusters[int(cid)] = np.where(cc_labels == cid)[0]

    n_cc_segs = len(cc_clusters)

    if n_gt_segs == 0 or n_cc_segs == 0:
        return n_gt_segs, 0, n_cc_segs, 0

    # ── For each CC cluster, build a hit → cc_id lookup ──
    hit_to_cc = np.full(len(particle_ids), -1, dtype=np.int64)
    for cc_id, hits in cc_clusters.items():
        hit_to_cc[hits] = cc_id

    # ── Efficiency: how many GT segments are found? ──
    n_found = 0
    for gt_hits in gt_segments.values():
        cc_ids_for_gt = hit_to_cc[gt_hits]
        for cc_id in np.unique(cc_ids_for_gt[cc_ids_for_gt >= 0]):
            overlap = np.sum(cc_ids_for_gt == cc_id)
            if overlap / len(gt_hits) >= match_fraction:
                n_found += 1
                break

    # ── Purity: how many CC clusters are predominantly one GT segment? ──
    # hit_to_gt is -1 only for true noise hits (pid=0), never for signal hits
    n_pure = 0
    for cc_hits in cc_clusters.values():
        gt_ids_in_cc = hit_to_gt[cc_hits]
        signal_hits_in_cc = gt_ids_in_cc[gt_ids_in_cc >= 0]
        if len(signal_hits_in_cc) == 0:
            continue  # genuine all-noise cluster
        _, cnts = np.unique(signal_hits_in_cc, return_counts=True)
        if cnts.max() / len(signal_hits_in_cc) >= match_fraction:
            n_pure += 1

    # ── Contamination: GT segments whose best-match CC cluster has foreign hits ──
    # For each GT segment, find the CC cluster with the most overlap.
    # If that cluster also contains hits from a different GT segment, it's contaminated.
    n_contaminated = 0
    n_has_match = 0
    for seg_idx, gt_hits in gt_segments.items():
        cc_ids_for_gt = hit_to_cc[gt_hits]
        valid = cc_ids_for_gt[cc_ids_for_gt >= 0]
        if len(valid) == 0:
            continue  # GT segment has no CC cluster at all
        n_has_match += 1
        unique_ids, cnts = np.unique(valid, return_counts=True)
        best_cc_id = int(unique_ids[np.argmax(cnts)])
        best_cc_hits = cc_clusters[best_cc_id]
        # Check for foreign hits: signal hits in the CC cluster not from this GT segment
        foreign = np.any((hit_to_gt[best_cc_hits] >= 0) & (hit_to_gt[best_cc_hits] != seg_idx))
        if foreign:
            n_contaminated += 1

    return n_gt_segs, n_found, n_cc_segs, n_pure, n_contaminated, n_has_match


# ─── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Sweep score_cut and evaluate CC clustering at the segment level.'
    )
    parser.add_argument(
        '--range', nargs=2, type=float, metavar=('A', 'B'), default=[0.5, 0.95],
        help='Score cut range [A, B] inclusive (default: 0.5 0.95)',
    )
    parser.add_argument(
        '--stepsize', type=float, default=0.05,
        help='Step size for score_cut sweep (default: 0.05)',
    )
    parser.add_argument(
        '--num_events', type=int, default=None,
        help='Number of validation events to use (default: all)',
    )
    parser.add_argument(
        '--match_fraction', type=float, default=0.5,
        help='Minimum hit-overlap fraction to count a segment as matched (default: 0.5)',
    )
    parser.add_argument(
        '--dataset', type=str, default='valset',
        help='Dataset split to evaluate (default: valset)',
    )
    args = parser.parse_args()

    score_cuts = np.arange(args.range[0], args.range[1] + args.stepsize * 0.5, args.stepsize)
    score_cuts = np.round(score_cuts, 6)

    input_dir = PIPELINE_ROOT / 'data' / 'gnn_stage'
    all_paths = sorted(load_datafiles_in_dir(str(input_dir), args.dataset))

    if args.num_events is not None:
        all_paths = all_paths[:args.num_events]

    n_events = len(all_paths)
    print(f"Evaluating {n_events} events from {input_dir / args.dataset}")
    print(f"Score cuts: {score_cuts.tolist()}")
    print(f"Match fraction: {args.match_fraction}\n")

    # Load all graphs once to avoid re-loading for every score cut
    print("Loading graphs...")
    graphs = []
    for path in tqdm(all_paths):
        g = torch.load(path, map_location='cpu', weights_only=False)
        graphs.append(g)

    # Check attributes
    sample = graphs[0]
    if not hasattr(sample, 'edge_scores') and not hasattr(sample, 'scores'):
        print("ERROR: graphs have no 'edge_scores' or 'scores' attribute.")
        print("       Run infer_gnn.py first to generate edge scores.")
        sys.exit(1)
    if not hasattr(sample, 'hit_segment_id'):
        print("ERROR: graphs have no 'hit_segment_id' attribute.")
        print("       Run clean_loops_and_attribute_segments.py first.")
        sys.exit(1)

    # ── One-event debug: what's inside impure CC clusters? ──────────────────
    g0 = graphs[0]
    pid0  = g0.hit_particle_id.numpy().astype(np.int64)
    sid0  = g0.hit_segment_id.numpy().astype(np.int64)
    sig0  = pid0 > 0

    pair_to_seg0: dict[tuple, int] = {}
    gt0: dict[int, list] = {}
    hit_to_gt0 = np.full(len(pid0), -1, dtype=np.int64)
    for idx in np.where(sig0)[0]:
        key = (int(pid0[idx]), int(sid0[idx]))
        if key not in pair_to_seg0:
            pair_to_seg0[key] = len(pair_to_seg0)
            gt0[pair_to_seg0[key]] = []
        si = pair_to_seg0[key]
        gt0[si].append(int(idx))
        hit_to_gt0[idx] = si

    cc0 = run_cc(g0, 0.5)
    u0, c0 = np.unique(cc0[cc0 >= 0], return_counts=True)
    cc_clusters0 = {int(cid): np.where(cc0 == cid)[0] for cid, cnt in zip(u0, c0) if cnt >= 2}

    all_noise_count = 0
    impure_counts = []
    for cc_hits in cc_clusters0.values():
        gt_ids = hit_to_gt0[cc_hits]
        sig_ids = gt_ids[gt_ids >= 0]
        if len(sig_ids) == 0:
            all_noise_count += 1
            continue
        _, cnts = np.unique(sig_ids, return_counts=True)
        purity = cnts.max() / len(sig_ids)
        if purity < 0.5:
            n_segs_in_cluster = len(cnts)
            impure_counts.append((len(cc_hits), len(sig_ids), n_segs_in_cluster, round(purity, 3)))

    print(f"\nDebug (event 0, score_cut=0.5): {len(gt0)} GT segs, {len(cc_clusters0)} CC clusters")
    n_true_noise = int((pid0 == 0).sum())
    n_signal_hits = int(sig0.sum())
    unassigned = int((hit_to_gt0 == -1).sum()) - n_true_noise  # signal hits missing from GT

    # Sample what's in an all-noise cluster
    all_noise_sample_pids = []
    for cc_hits in cc_clusters0.values():
        gt_ids = hit_to_gt0[cc_hits]
        sig_ids = gt_ids[gt_ids >= 0]
        if len(sig_ids) == 0:
            all_noise_sample_pids.append(pid0[cc_hits].tolist())
            if len(all_noise_sample_pids) >= 3:
                break

    print(f"  All-noise CC clusters: {all_noise_count}")
    print(f"  Impure CC clusters: {len(impure_counts)}")
    print(f"  Total hits: {len(pid0)}  signal (pid>0): {n_signal_hits}  true noise (pid=0): {n_true_noise}")
    print(f"  Signal hits NOT in any GT segment: {unassigned}  (bug if > 0)")
    print(f"  Sample pid values inside all-noise clusters: {all_noise_sample_pids}")
    if impure_counts:
        impure_counts.sort(key=lambda x: -x[2])
        print(f"  {'total_hits':>10}  {'signal_hits':>11}  {'n_gt_segs':>9}  {'purity':>6}")
        for row in impure_counts[:15]:
            print(f"  {row[0]:>10}  {row[1]:>11}  {row[2]:>9}  {row[3]:>6}")
    print(f"  Unique segment_ids in event 0: {np.unique(sid0[sig0]).tolist()}")
    print(f"  Particles with >1 segment: "
          f"{sum(1 for p in np.unique(pid0[sig0]) if len(np.unique(sid0[pid0 == p])) > 1)}"
          f" / {len(np.unique(pid0[sig0]))} total particles")
    print()

    # ── Edge score diagnostics: intra-segment vs cross-segment ──────────────
    intra_all, cross_all, noise_all = [], [], []
    for g in graphs:
        scores = (g.edge_scores if hasattr(g, 'edge_scores') else g.scores).numpy()
        pid = g.hit_particle_id.numpy()
        sid = g.hit_segment_id.numpy()
        ei  = g.edge_index.numpy()

        src_pid, dst_pid = pid[ei[0]], pid[ei[1]]
        src_sid, dst_sid = sid[ei[0]], sid[ei[1]]

        same_particle = src_pid == dst_pid
        signal        = same_particle & (src_pid > 0)
        same_segment  = src_sid == dst_sid

        intra_all.append(scores[signal & same_segment])
        cross_all.append(scores[signal & ~same_segment])
        noise_all.append(scores[~same_particle])

    intra = np.concatenate(intra_all)
    cross = np.concatenate(cross_all)
    noise = np.concatenate(noise_all)

    def score_stats(arr, label):
        if len(arr) == 0:
            print(f"  {label}: no edges")
            return
        def pct(x): return f"{x:.4%}"
        print(f"  {label} ({len(arr):,} edges): "
              f"mean={arr.mean():.3f}  median={np.median(arr):.3f}  "
              f">0.5: {pct((arr > 0.5).mean())}  >0.7: {pct((arr > 0.7).mean())}  "
              f">0.9: {pct((arr > 0.9).mean())}  "
              f"  abs>0.5: {int((arr > 0.5).sum() / n_events)}/ev  "
              f">0.9: {int((arr > 0.9).sum() / n_events)}/ev")

    print("\nEdge score diagnostics (averaged over all events):")
    score_stats(intra, "intra-segment (same particle, same segment)")
    score_stats(cross, "cross-segment (same particle, diff segment)")
    score_stats(noise, "noise         (different particles)         ")
    print()

    # Header
    col_w = 10
    print(f"{'score_cut':>{col_w}}  {'efficiency':>{col_w}}  {'precision':>{col_w}}"
          f"  {'contam_rate':>{col_w}}  {'gt_segs':>{col_w}}  {'cc_segs':>{col_w}}")
    print("-" * (col_w * 6 + 10))

    efficiencies, precisions, contam_rates = [], [], []
    for sc in score_cuts:
        total_gt = total_found = total_cc = total_pure = total_contam = total_has_match = 0
        for g in graphs:
            n_gt, n_found, n_cc, n_pure, n_contam, n_has_match = evaluate_event(g, float(sc), args.match_fraction)
            total_gt        += n_gt
            total_found     += n_found
            total_cc        += n_cc
            total_pure      += n_pure
            total_contam    += n_contam
            total_has_match += n_has_match

        efficiency  = total_found   / total_gt        if total_gt > 0        else float('nan')
        precision   = total_pure    / total_cc        if total_cc > 0        else float('nan')
        contam_rate = total_contam  / total_has_match if total_has_match > 0 else float('nan')

        efficiencies.append(efficiency)
        precisions.append(precision)
        contam_rates.append(contam_rate)
        print(f"{sc:>{col_w}.4f}  {efficiency:>{col_w}.4f}  {precision:>{col_w}.4f}"
              f"  {contam_rate:>{col_w}.4f}  {total_gt / n_events:>{col_w}.1f}  {total_cc / n_events:>{col_w}.1f}")

    # ── Plot efficiency vs contamination rate trade-off ──────────────────────
    valid       = np.array([not (np.isnan(e) or np.isnan(c)) for e, c in zip(efficiencies, contam_rates)])
    plot_cuts   = score_cuts[valid]
    plot_eff    = np.array(efficiencies)[valid]
    plot_prec   = np.array(precisions)[valid]
    plot_contam = np.array(contam_rates)[valid]

    fig, ax1 = plt.subplots(figsize=(8, 5))
    color_eff    = '#1f77b4'  # blue
    color_prec   = '#d62728'  # red
    color_contam = '#ff7f0e'  # orange

    ax1.plot(plot_cuts, plot_eff,   color=color_eff,    marker='o', label='Efficiency (recall)')
    ax1.set_xlabel('score_cut')
    ax1.set_ylabel('Efficiency (segment recall)', color=color_eff)
    ax1.tick_params(axis='y', labelcolor=color_eff)
    ax1.set_ylim(0, 1.05)

    ax2 = ax1.twinx()
    ax2.plot(plot_cuts, plot_prec,   color=color_prec,   marker='s', label='Precision (cluster purity)')
    ax2.plot(plot_cuts, plot_contam, color=color_contam, marker='^', label='Contamination rate')
    ax2.set_ylabel('Precision / Contamination rate', color='black')
    ax2.tick_params(axis='y', labelcolor='black')
    ax2.set_ylim(0, 1.05)

    ax1.set_title(f'CC segment metrics vs score_cut  ({n_events} events, match_fraction={args.match_fraction})')
    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [l.get_label() for l in lines], loc='center left')

    plot_path = PIPELINE_ROOT / 'data' / 'visuals' / 'cc_threshold_tradeoff.png'
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=150)
    print(f"\nPlot saved to {plot_path}")

    # ── TP vs FP edge confidence summary ────────────────────────────────────
    # intra/cross/noise arrays already built above; FP = cross-segment + noise
    fp = np.concatenate([cross, noise])

    print("\nTP vs FP edge confidence summary (across all evaluated events):")
    col_w2 = 12
    print(f"{'score_cut':>{col_w2}}  {'TP_passing':>{col_w2}}  {'FP_passing':>{col_w2}}"
          f"  {'TP_mean_score':>{col_w2}}  {'FP_mean_score':>{col_w2}}  {'FP_rate':>{col_w2}}")
    print("-" * (col_w2 * 6 + 10))

    for sc in score_cuts:
        tp_pass = intra[intra > sc]
        fp_pass = fp[fp > sc]
        n_tp = len(tp_pass)
        n_fp = len(fp_pass)
        tp_mean = tp_pass.mean() if n_tp > 0 else float('nan')
        fp_mean = fp_pass.mean() if n_fp > 0 else float('nan')
        fp_rate = n_fp / (n_tp + n_fp) if (n_tp + n_fp) > 0 else float('nan')
        print(f"{sc:>{col_w2}.4f}  {n_tp / n_events:>{col_w2}.1f}  {n_fp / n_events:>{col_w2}.2f}"
              f"  {tp_mean:>{col_w2}.4f}  {fp_mean:>{col_w2}.4f}  {fp_rate:>{col_w2}.4%}")


if __name__ == '__main__':
    main()
