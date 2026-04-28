#!/usr/bin/env python3
"""
Sweep r_max from 0.05 to 0.4 and plot precision/recall vs r_max.

Example:
    python sweep_rmax.py saved_models/low_pt_latentmodel_100_mixed_f1=0.2752.ckpt --num-events 50 --step 0.05 --segmented true --dr-same-layer-cut 10
    python sweep_rmax.py saved_models/NOT_SEGMENTED_latentmodel_f1=0.3062.ckpt --num-events 50 --step 0.05 --segmented false --r-max-geometric 1000 --dr-same-layer-cut 10
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent


def run_test(checkpoint_path, knn, r_max, r_max_geometric=None, dr_same_layer_cut=None, num_events=None, segmented=None):
    cmd = [
        sys.executable, str(SCRIPT_DIR / "test_my_latent_model.py"),
        str(checkpoint_path),
        "--knn", str(knn),
        "--r-max", str(r_max),
        "--device", "auto",
    ]
    if r_max_geometric is not None:
        cmd += ["--r-max-geometric", str(r_max_geometric)]
    cmd += ["--dr-same-layer-cut", str(dr_same_layer_cut or 0)]
    if num_events is not None:
        cmd += ["--num-events", str(num_events)]
    if segmented is not None:
        cmd += ["--segmented", str(segmented)]
    print(f"\n{'='*60}")
    print(f"  r_max={r_max}  k_max={knn}  r_max_geometric={r_max_geometric}  dr_same_layer_cut={dr_same_layer_cut}  segmented={segmented}")
    print(f"{'='*60}")

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR (exit {result.returncode})")
        if result.stdout.strip():
            print("--- stdout ---")
            print(result.stdout[-3000:])
        if result.stderr.strip():
            print("--- stderr ---")
            print(result.stderr[-2000:])
        return None, None, None

    precision = recall = seg_parsed = None
    for line in result.stdout.splitlines():
        m = re.match(r"\s*Precision:\s+([\d.]+)", line)
        if m:
            precision = float(m.group(1))
        m = re.match(r"\s*Recall:\s+([\d.]+)", line)
        if m:
            recall = float(m.group(1))
        m = re.search(r"segmented[=:]\s*(True|False)", line, re.IGNORECASE)
        if m:
            seg_parsed = m.group(1).lower() == "true"

    if precision is None or recall is None:
        print("  WARNING: could not parse precision/recall from output")
        print(result.stdout[-3000:])
        return None, None, None

    print(f"  => precision={precision:.4f} ({precision*100:.2f}%)  recall={recall:.4f} ({recall*100:.2f}%)  segmented={seg_parsed}")
    return precision, recall, seg_parsed


def main():
    parser = argparse.ArgumentParser(description="Sweep r_max and plot precision/recall")
    parser.add_argument("checkpoint_path", type=str, help="Path to model checkpoint (.ckpt)")
    parser.add_argument("--r-max-min", type=float, default=0.02, help="Minimum r_max (default: 0.05)")
    parser.add_argument("--r-max-max", type=float, default=0.4,  help="Maximum r_max (default: 0.8)")
    parser.add_argument("--step", type=float, default=0.02, help="Step size for r_max sweep (default: 0.05)")
    parser.add_argument("--knn", type=int, default=1200, help="Fixed k_max (default: 1000)")
    parser.add_argument("--r-max-geometric", type=float, default=None, help="Fixed r_max_geometric in mm (default: None = disabled)")
    parser.add_argument("--dr-same-layer-cut", type=float, default=None, help="Remove edges between hits with |Δr| < this threshold (mm); 0 = disabled (default: 0)")
    parser.add_argument("--num-events", type=int, default=None, help="Number of test events to evaluate per sweep point (default: all in testset)")
    parser.add_argument("--output-dir", type=str, default=None, help="Directory to save plots (default: data/visuals/rmax_sweep/)")
    parser.add_argument("--segmented", type=lambda x: x.lower() not in ('false', '0', 'no'), default=None, metavar='BOOL',
                        help="Override segmented flag passed to test script (default: read from graph_construction_latent.yaml)")
    args = parser.parse_args()

    checkpoint_path = Path(args.checkpoint_path)
    if not checkpoint_path.is_absolute():
        checkpoint_path = PIPELINE_ROOT / checkpoint_path

    # Build r_max values with float step using numpy to avoid float accumulation errors
    r_max_values = [round(v, 6) for v in np.arange(args.r_max_min, args.r_max_max + args.step * 0.5, args.step)]

    ckpt_stem = checkpoint_path.stem
    base_dir = Path(args.output_dir) if args.output_dir else PIPELINE_ROOT / "data" / "visuals" / "rmax_sweep"
    output_dir = base_dir / ckpt_stem
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Checkpoint: {checkpoint_path}")
    print(f"r_max sweep: {r_max_values}")
    layer_str = f"{args.dr_same_layer_cut} mm" if args.dr_same_layer_cut is not None else "disabled"
    events_str = str(args.num_events) if args.num_events is not None else "all (from config)"
    print(f"k_max={args.knn}, r_max_geometric={args.r_max_geometric} mm, dr_same_layer_cut={layer_str}, num_events={events_str}, segmented={args.segmented if args.segmented is not None else '(from config)'}")
    print(f"Output: {output_dir}")

    precisions = []
    recalls = []
    valid_rmax = []
    segmented = None  # resolved from first successful run

    for r_max in r_max_values:
        p, r, seg = run_test(checkpoint_path, args.knn, r_max, args.r_max_geometric, dr_same_layer_cut=args.dr_same_layer_cut, num_events=args.num_events, segmented=args.segmented)
        if p is not None:
            precisions.append(p)
            recalls.append(r)
            valid_rmax.append(r_max)
            if segmented is None:
                segmented = seg

    if not valid_rmax:
        print("\nERROR: no successful runs, cannot plot.")
        sys.exit(1)

    print(f"\n{'='*60}")
    print("SWEEP SUMMARY")
    print(f"{'='*60}")
    print(f"  segmented={segmented}  k_max={args.knn}  r_max_geometric={args.r_max_geometric} mm")
    print(f"{'r_max':>8}  {'precision':>10}  {'recall':>10}")
    print(f"{'-'*34}")
    for rv, p, r in zip(valid_rmax, precisions, recalls):
        print(f"{rv:>8.3f}  {p:>10.4f}  {r:>10.4f}")

    # --- Save raw data ---
    raw_data_path = output_dir / "sweep_data.tsv"
    with open(raw_data_path, 'w') as f:
        f.write(f"# segmented={segmented}  k_max={args.knn}  r_max_geometric={args.r_max_geometric} mm\n")
        f.write("r_max\tprecision\trecall\n")
        for rv, p, r in zip(valid_rmax, precisions, recalls):
            f.write(f"{rv}\t{p:.6f}\t{r:.6f}\n")
    print(f"\nSaved raw data: {raw_data_path}")

    def padded_ylim(values, pad=0.05):
        lo, hi = min(values), max(values)
        span = hi - lo or 0.01
        return lo - pad * span, hi + pad * span

    subtitle = f"k_max={args.knn}, r_max_geometric={args.r_max_geometric} mm, segmented={segmented}"

    # --- Precision plot ---
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(valid_rmax, precisions, marker='o', linewidth=2, color='steelblue')
    ax.set_xlabel("r_max (latent space radius)", fontsize=12)
    ax.set_ylabel("Precision", fontsize=12)
    ax.set_title(f"Precision vs r_max\n{subtitle}", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(*padded_ylim(precisions))
    fig.tight_layout()
    precision_path = output_dir / "precision_vs_rmax.png"
    fig.savefig(precision_path, dpi=150)
    plt.close(fig)
    print(f"\nSaved: {precision_path}")

    # --- Recall plot ---
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(valid_rmax, recalls, marker='o', linewidth=2, color='darkorange')
    ax.set_xlabel("r_max (latent space radius)", fontsize=12)
    ax.set_ylabel("Recall", fontsize=12)
    ax.set_title(f"Recall vs r_max\n{subtitle}", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(*padded_ylim(recalls))
    fig.tight_layout()
    recall_path = output_dir / "recall_vs_rmax.png"
    fig.savefig(recall_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {recall_path}")

    # --- Combined plot (dual y-axes) ---
    fig, ax1 = plt.subplots(figsize=(7, 4))
    ax1.plot(valid_rmax, precisions, marker='o', linewidth=2, color='steelblue', label='Precision')
    ax1.set_xlabel("r_max (latent space radius)", fontsize=12)
    ax1.set_ylabel("Precision", fontsize=12, color='steelblue')
    ax1.tick_params(axis='y', labelcolor='steelblue')
    ax1.set_ylim(*padded_ylim(precisions))
    ax2 = ax1.twinx()
    ax2.plot(valid_rmax, recalls, marker='s', linewidth=2, color='darkorange', label='Recall')
    ax2.set_ylabel("Recall", fontsize=12, color='darkorange')
    ax2.tick_params(axis='y', labelcolor='darkorange')
    ax2.set_ylim(*padded_ylim(recalls))
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=11)
    ax1.set_title(f"Precision & Recall vs r_max\n{subtitle}", fontsize=11)
    ax1.grid(True, alpha=0.3)
    fig.tight_layout()
    combined_path = output_dir / "precision_recall_vs_rmax.png"
    fig.savefig(combined_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {combined_path}")


if __name__ == "__main__":
    main()
