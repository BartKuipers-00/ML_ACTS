#!/usr/bin/env python3
"""
Sweep k_max from 50 to 400 un step ... and plot precision/recall vs k_max.

Usage:
  

Example:
    python sweep_knn.py saved_models/low_pt_latentmodel_100_mixed_f1=0.2752.ckpt --num-events 50 --step 25 --segmented true --dr-same-layer-cut 10
    python sweep_knn.py saved_models/NOT_SEGMENTED_latentmodel_f1=0.3062.ckpt --num-events 50 --step 25 --segmented false --r-max-geometric 1000
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt

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
    print(f"  k_max={knn}  r_max={r_max}  r_max_geometric={r_max_geometric}  dr_same_layer_cut={dr_same_layer_cut}  segmented={segmented}")
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

    # Parse aggregate precision, recall, and segmented from output
    precision = recall = segmented = None
    for line in result.stdout.splitlines():
        m = re.match(r"\s*Precision:\s+([\d.]+)", line)
        if m:
            precision = float(m.group(1))
        m = re.match(r"\s*Recall:\s+([\d.]+)", line)
        if m:
            recall = float(m.group(1))
        m = re.search(r"segmented[=:]\s*(True|False)", line, re.IGNORECASE)
        if m:
            segmented = m.group(1).lower() == "true"

    if precision is None or recall is None:
        print("  WARNING: could not parse precision/recall from output")
        print(result.stdout[-3000:])
        return None, None, None

    print(f"  => precision={precision:.4f} ({precision*100:.2f}%)  recall={recall:.4f} ({recall*100:.2f}%)  segmented={segmented}")
    return precision, recall, segmented


def main():
    parser = argparse.ArgumentParser(description="Sweep k_max and plot precision/recall")
    parser.add_argument("checkpoint_path", type=str, help="Path to model checkpoint (.ckpt)")
    parser.add_argument("--knn-min", type=int, default=10, help="Minimum k_max (default: 10)")
    parser.add_argument("--knn-max", type=int, default=800, help="Maximum k_max (default: 400)")
    parser.add_argument("--step", type=int, default=25, help="Step size for k_max sweep (default: 25)")
    parser.add_argument("--r-max", type=float, default=1, help="Fixed r_max (default: 1)")
    parser.add_argument("--r-max-geometric", type=float, default=None, help="Fixed r_max_geometric in mm (default: None = disabled)")
    parser.add_argument("--dr-same-layer-cut", type=float, default=None, help="Remove edges between hits with |Δr| < this threshold (mm); 0 = disabled (default: 0)")
    parser.add_argument("--num-events", type=int, default=None, help="Number of test events to evaluate per sweep point (default: all in testset)")
    parser.add_argument("--output-dir", type=str, default=None, help="Directory to save plots (default: data/visuals/knn_sweep/)")
    parser.add_argument("--segmented", type=lambda x: x.lower() not in ('false', '0', 'no'), default=None, metavar='BOOL',
                        help="Override segmented flag passed to test script (default: read from graph_construction_latent.yaml)")
    args = parser.parse_args()

    checkpoint_path = Path(args.checkpoint_path)
    if not checkpoint_path.is_absolute():
        checkpoint_path = PIPELINE_ROOT / checkpoint_path

    knn_values = list(range(args.knn_min, args.knn_max + 1, args.step))

    ckpt_stem = checkpoint_path.stem
    base_dir = Path(args.output_dir) if args.output_dir else PIPELINE_ROOT / "data" / "visuals" / "knn_sweep"
    output_dir = base_dir / ckpt_stem
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Checkpoint: {checkpoint_path}")
    print(f"k_max sweep: {knn_values}")
    layer_str = f"{args.dr_same_layer_cut} mm" if args.dr_same_layer_cut is not None else "disabled"
    events_str = str(args.num_events) if args.num_events is not None else "all (from config)"
    print(f"r_max={args.r_max}, r_max_geometric={args.r_max_geometric} mm, dr_same_layer_cut={layer_str}, num_events={events_str}, segmented={args.segmented if args.segmented is not None else '(from config)'}")
    print(f"Output: {output_dir}")

    precisions = []
    recalls = []
    valid_knn = []
    segmented = None  # resolved from first successful run (either args.segmented or parsed from output)

    for knn in knn_values:
        p, r, seg = run_test(checkpoint_path, knn, args.r_max, args.r_max_geometric, dr_same_layer_cut=args.dr_same_layer_cut, num_events=args.num_events, segmented=args.segmented)
        if p is not None:
            precisions.append(p)
            recalls.append(r)
            valid_knn.append(knn)
            if segmented is None:
                segmented = seg

    if not valid_knn:
        print("\nERROR: no successful runs, cannot plot.")
        sys.exit(1)

    print(f"\n{'='*60}")
    print("SWEEP SUMMARY")
    print(f"{'='*60}")
    print(f"  segmented={segmented}  r_max={args.r_max}  r_max_geometric={args.r_max_geometric} mm")
    print(f"{'k_max':>8}  {'precision':>10}  {'recall':>10}")
    print(f"{'-'*34}")
    for k, p, r in zip(valid_knn, precisions, recalls):
        print(f"{k:>8}  {p:>10.4f}  {r:>10.4f}")

    # --- Save raw data ---
    raw_data_path = output_dir / "sweep_data.tsv"
    with open(raw_data_path, 'w') as f:
        f.write(f"# segmented={segmented}  r_max={args.r_max}  r_max_geometric={args.r_max_geometric} mm\n")
        f.write("k_max\tprecision\trecall\n")
        for k, p, r in zip(valid_knn, precisions, recalls):
            f.write(f"{k}\t{p:.6f}\t{r:.6f}\n")
    print(f"\nSaved raw data: {raw_data_path}")

    def padded_ylim(values, pad=0.05):
        lo, hi = min(values), max(values)
        span = hi - lo or 0.01
        return lo - pad * span, hi + pad * span

    # --- Precision plot ---
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(valid_knn, precisions, marker='o', linewidth=2, color='steelblue')
    ax.set_xlabel("k_max (KNN neighbors)", fontsize=12)
    ax.set_ylabel("Precision", fontsize=12)
    ax.set_title(f"Precision vs k_max\nr_max={args.r_max}, r_max_geometric={args.r_max_geometric} mm, segmented={segmented}", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(*padded_ylim(precisions))
    fig.tight_layout()
    precision_path = output_dir / f"precision_vs_knn.png"
    fig.savefig(precision_path, dpi=150)
    plt.close(fig)
    print(f"\nSaved: {precision_path}")

    # --- Recall plot ---
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(valid_knn, recalls, marker='o', linewidth=2, color='darkorange')
    ax.set_xlabel("k_max (KNN neighbors)", fontsize=12)
    ax.set_ylabel("Recall", fontsize=12)
    ax.set_title(f"Recall vs k_max\nr_max={args.r_max}, r_max_geometric={args.r_max_geometric} mm, segmented={segmented}", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(*padded_ylim(recalls))
    fig.tight_layout()
    recall_path = output_dir / f"recall_vs_knn.png"
    fig.savefig(recall_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {recall_path}")

    # --- Combined plot (dual y-axes since precision and recall live in different ranges) ---
    fig, ax1 = plt.subplots(figsize=(7, 4))
    ax1.plot(valid_knn, precisions, marker='o', linewidth=2, color='steelblue', label='Precision')
    ax1.set_xlabel("k_max (KNN neighbors)", fontsize=12)
    ax1.set_ylabel("Precision", fontsize=12, color='steelblue')
    ax1.tick_params(axis='y', labelcolor='steelblue')
    ax1.set_ylim(*padded_ylim(precisions))
    ax2 = ax1.twinx()
    ax2.plot(valid_knn, recalls, marker='s', linewidth=2, color='darkorange', label='Recall')
    ax2.set_ylabel("Recall", fontsize=12, color='darkorange')
    ax2.tick_params(axis='y', labelcolor='darkorange')
    ax2.set_ylim(*padded_ylim(recalls))
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=11)
    ax1.set_title(f"Precision & Recall vs k_max\nr_max={args.r_max}, r_max_geometric={args.r_max_geometric} mm, segmented={segmented}", fontsize=11)
    ax1.grid(True, alpha=0.3)
    fig.tight_layout()
    combined_path = output_dir / f"precision_recall_vs_knn.png"
    fig.savefig(combined_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {combined_path}")


if __name__ == "__main__":
    main()
