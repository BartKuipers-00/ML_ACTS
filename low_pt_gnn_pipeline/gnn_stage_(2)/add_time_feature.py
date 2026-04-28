#!/usr/bin/env python3
"""
Add smeared time feature to graph_constructed_latent graphs for GNN training.

Adds new node attributes to each graph:
  hit_t_smeared : Gaussian-smeared truth time for pixel barrel hits (hit_region == 2,
                  i.e. volume_id 8). Zero for all other hits.
  hit_has_time  : Binary flag — 1.0 for pixel barrel hits, 0.0 for all others.
                  Omitted when --all is passed (every hit gets smeared time).

Usage:
    python add_time_feature.py --sigma 25
    python add_time_feature.py --sigma 25 --all   # smear all hits, no hit_has_time

Output directory: data/graph_constructed_latent_tsmear_{sigma}ps/
                  data/graph_constructed_latent_tsmear_{sigma}ps_all/  (with --all)

After running, update gnn_train.yaml:
    input_dir:     data/graph_constructed_latent_tsmear_{sigma}ps
    node_features: [hit_r, hit_phi, hit_z, hit_t_smeared, hit_has_time]
    node_scales:   [1000, 3.14, 500, 4000, 1]

    # with --all:
    node_features: [hit_r, hit_phi, hit_z, hit_t_smeared]
    node_scales:   [1000, 3.14, 500, 4000]
"""

import argparse
import torch
import numpy as np
from pathlib import Path
from tqdm import tqdm

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent

# hit_region == 2 corresponds to volume_id 8 (pixel barrel) — only these hits have timing
PIXEL_BARREL_REGION = 2


def add_time_feature(graph, sigma_ps: float, rng: np.random.Generator, all_hits: bool = False):
    """
    Add hit_t_smeared (and optionally hit_has_time) to a single PyG graph.

    When all_hits=False (default):
      Pixel barrel hits (hit_region == 2) receive Gaussian-smeared truth time.
      All other hits get hit_t_smeared = 0 and hit_has_time = 0.

    When all_hits=True:
      Every hit receives Gaussian-smeared truth time. hit_has_time is not added.
    """
    hit_t = graph.hit_t.numpy()

    if all_hits:
        noise = rng.normal(0.0, sigma_ps, size=len(hit_t))
        graph.hit_t_smeared = torch.tensor(hit_t + noise, dtype=torch.float32)
    else:
        hit_region = graph.hit_region.numpy()
        pixel_mask = hit_region == PIXEL_BARREL_REGION

        hit_has_time = torch.zeros(len(hit_region), dtype=torch.float32)
        hit_has_time[pixel_mask] = 1.0

        hit_t_smeared = np.zeros(len(hit_region), dtype=np.float32)
        noise = rng.normal(0.0, sigma_ps, size=int(pixel_mask.sum()))
        hit_t_smeared[pixel_mask] = hit_t[pixel_mask] + noise

        graph.hit_t_smeared = torch.tensor(hit_t_smeared, dtype=torch.float32)
        graph.hit_has_time = hit_has_time

    return graph


def main():
    parser = argparse.ArgumentParser(
        description="Add Gaussian-smeared time feature to GNN input graphs"
    )
    parser.add_argument(
        "--sigma", type=float, required=True,
        help="Time smearing sigma in picoseconds, e.g. --sigma 25"
    )
    parser.add_argument(
        "--input-dir", type=str, default=None,
        help="Input graph directory (default: data/graph_constructed_latent)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility (default: 42)"
    )
    parser.add_argument(
        "--all", dest="all_hits", action="store_true",
        help="Smear all hits (not just pixel barrel); omits hit_has_time feature"
    )
    args = parser.parse_args()

    if args.input_dir is None:
        input_base = PIPELINE_ROOT / "data" / "graph_constructed_latent"
    else:
        p = Path(args.input_dir)
        input_base = p if p.is_absolute() else PIPELINE_ROOT / p

    if not input_base.exists():
        print(f"ERROR: input directory not found: {input_base}")
        raise SystemExit(1)

    sigma_ps = args.sigma
    all_hits = args.all_hits
    suffix = f"_all" if all_hits else ""
    output_base = PIPELINE_ROOT / "data" / f"graph_constructed_latent_tsmear_{int(sigma_ps)}ps{suffix}"
    rng = np.random.default_rng(args.seed)

    print("=" * 70)
    print("ADD SMEARED TIME FEATURE")
    print("=" * 70)
    if all_hits:
        print("Mode:         all hits (no hit_has_time)")
    else:
        print(f"Pixel barrel: hit_region == {PIXEL_BARREL_REGION}  (volume_id 8)")
    print(f"Sigma:        {sigma_ps} ps")
    print(f"Input:        {input_base}")
    print(f"Output:       {output_base}")
    print()

    for split in ["trainset", "valset", "testset"]:
        in_dir = input_base / split
        out_dir = output_base / split

        if not in_dir.exists():
            print(f"Skipping {split} (not found)")
            continue

        out_dir.mkdir(parents=True, exist_ok=True)
        files = sorted(in_dir.glob("*.pyg"))

        for f in tqdm(files, desc=split, unit="graph"):
            graph = torch.load(f, map_location="cpu", weights_only=False)
            graph = add_time_feature(graph, sigma_ps, rng, all_hits=all_hits)
            torch.save(graph, out_dir / f.name)

    out_dir_name = f"graph_constructed_latent_tsmear_{int(sigma_ps)}ps{suffix}"
    print("\nDone. Update gnn_train.yaml:")
    print(f"  input_dir:     data/{out_dir_name}")
    if all_hits:
        print(f"  node_features: [hit_r, hit_phi, hit_z, hit_t_smeared]")
        print(f"  node_scales:   [1000, 3.14, 500, 4000]")
    else:
        print(f"  node_features: [hit_r, hit_phi, hit_z, hit_t_smeared, hit_has_time]")
        print(f"  node_scales:   [1000, 3.14, 500, 4000, 1]")


if __name__ == "__main__":
    main()
