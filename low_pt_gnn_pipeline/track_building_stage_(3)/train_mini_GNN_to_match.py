#!/usr/bin/env python3
"""
Train Mini-GNN Segment Embedder  [CPU-only]

Trains a small GNN to embed detector hit segments into a latent space where
segments from the same particle have high cosine similarity.

Uses a streaming DataLoader: only `batch_events` events are in memory at a
time — memory usage is constant regardless of dataset size.

This script is designed for CPU execution. The model is small (~100k params)
and the per-step forward pass is fast; GPU overhead (data transfer, launch
latency) is not worthwhile here. Keep batch_events small (e.g. 16-32).

Prerequisites:
    Run mini_gnn_data_mining.py to produce per-event .pyg files, then
    combine_mini_gnn_chunks.py to organise them into split directories.

Usage:
    python train_mini_GNN_to_match.py
    python train_mini_GNN_to_match.py --config path/to/config.yaml
    python train_mini_GNN_to_match.py --set hidden_dim=64 --set lr=0.0005
"""

import argparse
import math
import random
import sys
from pathlib import Path
from time import perf_counter

import numpy as np
import torch
import yaml
from torch.utils.data import DataLoader, Dataset
from torch_geometric.data import Batch
from tqdm import tqdm

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / "acorn"))
sys.path.insert(0, str(PIPELINE_ROOT))

from low_pt_custom_utils.mini_gnn_segment_embedding import (
    SegmentGNN,
    supcon_loss,
    NODE_FEATURE_DIM,
)


# ─── Dataset & Collation ─────────────────────────────────────────────────────


class SegmentEventDataset(Dataset):
    """
    Lazy-loading dataset over mined per-event segment files.

    Each .pyg file (produced by mini_gnn_data_mining.py) stores:
        {"seg_data_list": [PyG Data, ...], "particle_ids": [int, ...]}

    One item = one event. Files are loaded on demand, so only the current
    batch occupies memory.
    """

    def __init__(self, split_dir: Path, n_events: int = None):
        paths = sorted(split_dir.glob("*.pyg"))
        if not paths:
            raise FileNotFoundError(
                f"No per-event .pyg files found in {split_dir}\n"
                f"Run mini_gnn_data_mining.py first.\n"
                f"Do NOT run combine_mini_gnn_chunks.py — the streaming loader "
                f"reads individual files directly."
            )
        if n_events is not None:
            paths = paths[:n_events]
        self.paths = paths

    def __len__(self):
        return len(self.paths)

    def __getitem__(self, idx):
        cached = torch.load(self.paths[idx], map_location="cpu", weights_only=False)
        return cached["seg_data_list"], cached["particle_ids"]


def collate_events(batch):
    """
    Merge a list of (seg_data_list, particle_ids) tuples into one flat batch.

    Particle IDs are offset per event so that pid=7 in event A never collides
    with pid=7 in event B — they are separate particles.

    Returns:
        all_seg_data: flat list of all segment PyG Data objects
        pid_tensor:   (N,) long tensor of offset-adjusted particle IDs
        boundaries:   list of (start, end) index pairs, one per event
    """
    all_seg_data = []
    all_pids = []
    boundaries = []
    pid_offset = 0

    for seg_data_list, particle_ids in batch:
        start = len(all_seg_data)
        all_seg_data.extend(seg_data_list)
        for pid in particle_ids:
            all_pids.append(pid + pid_offset if pid > 0 else 0)
        if particle_ids:
            pid_offset += max(particle_ids) + 1
        boundaries.append((start, len(all_seg_data)))

    return all_seg_data, torch.tensor(all_pids, dtype=torch.long), boundaries


def make_loader(split_dir, n_events, batch_events, shuffle, num_workers):
    dataset = SegmentEventDataset(split_dir, n_events)
    return DataLoader(
        dataset,
        batch_size=batch_events,
        shuffle=shuffle,
        num_workers=num_workers,
        collate_fn=collate_events,
        pin_memory=False,
    ), len(dataset)


# ─── Training / Validation ───────────────────────────────────────────────────


def run_batches(model, batches, config, device):
    """
    Compute mean SupCon loss over an iterable of pre-collated batches (no grad).

    Args:
        batches: iterable of (all_seg_data, pid_tensor, boundaries)

    Returns:
        mean loss over all batches
    """
    temperature = config.get("temperature", 0.1)
    model.eval()
    losses = []

    with torch.no_grad():
        for all_seg_data, pid_tensor, _boundaries in batches:
            if len(all_seg_data) < 2:
                continue

            batch_data = Batch.from_data_list(all_seg_data).to(device)
            pid_tensor = pid_tensor.to(device)

            embeddings = model(batch_data.x, batch_data.edge_index, batch_data.batch)
            loss = supcon_loss(embeddings, pid_tensor, temperature=temperature)
            losses.append(loss.item())

    return float(np.mean(losses)) if losses else 0.0


def preload_val(val_loader):
    """
    Load all validation batches into RAM once.
    Stores the full collated tuples (all_seg_data, pid_tensor, boundaries)
    so they can be passed directly to both run_batches and compute_val_matching_accuracy.
    """
    return list(val_loader)


def compute_val_matching_accuracy(model, val_preloaded, config, device):
    """
    Greedy per-event segment matching accuracy on the pre-loaded validation batches.

    Vectorized: uses torch.triu_indices + argsort instead of Python list comprehension.
    """
    threshold = config.get("val_cos_sim_threshold", 0.5)
    n_correct = 0
    n_total = 0

    model.eval()
    with torch.no_grad():
        for all_seg_data, pid_tensor, boundaries in val_preloaded:
            if len(all_seg_data) < 2:
                continue

            batch_data = Batch.from_data_list(all_seg_data).to(device)
            embeddings = model(
                batch_data.x, batch_data.edge_index, batch_data.batch
            ).cpu()

            for start, end in boundaries:
                n = end - start
                if n < 2:
                    continue
                event_embs = embeddings[start:end]
                event_pids = pid_tensor[start:end]

                # Vectorized upper-triangle similarity + threshold + sort
                sim = event_embs @ event_embs.T
                rows, cols = torch.triu_indices(n, n, offset=1)
                sims_flat = sim[rows, cols]
                above = sims_flat >= threshold
                if not above.any():
                    continue
                order = torch.argsort(sims_flat[above], descending=True)
                r = rows[above][order].tolist()
                c = cols[above][order].tolist()

                matched = set()
                for i, j in zip(r, c):
                    if i not in matched and j not in matched:
                        matched |= {i, j}
                        if event_pids[i] == event_pids[j] and event_pids[i] > 0:
                            n_correct += 1
                        n_total += 1

    return n_correct / n_total if n_total > 0 else 0.0


# ─── Training Entry Point ────────────────────────────────────────────────────


def run_training(config):
    """Full training loop with streaming DataLoader."""

    model_save_dir = Path(config.get("model_save_dir", "saved_models"))
    if not model_save_dir.is_absolute():
        model_save_dir = PIPELINE_ROOT / model_save_dir
    model_save_dir.mkdir(parents=True, exist_ok=True)

    data_split = config.get("data_split", [500, 100, 0])
    n_train, n_val = data_split[0], data_split[1]

    gnn_config = config.get("gnn", {})
    hidden_dim = gnn_config.get("hidden_dim", 64)
    emb_dim = gnn_config.get("emb_dim", 32)
    n_layers = gnn_config.get("n_layers", 3)
    dropout = gnn_config.get("dropout", 0.0)
    proj_dim = gnn_config.get("proj_dim", None)
    proj_layers = gnn_config.get("proj_layers", 1)

    lr = config.get("lr", 1e-3)
    max_epochs = config.get("max_epochs", 100)
    patience = config.get("patience", 10)
    scheduler_patience = config.get("scheduler_patience", 5)
    scheduler_factor = config.get("scheduler_factor", 0.5)
    val_check_interval = config.get("val_check_interval", 1.0)
    batch_events = config.get("batch_events", 16)
    num_workers = config.get("num_workers", 0)

    random.seed(42)
    np.random.seed(42)
    torch.manual_seed(42)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    precomputed_dir = Path(config.get("precomputed_dir", "data/track_building/mini_gnn_segments"))
    if not precomputed_dir.is_absolute():
        precomputed_dir = PIPELINE_ROOT / precomputed_dir

    print(f"\nLoading mined segments from {precomputed_dir}")
    train_loader, n_train_events = make_loader(
        precomputed_dir / "trainset", n_train, batch_events, shuffle=True, num_workers=num_workers
    )
    val_loader, n_val_events = make_loader(
        precomputed_dir / "valset", n_val, batch_events, shuffle=False, num_workers=num_workers
    )
    print(f"  Pre-loading valset into RAM...", end="", flush=True)
    t0 = perf_counter()
    val_preloaded = preload_val(val_loader)
    print(f" done ({perf_counter() - t0:.1f}s)")

    n_train_steps = math.ceil(n_train_events / batch_events)
    validate_every = max(1, round(n_train_steps * val_check_interval))
    print(f"  trainset: {n_train_events} events → {n_train_steps} steps/epoch")
    print(f"  valset:   {n_val_events} events (pre-loaded)")
    print(f"  Validating every {validate_every} steps "
          f"({int(1 / val_check_interval)} times/epoch)")

    # ── Build model ──────────────────────────────────────────────────────────
    model = SegmentGNN(
        node_in_dim=NODE_FEATURE_DIM,
        hidden_dim=hidden_dim,
        emb_dim=emb_dim,
        n_layers=n_layers,
        dropout=dropout,
        proj_dim=proj_dim,
        proj_layers=proj_layers,
    ).to(device)

    n_params = sum(p.numel() for p in model.parameters())
    print(f"\nModel: hidden_dim={hidden_dim}, emb_dim={emb_dim}, n_layers={n_layers}, "
          f"proj_dim={proj_dim or hidden_dim}, proj_layers={proj_layers}, "
          f"dropout={dropout}, params={n_params:,}")
    print(f"Loss:  SupCon (τ={config.get('temperature', 0.1)})")
    print(f"Training: lr={lr}, max_epochs={max_epochs}, patience={patience}, "
          f"batch_events={batch_events}")
    print(f"Saving to: {model_save_dir}/")
    print()

    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", patience=scheduler_patience, factor=scheduler_factor, verbose=False
    )

    best_match_acc = -1.0
    patience_counter = 0
    model_save_path = None
    early_stopped = False

    print(f"{'Step':>8}  {'train_loss':>12}  {'val_loss':>12}  {'match_acc':>10}  {'lr':>10}  "
          f"{'t_load':>8}  {'t_collate':>9}  {'t_fwd':>8}  {'t_bwd':>8}  {'t_vloss':>8}  {'t_macc':>8}  {'status'}")
    print("-" * 140)

    for epoch in range(1, max_epochs + 1):
        model.train(True)
        temperature = config.get("temperature", 0.1)
        train_losses = []
        step_in_epoch = 0

        # Timing accumulators (seconds)
        t_load = t_collate = t_fwd = t_bwd = 0.0
        n_timed = 0

        pbar = tqdm(train_loader, desc=f"Epoch {epoch:>3} train", unit="batch", leave=False)

        t_iter_start = perf_counter()
        for all_seg_data, pid_tensor, boundaries in pbar:
            t_load += perf_counter() - t_iter_start  # time spent in DataLoader

            if len(all_seg_data) < 2:
                t_iter_start = perf_counter()
                continue

            t0 = perf_counter()
            batch_data = Batch.from_data_list(all_seg_data).to(device)
            pid_tensor = pid_tensor.to(device)
            t_collate += perf_counter() - t0

            t0 = perf_counter()
            optimizer.zero_grad()
            embeddings = model(batch_data.x, batch_data.edge_index, batch_data.batch)
            loss = supcon_loss(embeddings, pid_tensor, temperature=temperature)
            t_fwd += perf_counter() - t0

            if not loss.requires_grad:
                t_iter_start = perf_counter()
                continue

            t0 = perf_counter()
            loss.backward()
            optimizer.step()
            t_bwd += perf_counter() - t0

            train_losses.append(loss.item())
            step_in_epoch += 1
            n_timed += 1
            pbar.set_postfix(loss=f"{np.mean(train_losses):.4f}")

            # Intra-epoch validation
            if step_in_epoch % validate_every == 0 or step_in_epoch == n_train_steps:
                frac = step_in_epoch / n_train_steps
                step_label = f"{epoch - 1 + frac:.2f}" if frac < 1.0 else f"{epoch}"
                train_loss = float(np.mean(train_losses))
                train_losses = []

                t0 = perf_counter()
                val_loss = run_batches(model, val_preloaded, config, device)
                t_val_loss = perf_counter() - t0

                t0 = perf_counter()
                match_acc = compute_val_matching_accuracy(model, val_preloaded, config, device)
                t_match_acc = perf_counter() - t0

                scheduler.step(val_loss)
                current_lr = optimizer.param_groups[0]["lr"]

                n = max(1, n_timed)
                status = ""
                if match_acc > best_match_acc:
                    if model_save_path is not None and model_save_path.exists():
                        model_save_path.unlink()
                    best_match_acc = match_acc
                    model_save_path = (
                        model_save_dir
                        / f"mini_gnn_segment_embedder_match_acc={best_match_acc:.4f}.pt"
                    )
                    torch.save({
                        "state_dict": model.state_dict(),
                        "hidden_dim": hidden_dim,
                        "emb_dim": emb_dim,
                        "n_layers": n_layers,
                        "dropout": dropout,
                        "proj_dim": proj_dim,
                        "proj_layers": proj_layers,
                        "node_scales": gnn_config.get("node_scales", [1000.0, 1000.0, 500.0, 1000.0]),
                    }, model_save_path)
                    patience_counter = 0
                    status = f"* saved → {model_save_path.name}"
                else:
                    patience_counter += 1

                print(f"{step_label:>8}  {train_loss:>12.6f}  {val_loss:>12.6f}  "
                      f"{match_acc:>10.4f}  {current_lr:>10.2e}  "
                      f"{t_load/n:>7.3f}s  {t_collate/n:>8.3f}s  "
                      f"{t_fwd/n:>7.3f}s  {t_bwd/n:>7.3f}s  "
                      f"{t_val_loss:>7.1f}s  {t_match_acc:>7.1f}s  {status}")

                # Reset timing for next chunk
                t_load = t_collate = t_fwd = t_bwd = 0.0
                n_timed = 0

                if patience_counter >= patience:
                    print(f"\nEarly stopping: match_acc did not improve for {patience} checks.")
                    early_stopped = True
                    break

            t_iter_start = perf_counter()  # restart load timer for next batch

        if early_stopped:
            break

    print(f"\nTraining complete.")
    print(f"  Best match_acc: {best_match_acc:.4f}")
    print(f"  Model saved:    {model_save_path}")


# ─── Entry Point ─────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Train mini-GNN segment embedder with streaming DataLoader",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", type=str, default=None)
    parser.add_argument(
        "--set",
        action="append",
        metavar="KEY=VALUE",
        dest="overrides",
        default=[],
        help="Override a config key (e.g. --set hidden_dim=128)",
    )
    args = parser.parse_args()

    if args.config is None:
        config_file = (
            PIPELINE_ROOT
            / "acorn_configs"
            / "track_building_stage_(3)"
            / "mini_gnn_segment_matching_train.yaml"
        )
    else:
        config_file = Path(args.config)

    if not config_file.exists():
        raise FileNotFoundError(f"Config file not found: {config_file}")

    with open(config_file, "r") as f:
        raw_config = yaml.safe_load(f)

    # If this is a grid YAML (has fixed_params/grid_params), flatten fixed_params as the base config
    if "fixed_params" in raw_config:
        config = raw_config.get("fixed_params", {}).copy()
    else:
        config = raw_config

    # GNN arch params live under config["gnn"] — route them there so that
    # --set hidden_dim=128 (used by the grid search shell script) works correctly.
    _GNN_KEYS = {"hidden_dim", "emb_dim", "n_layers", "dropout", "proj_dim", "proj_layers"}
    for override in args.overrides:
        key, _, raw = override.partition("=")
        for cast in (int, float):
            try:
                raw = cast(raw)
                break
            except ValueError:
                pass
        if key in _GNN_KEYS:
            config.setdefault("gnn", {})[key] = raw
        else:
            config[key] = raw

    print("=" * 65)
    print("MINI-GNN SEGMENT EMBEDDER — TRAINING")
    print("=" * 65)
    print(f"Config: {config_file}")
    print()
    print(yaml.dump(config, default_flow_style=False, sort_keys=False).strip())
    print()

    run_training(config)


if __name__ == "__main__":
    main()
