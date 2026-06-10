"""
Performance / wall-clock timing helpers for the inference benchmark.

Shared by performance_analysis/benchmark_inference.py (the harness) and
performance_analysis/plot_inference_timing.py (the plots). Kept here so the timing
logic — CUDA-synchronised wall clock, scratch run-dir setup that keeps the original
datasets pristine, spacepoint counting, and the (yellow-free) stage colour map —
lives in one reusable place.
"""

import os
import re
import shutil
import time
from pathlib import Path

# NOTE: torch is imported lazily inside the functions that need it (cuda_sync,
# count_spacepoints, aggregate_graph_times) so the torch-free helpers below
# (scratch run dirs, parse_multiplicity, count_csv_events) can be reused from the
# ACTS LCG environment, where importing torch into a C++ reco job is an ABI risk.

PIPELINE_ROOT = Path(__file__).resolve().parent.parent

# ── Stage identity (no yellow — mirrors track_evaluation_utils.SPECIES_COLORS family) ──
# Stage 4 is split into "track_building" (CC score cut + Wrangler segment extraction)
# and "segment_matching" (miniGNN embedding + cosine matching).
STAGE_COLORS = {
    "csv_to_pyg":       "#1f77b4",  # blue
    "build_graphs":     "#ff7f0e",  # orange
    "gnn_inference":    "#2ca02c",  # green
    "track_building":   "#d62728",  # red
    "segment_matching": "#9467bd",  # purple
    "evaluate":         "#8c564b",  # brown
}
STAGE_LABELS = {
    "csv_to_pyg":       "1. CSV-PyG (load)",
    "build_graphs":     "2. Graph constr.",
    "gnn_inference":    "3. GNN inference",
    "track_building":   "4. Track building (CC+Wrangler)",
    "segment_matching": "5. Segment matching (miniGNN)",
    "evaluate":         "6. Evaluation",
}

# Stacked-bar order: (per-event-dict key in the JSON record) -> (stage identity key)
STACK_STAGES = [
    ("csv_to_pyg",            "csv_to_pyg"),
    ("build_graphs",          "build_graphs"),
    ("gnn_inference",         "gnn_inference"),
    ("build_tracks_extract",  "track_building"),
    ("build_tracks_match",    "segment_matching"),
    ("evaluate",              "evaluate"),
]

# Stages counted toward the headline "inference total" (graph constr. + GNN + track build).
INFERENCE_STAGE_KEYS = ["build_graphs", "gnn_inference", "build_tracks_wall"]


# ── Timing ────────────────────────────────────────────────────────────────────
def cuda_sync(device):
    """Block until all CUDA work is done, but only when actually on CUDA."""
    import torch
    if str(device) == "cuda" and torch.cuda.is_available():
        torch.cuda.synchronize()


class StageTimer:
    """Context manager measuring wall time (perf_counter) with CUDA sync at both ends.

    The leading sync drains any pending async work so it is not charged to this stage;
    the trailing sync ensures this stage's GPU kernels finish before the clock stops.
    """

    def __init__(self, device="cpu"):
        self.device = str(device)
        self.elapsed = None

    def __enter__(self):
        cuda_sync(self.device)
        self._t0 = time.perf_counter()
        return self

    def __exit__(self, *exc):
        cuda_sync(self.device)
        self.elapsed = time.perf_counter() - self._t0
        return False  # never suppress exceptions


# ── Dataset bookkeeping ─────────────────────────────────────────────────────────
def parse_multiplicity(dataset_name):
    """`data_<N>` -> 4*N particles/event (electron+pion+muon+proton). None if no N."""
    m = re.search(r"data_(\d+)", str(dataset_name))
    return 4 * int(m.group(1)) if m else None


def count_csv_events(dataset_dir):
    """Number of simulated events available in the dataset's csv/ directory."""
    return len(sorted((Path(dataset_dir) / "csv").glob("event*-hits.csv")))


def _load_graphs(stage_dir, dataset, max_events=None):
    d = Path(stage_dir) / dataset
    paths = sorted(d.glob("*.pyg"))
    if not paths:  # fall back to whatever files are present
        paths = sorted(p for p in d.glob("*") if p.is_file())
    if max_events:
        paths = paths[:max_events]
    return paths


def count_spacepoints(stage_dir, dataset, max_events=None):
    """Average / min / max number of hits (spacepoints) per event, and event count.

    Reads graphs from a finished stage dir (e.g. graph_constructed_latent/<dataset>).
    """
    import torch
    counts = []
    for p in _load_graphs(stage_dir, dataset, max_events):
        g = torch.load(p, map_location="cpu", weights_only=False)
        if getattr(g, "hit_r", None) is not None:
            counts.append(int(g.hit_r.shape[0]))
        else:
            counts.append(int(g.num_nodes))
        del g
    if not counts:
        return 0.0, 0, 0, 0
    return sum(counts) / len(counts), min(counts), max(counts), len(counts)


def _mean_std(vals):
    if not vals:
        return 0.0, 0.0
    n = len(vals)
    mean = sum(vals) / n
    var = sum((v - mean) ** 2 for v in vals) / n if n > 1 else 0.0
    return mean, var ** 0.5


def aggregate_graph_times(stage_dir, dataset, attrs, max_events=None):
    """Aggregate per-event timing attributes from saved track-building graphs in one pass.

    `attrs` is a list of graph attribute names (e.g. ["time_taken", "time_extract",
    "time_match"]). Returns {attr: {"sum", "mean", "std", "n"}} with per-event mean/std
    (seconds), so callers get both the total and the event-to-event spread.
    """
    import torch
    series = {a: [] for a in attrs}
    for p in _load_graphs(stage_dir, dataset, max_events):
        g = torch.load(p, map_location="cpu", weights_only=False)
        for a in attrs:
            v = getattr(g, a, None)
            if v is not None:
                series[a].append(float(v))
        del g
    out = {}
    for a, vals in series.items():
        mean, std = _mean_std(vals)
        out[a] = {"sum": sum(vals), "mean": mean, "std": std, "n": len(vals)}
    return out


# ── Scratch run-dir (keeps the original datasets pristine) ──────────────────────
def _assert_under(path, root):
    """Guard: refuse to create/delete anything that is not inside `root`."""
    p = Path(path).resolve()
    r = Path(root).resolve()
    if not (p == r or r in p.parents):
        raise RuntimeError(f"Refusing to touch {p}: not under scratch root {r}")


def setup_scratch_run_dir(dataset_dir, scratch_root):
    """Create <scratch_root>/bench_<name> with csv/ symlinked from the real dataset.

    Pipeline stages run with debug_files=True, so nothing is ever deleted; even if a
    cleanup did fire it would only unlink the csv *symlink*, never the real csv/.
    The real dataset is only ever read, never written.
    """
    dataset_dir = Path(dataset_dir)
    scratch_root = Path(scratch_root)
    run_dir = scratch_root / f"bench_{dataset_dir.name}"
    _assert_under(run_dir, scratch_root)

    if run_dir.exists() or run_dir.is_symlink():
        teardown_scratch_run_dir(run_dir, scratch_root)
    run_dir.mkdir(parents=True)

    src_csv = (dataset_dir / "csv").resolve()
    if not src_csv.exists():
        raise FileNotFoundError(f"No csv/ directory in dataset {dataset_dir}")
    os.symlink(src_csv, run_dir / "csv")
    return run_dir


def teardown_scratch_run_dir(run_dir, scratch_root):
    """Remove a scratch run dir. Unlinks the csv symlink first so rmtree never follows
    it into the real dataset. Only ever operates under `scratch_root`."""
    run_dir = Path(run_dir)
    _assert_under(run_dir, scratch_root)
    csv_link = run_dir / "csv"
    if csv_link.is_symlink():
        csv_link.unlink()
    if run_dir.exists():
        shutil.rmtree(run_dir)
