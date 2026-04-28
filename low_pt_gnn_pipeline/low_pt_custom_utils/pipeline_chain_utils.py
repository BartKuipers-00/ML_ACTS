"""
Pipeline chain utilities for the full GNN tracking chain.

Each function runs one stage with a configurable data directory,
wrapping the individual stage scripts' key functions so that
full_chain_gnn.py can stay minimal.
"""

import shutil
import sys
from pathlib import Path
import yaml
import torch

PIPELINE_ROOT = Path(__file__).resolve().parent.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / "acorn"))
sys.path.insert(0, str(PIPELINE_ROOT))
sys.path.insert(0, str(PIPELINE_ROOT / "latent_stage_(1)"))
sys.path.insert(0, str(PIPELINE_ROOT / "track_building_stage_(3)"))

from acts_custom_low_pt_reader import ActsCustomLowPTReader
from build_latent_graphs_fast import load_model as load_latent_model, run_graph_construction
from acorn.core.infer_stage import infer as acorn_infer
from acorn.utils.loading_utils import add_variable_name_prefix_in_config
from low_pt_custom_utils.mini_gnn_segment_embedding import load_segment_gnn
from GNN_segment_matching_track_builder import run_gnn_segment_matching
from low_pt_custom_utils.track_evaluation_utils import (
    run_evaluation,
    save_evaluation_results,
    run_plotting,
)


def _load_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f)


def _cleanup(directory):
    """Delete a directory tree and log it."""
    d = Path(directory)
    if d.exists():
        shutil.rmtree(d)
        print(f"  [cleanup] removed {d}")


# ── Stage functions ──────────────────────────────────────────────────────────


def stage_csv_to_pyg(run_dir, num_events, dataset):
    """Stage 1: ACTS CSV → PyTorch Geometric feature store."""
    cfg = _load_yaml(
        PIPELINE_ROOT / "acorn_configs" / "latent_stage_(1)" / "convert_csv_to_pyg_sets.yaml"
    )
    cfg["input_dir"]     = str(run_dir / "csv")
    cfg["stage_dir"]     = str(run_dir / "feature_store")
    cfg["detector_path"] = str(run_dir / "csv" / "detectors.csv")
    cfg["data_split"]    = [0, 0, num_events]
    cfg["input_sets"]    = [dataset]
    ActsCustomLowPTReader.infer(cfg)
    return run_dir / "feature_store"


def _resolve_device(device_str):
    """Resolve 'auto' to 'cuda' if a GPU is available, else 'cpu'."""
    if device_str in (None, "auto"):
        return "cuda" if torch.cuda.is_available() else "cpu"
    return device_str


def stage_build_graphs(run_dir, latent_ckpt, graph_config, dataset, debug_files=False):
    """Stage 2: KNN graph construction in the learned latent embedding space."""
    model, hparams = load_latent_model(str(latent_ckpt))
    cfg = {
        "input_dir":  str(run_dir / "feature_store"),
        "output_dir": str(run_dir / "graph_constructed_latent"),
        "datasets":   [dataset],
        **graph_config,
        "device": _resolve_device(graph_config.get("device")),
    }
    run_graph_construction(model, hparams, cfg)
    if not debug_files:
        _cleanup(run_dir / "csv")
        _cleanup(run_dir / "feature_store")
    return run_dir / "graph_constructed_latent"


def stage_gnn_inference(run_dir, gnn_ckpt, num_events, dataset, extra_config=None, debug_files=False):
    """Stage 3: GNN InteractionGNN edge classification inference."""
    cfg = _load_yaml(
        PIPELINE_ROOT / "acorn_configs" / "gnn_stage_(2)" / "gnn_infer.yaml"
    )
    cfg["input_dir"]     = str(run_dir / "graph_constructed_latent")
    cfg["stage_dir"]     = str(run_dir / "gnn_stage")
    cfg["data_split"]    = [0, 0, num_events]
    cfg["checkpoint"]    = str(gnn_ckpt)
    cfg["skip_existing"] = False
    if extra_config:
        cfg.update(extra_config)
    if not cfg.get("variable_with_prefix"):
        cfg = add_variable_name_prefix_in_config(cfg)

    tmp_cfg = run_dir / "_tmp_gnn_infer.yaml"
    with open(tmp_cfg, "w") as f:
        yaml.dump(cfg, f)

    acorn_infer(str(tmp_cfg), verbose=False, checkpoint=str(gnn_ckpt))
    tmp_cfg.unlink(missing_ok=True)
    if not debug_files:
        _cleanup(run_dir / "graph_constructed_latent")
    return run_dir / "gnn_stage"


def stage_build_tracks(run_dir, mini_gnn_ckpt, num_events, dataset, track_config, device=None, debug_files=False):
    """Stage 4: GNN segment matching track building."""
    device = _resolve_device(device)
    model = load_segment_gnn(str(mini_gnn_ckpt), device=device)

    cfg = {
        "input_dir":  str(run_dir / "gnn_stage"),
        "stage_dir":  str(run_dir / "track_building"),
        "data_split": [0, 0, num_events],
        **track_config,
    }
    cfg["gnn"] = {**cfg.get("gnn", {}), "model_path": str(mini_gnn_ckpt)}

    run_gnn_segment_matching(dataset, cfg, model, device)
    if not debug_files:
        _cleanup(run_dir / "gnn_stage")
    return run_dir / "track_building"


def stage_evaluate(run_dir, dataset, eval_config, debug_files=False):
    """Stage 5: Evaluate reconstructed tracks and generate efficiency plots."""
    cfg = {
        **eval_config,
        "input_dir": str(run_dir / "track_building"),
    }
    eval_output_dir = run_dir / "track_evaluation" / dataset
    plot_output_dir = run_dir / "visuals" / "track_metrics" / dataset

    evaluated_events, summary, summary_text = run_evaluation(dataset, cfg)
    save_evaluation_results(evaluated_events, summary, summary_text, dataset, eval_output_dir, debug_files=debug_files)
    run_plotting(evaluated_events, summary, dataset, plot_output_dir, eval_config["plots"])
    if not debug_files:
        _cleanup(run_dir / "track_building")

    print(summary_text)
    return summary
