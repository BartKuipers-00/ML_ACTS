#!/usr/bin/env python
"""
CPU smoke test for the GNN edge classifier on the migrated graph structure.

Reproduces the exact path that crashed the express GPU job — data loading +
construct_weighting (the track_particle_primary -> edge mapping via
track_to_edge_map) — then runs one forward pass + loss, all on CPU. No GPU,
no W&B, no trainer. Uses only trainset (data_split [2,0,0]) so it touches only
already-migrated graphs.

Run:
    cd gnn_stage_(2) && python test_cpu_forward_pass.py
Exit code 0 = the fix works end-to-end; non-zero = still broken.
"""
import sys
import yaml
import traceback
from pathlib import Path
import torch

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / "acorn"))
sys.path.insert(0, str(PIPELINE_ROOT))

from acorn.core.core_utils import get_stage_module
from acorn.utils.loading_utils import add_variable_name_prefix_in_config
from train_myGNN import WandbInteractionGNN


def main():
    config_file = PIPELINE_ROOT / "acorn_configs" / "gnn_stage_(2)" / "gnn_train.yaml"
    with open(config_file) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    for path_key in ["input_dir", "stage_dir"]:
        if path_key in config and not Path(config[path_key]).is_absolute():
            config[path_key] = str(PIPELINE_ROOT / config[path_key])

    # Force a cheap CPU-only, single-graph run on trainset alone.
    config["data_split"] = [2, 0, 0]
    config["accelerator"] = "cpu"
    config["devices"] = 1
    config["nodes"] = 1
    config["num_workers"] = [0, 0, 0]
    config["use_wandb"] = False
    config["logger"] = False

    if not config.get("variable_with_prefix"):
        config = add_variable_name_prefix_in_config(config)

    stage_module, *_ = get_stage_module(
        config, WandbInteractionGNN, checkpoint_path=None, checkpoint_resume_dir=None
    )

    print("=" * 60)
    print("[1/3] setup('fit') — data load + run_data_tests + weighting")
    print("=" * 60)
    stage_module.setup("fit")            # <-- this is where the express job crashed
    print("    OK: data loaded and weighting constructed without error")

    print("\n[2/3] pulling one batch from train_dataloader")
    batch = next(iter(stage_module.train_dataloader()))
    print(f"    nodes={batch.num_nodes}  edges={batch.edge_index.shape[1]}  "
          f"track_to_edge_map={tuple(batch.track_to_edge_map.shape)}  "
          f"edge_weights={tuple(batch.edge_weights.shape) if hasattr(batch,'edge_weights') else 'NA'}")

    print("\n[3/3] one forward pass + loss (no grad, CPU)")
    stage_module.eval()
    with torch.no_grad():
        output = stage_module(batch)
        loss, pos_loss, neg_loss = stage_module.loss_function(output, batch)
    print(f"    output edge-scores shape={tuple(output.shape)}")
    print(f"    loss={loss.item():.4f}  pos={pos_loss.item():.4f}  neg={neg_loss.item():.4f}")

    print("\n" + "=" * 60)
    print("SUCCESS — GNN trains-path works on the migrated data structure")
    print("=" * 60)


if __name__ == "__main__":
    try:
        main()
    except Exception:
        print("\n" + "=" * 60)
        print("FAILED — see traceback below")
        print("=" * 60)
        traceback.print_exc()
        sys.exit(1)
