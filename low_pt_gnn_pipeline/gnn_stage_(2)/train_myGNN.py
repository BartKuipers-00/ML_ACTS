
import argparse
import os
import sys
import yaml
from pathlib import Path
import torch
from pytorch_lightning import Trainer
from pytorch_lightning.callbacks import ModelCheckpoint, EarlyStopping
from pytorch_lightning.loggers import CSVLogger, WandbLogger


SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / 'acorn'))
sys.path.insert(0, str(PIPELINE_ROOT))

from acorn.core.core_utils import get_stage_module
from acorn.utils.loading_utils import add_variable_name_prefix_in_config
from acorn.stages.edge_classifier.models.interaction_gnn import InteractionGNN2
from acorn.stages.edge_classifier.edge_classifier_stage import GraphDataset
from torch_geometric.loader import DataLoader
from torch_geometric.data import Dataset
from class_resolver import ClassResolver


class BoolEdgeYGraphDataset(GraphDataset):
    """GraphDataset that casts edge_y to bool before preprocessing.

    Our graph builder stores edge_y as int64, but ACORN's handle_weighting builds
    the per-edge weight mask via torch.ones_like(event.edge_y). With an int64
    edge_y that mask is an integer tensor, so `weights[mask] = value` becomes
    index-assignment (only touching rows 0/1) and every weighting spec silently
    no-ops -- zeroing the weight of ALL true edges (pos_loss collapses to 0).
    Casting edge_y to bool makes the weighting block behave as intended. This
    fixes existing graphs at load time, so no re-migration is needed.
    """

    def preprocess_event(self, event):
        if getattr(event, "edge_y", None) is not None and event.edge_y.dtype != torch.bool:
            event.edge_y = event.edge_y.bool()
        return super().preprocess_event(event)


class WandbInteractionGNN(InteractionGNN2):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Use the bool-edge_y dataset so the weighting block is applied correctly.
        # gnn_train.yaml does not set dataset_class, so this replaces the default.
        self.dataset_resolver = ClassResolver(
            [BoolEdgeYGraphDataset], base=Dataset, default=BoolEdgeYGraphDataset
        )

    def train_dataloader(self):
        """Override to enable shuffling for better gradient accumulation."""
        if self.trainset is None:
            return None
        num_workers = self.hparams.get("num_workers", [1, 1, 1])[0]
        return DataLoader(self.trainset, batch_size=1, num_workers=num_workers, shuffle=True)

    def training_step(self, batch, batch_idx):
        """Override to enable batch-level logging (on_step=True) for W&B."""
        output = self(batch)
        loss, pos_loss, neg_loss = self.loss_function(output, batch)

        # Scale loss for gradient accumulation (maintains same effective LR)
        accum = self.trainer.accumulate_grad_batches
        scaled_loss = loss / accum

        # Log unscaled loss so values are comparable across different accumulation settings
        self.log("train_loss", loss, on_step=True, on_epoch=True, batch_size=1, sync_dist=True)
        self.log("train_pos_loss", pos_loss, on_step=True, on_epoch=True, batch_size=1, sync_dist=True)
        self.log("train_neg_loss", neg_loss, on_step=True, on_epoch=True, batch_size=1, sync_dist=True)

        return scaled_loss

def main():
    
    # Parse arguments
    parser = argparse.ArgumentParser(description='Train GNN edge classifier')
    parser.add_argument(
        '--config',
        type=str,
        default=None,
        help='Path to config file (default: acorn_configs/gnn_stage_(2)/gnn_train.yaml)'
    )
    args = parser.parse_args()
    
    # Enable Tensor Cores for faster matrix operations on L40S GPU if available.    (GPU optimization)
    if torch.cuda.is_available():
        torch.set_float32_matmul_precision('high')  # or 'high' for even more speed

    SCRIPT_DIR = Path(__file__).resolve().parent
    
    # Use provided config or default
    if args.config:
        config_file = Path(args.config)
        if not config_file.is_absolute():
            # If path starts with acorn_configs, make it relative to pipeline root
            if str(config_file).startswith('acorn_configs'):
                config_file = PIPELINE_ROOT / config_file
            else:
                config_file = SCRIPT_DIR / config_file
    else:
        config_file = PIPELINE_ROOT / 'acorn_configs' / 'gnn_stage_(2)' / 'gnn_train.yaml'
    
    if not config_file.exists():
        raise FileNotFoundError(f"Config file not found: {config_file}")
    
    print(f"Loading config from: {config_file}\n")
    
    # Load config
    with open(config_file, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    # Make paths absolute relative to pipeline root directory
    for path_key in ['input_dir', 'stage_dir']:
        if path_key in config and not Path(config[path_key]).is_absolute():
            config[path_key] = str(PIPELINE_ROOT / config[path_key])

    if not config.get("variable_with_prefix"):
        config = add_variable_name_prefix_in_config(config)
    
    print("="*70)
    print("TRAINING CONFIGURATION")
    print("="*70)
    print(yaml.dump(config))
    print("="*70)
    

    stage_module_class = WandbInteractionGNN
    
    # Setup stage directory
    os.makedirs(config["stage_dir"], exist_ok=True)
    
    # Get stage module
    stage_module, ckpt_config, default_root_dir, checkpoint = get_stage_module(
        config, stage_module_class, checkpoint_path=None, checkpoint_resume_dir=None
    )
    
    if (not config.get("variable_with_prefix")) or config.get(
        "add_variable_name_prefix_in_ckpt"
    ):
        stage_module._hparams = add_variable_name_prefix_in_config(
            stage_module._hparams
        )
    
    # Create checkpoint callback to save best model
    checkpoint_callback = ModelCheckpoint(
        dirpath=Path(config["stage_dir"]) / "checkpoints",
        filename="gnn_best_val_loss_{val_loss:.4f}",
        monitor="val_loss",
        mode="min",
        save_top_k=1,  # Only save the single best model
        save_last=False,  # Don't save last.ckpt
        verbose=True,
    )

    # Create early stopping callback
    early_stopping = EarlyStopping(
        monitor="val_loss",
        patience=config.get("early_stopping_patience", 2),
        mode="min",
        verbose=True,
        strict=True,  # Crash if monitored metric is not found
    )

    # Setup loggers
    loggers = [CSVLogger(save_dir=config["stage_dir"], name="logs")]
    
    # Add W&B logger if enabled
    if config.get("use_wandb", False):
        loggers.append(WandbLogger(
            project=config.get("project", "GNN_Training"),
            entity=config.get("wandb_entity"),
            config=config,
        ))
    
    trainer = Trainer(
        accelerator=config.get("accelerator", "cpu"),
        devices=config.get("devices", 1),
        num_nodes=config.get("nodes", 1),
        max_epochs=config["max_epochs"],
        callbacks=[checkpoint_callback, early_stopping],
        logger=loggers,
        log_every_n_steps=config.get("wandb_log_every_n_batches", 1),  # Control logging frequency
        check_val_every_n_epoch=config.get("check_val_every_n_epoch", 1),
        val_check_interval=config.get("val_check_interval"),  # Check validation within epochs (e.g., 0.5 = every 50%)
        accumulate_grad_batches=config.get("accumulate_grad_batches", 1),
        precision=config.get("precision", 32),
        enable_progress_bar=True,
        enable_model_summary=True,
    )
    
    print("\n" + "="*70)
    print("STARTING TRAINING")
    print("="*70 + "\n")
    
    trainer.fit(stage_module)
    
    print("\n" + "="*70)
    print("TRAINING COMPLETE!")
    print("="*70)
    print(f"Best model: {checkpoint_callback.best_model_path}")
    print(f"Best val_loss: {checkpoint_callback.best_model_score:.4f}")
    print(f"Checkpoints saved in: {config['stage_dir']}/checkpoints/")
    print(f"Logs saved in: {config['stage_dir']}/logs/")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()


