#!/usr/bin/env python3
"""
Train metric learning model on PyG feature store

Loads PyG graphs with spacepoints (r, φ, z) and particle IDs,
trains embedding model using contrastive learning.

Usage:
    python train_latent_cluster_learning.py
    
Configuration:
    acorn_configs/latent_cluster_learning_train.yaml
"""

import sys
import argparse
from pathlib import Path
import yaml
import torch
from torch.utils.data import DataLoader
from pytorch_lightning import Trainer
from pytorch_lightning.callbacks import ModelCheckpoint, EarlyStopping, Callback
from pytorch_lightning.loggers import CSVLogger, WandbLogger

# Enable Tensor Cores for L40S GPU
torch.set_float32_matmul_precision('medium')

# Add acorn to path
SCRIPT_DIR = Path(__file__).resolve().parent
WORKSPACE_ROOT = SCRIPT_DIR.parent
PIPELINE_ROOT = SCRIPT_DIR.parent  # Alias for consistency with other scripts
sys.path.insert(0, str(WORKSPACE_ROOT / 'acorn'))

from acorn.stages.graph_construction.models.metric_learning import MetricLearning
import torch.optim.lr_scheduler as lr_scheduler

# Override acorn's build_edges with proper KNN+radius (no FRNN needed, no PyG bias)
sys.path.insert(0, str(PIPELINE_ROOT))
from low_pt_custom_utils.graph_utils import patch_acorn_build_edges
patch_acorn_build_edges()


class MetricLearningWithReduceLROnPlateau(MetricLearning):
    """ use ReduceLROnPlateau and ensure validation metrics are prominently logged to W&B"""

    def _build_particle_only_edges(self, batch):
        """Build consecutive track edges per particle, ignoring segment boundaries.

        Used when segmented=False in config: all consecutive hits of the same particle
        (sorted by hit_t) are true edges, including across segment boundaries.
        """
        particle_ids = batch.hit_particle_id
        hit_t = batch.hit_t
        edges = []
        unique_particles = particle_ids.unique()
        unique_particles = unique_particles[unique_particles != 0]  # exclude noise
        for pid in unique_particles:
            indices = (particle_ids == pid).nonzero(as_tuple=True)[0]
            order = hit_t[indices].argsort()
            sorted_indices = indices[order]
            if len(sorted_indices) >= 2:
                edges.append(torch.stack([sorted_indices[:-1], sorted_indices[1:]], dim=0))
        if not edges:
            return torch.zeros(2, 0, dtype=torch.long, device=particle_ids.device)
        return torch.cat(edges, dim=1)

    def get_truth(self, batch, pred_edges, undirected=True):
        """Override to recompute track_edges when segmented=False. Always undirected."""
        if not self.hparams.get("segmented", True):
            batch.track_edges = self._build_particle_only_edges(batch)
        return super().get_truth(batch, pred_edges, undirected=True)

    def train_dataloader(self):
        """Override to add shuffle=True for better training"""
        if self.trainset is None:
            return None
        num_workers = (
            16 if (  "num_workers" not in self.hparams or self.hparams["num_workers"] is None )
            else self.hparams["num_workers"][0]
        )
        # With batch_size=1, collate_fn just returns the first (and only) element
        return DataLoader(
            self.trainset, 
            batch_size=1, 
            num_workers=num_workers, 
            shuffle=True,
            collate_fn=lambda lst: lst[0]  # Required for PyG Data objects with batch_size=1
        )
    
    def validation_step(self, batch, batch_idx):
        """Override to ensure metrics are logged with prog_bar=True for W&B visibility"""
        knn_val = self.hparams["knn_val"]
        # Call parent's shared_evaluation which computes and logs all metrics
        outputs = self.shared_evaluation(batch, self.hparams["r_train"], knn_val)
        return outputs["loss"]
    
    def log_metrics(self, batch, loss, pred_edges, true_edges, truth, weights):
        """Override parent's log_metrics to ensure all metrics show up in W&B with prog_bar"""
        from acorn.stages.graph_construction.utils import build_signal_edges
        
        signal_true_edges = build_signal_edges(
            batch, self.hparams["weighting"], true_edges
        )
        true_pred_edges = pred_edges[:, truth == 1]
        signal_true_pred_edges = pred_edges[:, (truth == 1) & (weights > 0)]

        total_eff = true_pred_edges.shape[1] / true_edges.shape[1]
        signal_eff = signal_true_pred_edges.shape[1] / signal_true_edges.shape[1]
        total_pur = true_pred_edges.shape[1] / pred_edges.shape[1]
        signal_pur = signal_true_pred_edges.shape[1] / pred_edges.shape[1]
        f1 = 2 * (signal_eff * signal_pur) / (signal_eff + signal_pur)

        current_lr = self.optimizers().param_groups[0]["lr"]
        
        # Log with prog_bar=True to ensure W&B picks them up prominently
        self.log_dict(
            {
                "val_loss": loss,
                "lr": current_lr,
                "total_eff": total_eff,
                "total_pur": total_pur,

                "f1": f1,
            },
            batch_size=1,
            on_epoch=True,
            on_step=True,
            prog_bar=True,  # This ensures metrics show in progress bar AND are prominent in W&B
            sync_dist=True,
        )
    
    def configure_optimizers(self):
        optimizer = [
            torch.optim.AdamW(
                self.parameters(),
                lr=(self.hparams["lr"]),
                betas=(0.9, 0.999),
                eps=1e-08,
                amsgrad=True,
            )
        ]
        
        # Use ReduceLROnPlateau (decreases LR only when validation metric doesn't improve)
        scheduler_type = self.hparams.get("scheduler", "ReduceLROnPlateau")
        
        if scheduler_type == "ReduceLROnPlateau":
            metric_mode = self.hparams.get("metric_mode", "min")
            metric_to_monitor = self.hparams.get("metric_to_monitor", "val_loss")
            scheduler = [
                {
                    "scheduler": lr_scheduler.ReduceLROnPlateau(
                        optimizer[0],
                        mode=metric_mode,
                        factor=self.hparams["factor"],
                        patience=self.hparams["patience"],
                        verbose=True,
                    ),
                    "interval": "epoch",
                    "frequency": 1,
                    "monitor": metric_to_monitor,
                }
            ]
        elif scheduler_type == "StepLR":
            scheduler = [
                {
                    "scheduler": lr_scheduler.StepLR(
                        optimizer[0],
                        step_size=self.hparams["patience"],
                        gamma=self.hparams["factor"],
                    ),
                    "interval": "epoch",
                    "frequency": 1,
                }
            ]
        else:
            raise ValueError(f"Unknown scheduler: {scheduler_type}")
        
        return optimizer, scheduler


class LossPrinterCallback(Callback):
    """Print clean one-line summary per epoch"""
    
    def __init__(self):
        self.train_loss = None
    
    def on_train_epoch_end(self, trainer, pl_module):
        """Store training loss"""
        metrics = trainer.callback_metrics
        if 'train_loss' in metrics:
            self.train_loss = metrics['train_loss'].item()
    
    def on_validation_epoch_end(self, trainer, pl_module):
        """Print epoch summary"""
        metrics = trainer.callback_metrics
        
        # Extract metrics
        epoch = trainer.current_epoch
        train_loss = self.train_loss if self.train_loss is not None else float('nan')
        val_loss = metrics.get('val_loss', torch.tensor(float('nan'))).item()
        f1 = metrics.get('f1', torch.tensor(float('nan'))).item()
        signal_eff = metrics.get('signal_eff', torch.tensor(float('nan'))).item()
        signal_pur = metrics.get('signal_pur', torch.tensor(float('nan'))).item()
        
        # Print clean one-liner
        print(f"Epoch {epoch}: train_loss={train_loss:.4f}, val_loss={val_loss:.4f}, "
              f"f1={f1:.4f}, eff={signal_eff:.4f}, pur={signal_pur:.4f}")


def main():
    """Train metric learning model"""
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Train latent space metric learning model')
    parser.add_argument('--config', type=str, default=None,
                        help='Path to config file (default: acorn_configs/latent_cluster_learning_train.yaml)')
    args = parser.parse_args()
    
    # Load configuration
    if args.config:
        config_path = Path(args.config)
        if not config_path.is_absolute():
            config_path = WORKSPACE_ROOT / config_path
    else:
        config_path = WORKSPACE_ROOT / 'acorn_configs' / 'latent_stage_(1)' / 'latent_cluster_learning_train.yaml'
    
    print("="*80)
    print("METRIC LEARNING TRAINING - LATENT SPACE EMBEDDING")
    print("="*80)
    print()
    print(f"Loading config from: {config_path.relative_to(WORKSPACE_ROOT) if config_path.is_relative_to(WORKSPACE_ROOT) else config_path.name}")
    print()
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Convert relative paths to absolute paths (relative to pipeline root)
    if 'input_dir' in config and not Path(config['input_dir']).is_absolute():
        config['input_dir'] = str(PIPELINE_ROOT / config['input_dir'])
    if 'stage_dir' in config and not Path(config['stage_dir']).is_absolute():
        config['stage_dir'] = str(PIPELINE_ROOT / config['stage_dir'])
    
    print("="*80)
    print("TRAINING CONFIGURATION")
    print("="*80)
    print(yaml.dump(config, default_flow_style=False, sort_keys=False))
    print("="*80)
    print()
    
    # Create model (use subclass that supports ReduceLROnPlateau)
    stage_module = MetricLearningWithReduceLROnPlateau(config)

    # Setup output directory (config['stage_dir'] is now absolute)
    output_dir = Path(config['stage_dir'])
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get accelerator setting
    accelerator = config.get('accelerator', 'gpu')
    if accelerator == 'gpu' and not torch.cuda.is_available():
        print("WARNING: GPU requested but not available, falling back to CPU")
        accelerator = 'cpu'
    
    # Setup callbacks
    checkpoint_callback = ModelCheckpoint(
        dirpath=output_dir / 'checkpoints',
        filename=config.get('checkpoint_filename', 'metric-{epoch:02d}-{val_loss:.4f}'),
        monitor=config.get('metric_to_monitor', 'val_loss'),
        mode=config.get('metric_mode', 'min'),
        save_top_k=config.get('save_top_k', 1),
        save_last=config.get('save_last', False),
    )
    
    early_stop_callback = EarlyStopping(
        monitor=config.get('metric_to_monitor', 'val_loss'),
        patience=config.get('early_stopping_patience', 10),
        mode=config.get('metric_mode', 'min'),
        verbose=True,
    )
    
    loss_printer = LossPrinterCallback()
    
    # Setup loggers
    loggers = [CSVLogger(save_dir=config["stage_dir"], name="logs")]
    
    # Add W&B logger if enabled
    if config.get("use_wandb", False):
        wandb_kwargs = {
            "project": config.get("project", "Low_pt_latent_map_MLP"),
            "entity": config.get("wandb_entity"),
            "config": config,
        }
        # Add run name if specified (for grid search)
        if "run_name" in config:
            wandb_kwargs["name"] = config["run_name"]
        
        loggers.append(WandbLogger(**wandb_kwargs))
    
    print("="*80)
    print("MODEL ARCHITECTURE")
    print("="*80)
    print(f"Input features: {config['node_features']}")
    print(f"Hidden layers: {config['nb_layer']} x {config['emb_hidden']} neurons")
    print(f"Embedding dimension: {config['emb_dim']}D latent space")
    print(f"Activation: {config['activation']}")
    print("="*80)
    print()
    
    print("="*80)
    print("TRAINING STRATEGY")
    print("="*80)
    print(f"Training approach: Contrastive learning with margin={config['margin']}")
    print(f"Points per batch: {config['points_per_batch']}")
    print(f"Training radius (normalized): {config['r_train']}")
    print(f"Max neighbors (train/val): {config['knn']}/{config['knn_val']}")
    print(f"Randomisation: {config['randomisation']}x random negative pairs")
    print("="*80)
    print()
    
    # Setup trainer
    trainer = Trainer(
        max_epochs=config['max_epochs'],
        accelerator=accelerator,
        devices=config.get('devices', 1),
        callbacks=[checkpoint_callback, early_stop_callback, loss_printer],
        logger=loggers,
        log_every_n_steps=config.get('wandb_log_every_n_batches', config.get('log_every_n_steps', 50)),
        check_val_every_n_epoch=config.get('check_val_every_n_epoch', 1),
        val_check_interval=config.get('val_check_interval'),
        enable_progress_bar=True,  # Show progress bar
        enable_model_summary=True,
    )
    
    # Train
    print("="*80)
    print("STARTING TRAINING")
    print("="*80)
    print()
    trainer.fit(stage_module)
    
    print("="*80)
    print("✓ TRAINING COMPLETE!")
    print("="*80)
    print(f"Best model: {checkpoint_callback.best_model_path}")
    print(f"Best {config['metric_to_monitor']}: {checkpoint_callback.best_model_score:.4f}")
    print()
    print("Next steps:")
    print(f"  python build_latent_graphs.py {Path(checkpoint_callback.best_model_path).stem}")
    print()


if __name__ == "__main__":
    main()
