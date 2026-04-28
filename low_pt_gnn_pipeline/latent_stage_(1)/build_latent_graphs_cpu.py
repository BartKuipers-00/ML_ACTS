#!/usr/bin/env python3
"""
Graph construction using learned latent space embeddings

This script uses a trained metric learning model to map hits into a learned
latent space, then builds edges using KNN in that space. This is more intelligent
than simple radius + KNN because the model learns which hits are likely to be
connected.

Usage:
    python build_latent_graphs.py <model_name>
    
Example:
    python build_latent_graphs.py latent_builder_10epochs
    python build_latent_graphs.py last

The script looks for <model_name>.ckpt in saved_models/ directory.
Configuration is loaded from acorn_configs/graph_construction_latent.yaml
"""

import sys
import os
from pathlib import Path
import yaml
import torch
import pandas as pd
import numpy as np
from tqdm import tqdm

# Add acorn to path
SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / 'acorn'))
sys.path.insert(0, str(PIPELINE_ROOT))

from acorn.stages.graph_construction.models.metric_learning import MetricLearning
from low_pt_custom_utils.graph_utils import compute_edge_y


def load_config(config_path=None):
    """Load configuration from YAML file"""
    config_path = PIPELINE_ROOT / 'acorn_configs' / 'latent_stage_(1)' / 'graph_construction_latent.yaml'
    
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Make paths absolute relative to pipeline root directory
    for path_key in ['input_dir', 'output_dir']:
        if path_key in config:
            config[path_key] = str(PIPELINE_ROOT / config[path_key])
    
    return config


def load_model(checkpoint_path):
    """
    Load trained metric learning model from checkpoint
    
    Args:
        checkpoint_path: Path to the .ckpt file
    
    Returns:
        model: Loaded MetricLearning model in eval mode
        config: Model configuration from checkpoint
    """
    print(f"Loading model from: {checkpoint_path}")
    
    # Load checkpoint
    checkpoint = torch.load(checkpoint_path, map_location='cpu')
    
    # Extract hyperparameters
    hparams = checkpoint['hyper_parameters']
    
    print(f"  Model architecture:")
    print(f"    Input features: {hparams['node_features']}")
    print(f"    Hidden layers: {hparams['nb_layer']} x {hparams['emb_hidden']}")
    print(f"    Embedding dim: {hparams['emb_dim']}D")
    print(f"    Activation: {hparams['activation']}")
    print()
    
    # Create model instance
    model = MetricLearning(hparams)
    # Load weights
    model.load_state_dict(checkpoint['state_dict'])
    model.eval()
    
    return model, hparams


def build_edges_latent(model, graph, node_features, k_max=500, r_max=0.4, device='cpu'):
    """
    Build edges using learned embeddings + KNN

    Args:
        model: Trained MetricLearning model
        graph: PyG graph with node features
        node_features: List of feature names to use (e.g., ['hit_r', 'hit_phi', 'hit_z'])
        k_max: Maximum number of neighbors per hit
        r_max: Maximum radius for edge construction in latent space
        device: Device to run model on

    Returns:
        edge_index: Tensor of shape [2, num_edges] with candidate edges
        embeddings: Tensor of shape [num_hits, emb_dim] with learned embeddings
    """
    from acorn.stages.graph_construction.models.utils import build_edges

    model = model.to(device)

    # Extract features from graph
    feature_list = []
    for feat_name in node_features:
        if hasattr(graph, feat_name):
            feature_list.append(getattr(graph, feat_name))
        else:
            raise ValueError(f"Feature {feat_name} not found in graph. Available: {graph.keys}")

    # Stack features into tensor [num_hits, num_features]
    x = torch.stack(feature_list, dim=1).float().to(device)

    with torch.no_grad():
        # Get embeddings from model - pass feature tensor, not graph
        embeddings = model(x)  # Shape: [num_hits, emb_dim]

    # Build edges using KNN in latent space (using ACORN's build_edges with FRNN fallback)
    edge_index = build_edges(
        query=embeddings,
        database=embeddings,
        indices=None,
        r_max=r_max,
        k_max=k_max,
        backend="FRNN",
    )

    return edge_index.cpu(), embeddings.cpu()



def process_event(model, node_features, input_path, output_path, truth_csv_path, k_max=500, r_max=0.4, device='cpu', config=None):
    """Process a single event file using learned embeddings"""
    graph = torch.load(input_path)

    # Load hit coordinates and particle IDs from truth CSV (ground truth)
    truth_df = pd.read_csv(truth_csv_path)
    hit_particle_ids = torch.tensor(truth_df['particle_id'].values, dtype=torch.long)

    # Extract Cartesian coordinates from CSV
    x = truth_df['x'].values
    y = truth_df['y'].values
    z_cart = truth_df['z'].values

    # Compute cylindrical coordinates (r, phi, z)
    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)

    # Convert to tensors
    r_tensor = torch.tensor(r, dtype=torch.float32)
    phi_tensor = torch.tensor(phi, dtype=torch.float32)
    z_tensor = torch.tensor(z_cart, dtype=torch.float32)

    # Store coordinates on graph for later use (with 'hit_' prefix for model)
    graph.hit_r = r_tensor.double()
    graph.hit_phi = phi_tensor.double()
    graph.hit_z = z_tensor.double()
    graph.r = r_tensor.double()
    graph.phi = phi_tensor.double()
    graph.z = z_tensor.double()
    graph.x = torch.tensor(x, dtype=torch.float64)

    # Build candidate edges using learned embeddings
    edge_index, embeddings = build_edges_latent(model, graph, node_features, k_max, r_max, device)

    # Remove same-layer edges: hits with |Δr| < dr_same_layer_cut are on the same detector layer
    dr_same_layer_cut = config.get('dr_same_layer_cut', None) if config else None
    if dr_same_layer_cut is not None:
        hit_r = torch.tensor(np.sqrt(truth_df['x'].values**2 + truth_df['y'].values**2), dtype=torch.float32)
        dr = torch.abs(hit_r[edge_index[0]] - hit_r[edge_index[1]])
        edge_index = edge_index[:, dr >= dr_same_layer_cut]

    # Compute truth labels
    segmented = config.get('segmented', True) if config else True
    one_in_one_out = config.get('one_in_one_out', False) if config else False
    hit_segment_ids = graph.hit_segment_id.long() if (hasattr(graph, 'hit_segment_id') and segmented) else None
    edge_y = compute_edge_y(
        edge_index,
        hit_particle_ids,
        graph.track_edges,
        segment_id=hit_segment_ids,
        one_in_one_out=one_in_one_out,
    )

    # Add to graph
    graph.edge_index = edge_index
    graph.edge_y = edge_y
    graph.particle_id = hit_particle_ids
    graph.embeddings = embeddings  # Store embeddings for analysis
    
    # Save
    torch.save(graph, output_path)
    
    return graph.num_nodes, edge_index.shape[1], edge_y.sum().item()


def run_graph_construction(model, hparams, config):
    """
    Build edges using learned embeddings + KNN
    """
    print("="*80)
    print("GRAPH CONSTRUCTION - LEARNED LATENT SPACE")
    print("="*80)
    
    # Extract parameters from config
    input_dir = Path(config['input_dir'])
    output_dir = Path(config['output_dir'])
    k_max = config['k_max']
    r_max = config['r_max']
    device = config.get('device', 'cpu')
    datasets_config = config.get('datasets', ['trainset', 'valset', 'testset'])
    if isinstance(datasets_config, dict):
        datasets = list(datasets_config.keys())
        dataset_limits = datasets_config
    else:
        datasets = datasets_config
        dataset_limits = {}
    node_features = hparams['node_features']

    print(f"Input:  {input_dir}")
    print(f"Output: {output_dir}")
    print(f"Max neighbors (KNN): {k_max}")
    print(f"Max radius: {r_max}")
    print(f"Device: {device}")
    print(f"Node features: {node_features}")
    print()
    
    for dataset_name in datasets:
        input_dataset_dir = input_dir / dataset_name
        output_dataset_dir = output_dir / dataset_name
        
        if not input_dataset_dir.exists():
            print(f"Skipping {dataset_name} - directory not found")
            continue
            
        output_dataset_dir.mkdir(parents=True, exist_ok=True)
        
        # Get all event files
        event_files = sorted([f.name for f in input_dataset_dir.glob('*-graph.pyg')])

        if len(event_files) == 0:
            print(f"Skipping {dataset_name} - no graph files found")
            continue

        limit = dataset_limits.get(dataset_name)
        if limit is not None and len(event_files) > limit:
            event_files = event_files[:limit]
            print(f"\nProcessing {dataset_name}: {len(event_files)} events (limited to {limit})")
        else:
            print(f"\nProcessing {dataset_name}: {len(event_files)} events")
        
        total_nodes = 0
        total_edges = 0
        total_true_edges = 0
        
        for event_file in tqdm(event_files, desc=f"  {dataset_name}"):
            input_path = input_dataset_dir / event_file
            output_path = output_dataset_dir / event_file
            # Convert Path to string for replace operation
            truth_csv_path = str(input_path).replace('-graph.pyg', '-truth.csv')
            
            num_nodes, num_edges, num_true = process_event(
                model, node_features, input_path, output_path, truth_csv_path, k_max, r_max, device, config=config
            )
            total_nodes += num_nodes
            total_edges += num_edges
            total_true_edges += num_true
        
        # Print statistics
        print(f"  ✓ Total nodes: {total_nodes}")
        print(f"  ✓ Total edges: {total_edges}")
        print(f"  ✓ True edges: {total_true_edges} ({100*total_true_edges/total_edges:.1f}%)")
        print(f"  ✓ Avg edges/node: {total_edges/total_nodes:.1f}")
    
    print("\n✓ Graph construction complete!\n")


def main():
    """Run graph construction with learned embeddings"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Build graphs using learned latent space",
        epilog="Example: python build_latent_graphs.py latent_builder_10epochs"
    )
    parser.add_argument(
        "model_name",
        type=str,
        nargs='?',
        default=None,
        help="Model name (without .ckpt extension) from saved_models/ directory"
    )

    args = parser.parse_args()
    
    print("\n" + "="*80)
    print("GNN TRAINING PIPELINE - LEARNED GRAPH CONSTRUCTION")
    print("="*80 + "\n")
    
    # Load configuration
    config = load_config(SCRIPT_DIR / 'acorn_configs' / 'graph_construction_latent.yaml')
    print("="*80)
    print("CONFIGURATION")
    print("="*80)
    print(yaml.dump(config))
    print("="*80 + "\n")
    
    # Get checkpoint name (command line overrides config)
    checkpoint_name = args.model_name or config.get('checkpoint')
    if not checkpoint_name:
        print("ERROR: No model specified.")
        print("Usage: python build_latent_graphs.py <model_name>")
        print("\nExample: python build_latent_graphs.py latent_builder_10epochs")
        sys.exit(1)
    
    # Add .ckpt extension if not present
    if not checkpoint_name.endswith('.ckpt'):
        checkpoint_name = checkpoint_name + '.ckpt'
    
    checkpoint_path = PIPELINE_ROOT / 'saved_models' / checkpoint_name

    if not checkpoint_path.exists():
        print(f"ERROR: Checkpoint not found: {checkpoint_path}")
        print(f"\nAvailable checkpoints in saved_models/:")
        saved_models_dir = PIPELINE_ROOT / 'saved_models'
        if saved_models_dir.exists():
            for ckpt in saved_models_dir.glob('*.ckpt'):
                print(f"  - {ckpt.name}")
        sys.exit(1)
    
    model, hparams = load_model(checkpoint_path)
    
    # Run graph construction
    run_graph_construction(model, hparams, config)
    
    print("="*80)
    print("✓ PIPELINE COMPLETE!")
    print("="*80)
    print(f"\nProcessed graphs are in: {config['output_dir']}")
    print(f"Ready for training with: python train_with_loss_logging.py")
    print(f"  (Update config to use input_dir: {config['output_dir']})\n")


if __name__ == "__main__":
    main()
