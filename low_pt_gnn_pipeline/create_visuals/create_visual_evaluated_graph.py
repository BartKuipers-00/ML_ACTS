#!/usr/bin/env python3
"""
create_evaluated_visual.py - Visualize GNN edge classification results

Creates an interactive 3D HTML visualization showing how the GNN model classified
edges in a test graph. Edges are colored by classification correctness:
  - Dark green: Correctly classified true edges (True Positives)
  - Light gray: Correctly classified false edges (True Negatives)
  - Red: Incorrectly classified true edges (False Negatives - missed tracks)
  - Orange: Incorrectly classified false edges (False Positives - false connections)

Usage:
  python create_evaluated_visual.py <model_name> <dataset> <index> [--edge-cut THRESHOLD] [--output FILE]

Examples:
  python create_evaluated_visual.py epoch9_900 testset 1
    # Evaluates 1st test graph with model epoch9_900
  
  python create_evaluated_visual.py epoch9_900 valset 3
    # Evaluates 3rd validation graph
  
  python create_evaluated_visual.py epoch23_900 trainset 5 --edge-cut 0.3
    # Evaluates 5th training graph with threshold 0.3
  
  python create_evaluated_visual.py epoch9_900 testset 1 --output results/graph1.html
    # Custom output location
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import torch
import yaml

try:
    import plotly.graph_objects as go
except ImportError:
    print("Error: plotly is required. Install with: pip install plotly")
    sys.exit(1)

from acorn.stages.edge_classifier.models.interaction_gnn import InteractionGNN, InteractionGNN2


def gnn_class_for_ckpt(checkpoint_path):
    """Pick IGNN1 vs IGNN2 from checkpoint hparams (in_out_diff_agg is IGNN2-only)."""
    hp = torch.load(checkpoint_path, map_location="cpu", weights_only=False)["hyper_parameters"]
    return InteractionGNN2 if "in_out_diff_agg" in hp else InteractionGNN


from visual_utils import (
    validate_dataset_name,
    cylindrical_to_cartesian,
    get_standard_scene_layout,
    build_edge_coordinates,
    create_node_hover_text
)


def find_checkpoint(model_name, base_dir=None):
    """Find checkpoint file in saved_models directory."""
    if base_dir is None:
        base_dir = Path(__file__).resolve().parent.parent / "saved_models"

    checkpoint_path = base_dir / f"{model_name}.ckpt"
    if checkpoint_path.exists():
        return checkpoint_path

    raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}")


def evaluate_and_visualize(model_name, dataset_name, index, edge_cut=0.5, config_file=None, output_path=None):
    """
    Load model, evaluate on a specific graph, and create visualization.
    
    Args:
        model_name: Name of model checkpoint
        dataset_name: Name of dataset ('trainset', 'valset', or 'testset')
        index: Index of graph in dataset (1-based)
        edge_cut: Classification threshold
        config_file: Path to config file
        output_path: Output HTML file path
    """
    print("=" * 70)
    print(f"EVALUATING MODEL: {model_name} on {dataset_name} graph #{index}")
    print("=" * 70)
    print()
    
    # Validate dataset name
    validate_dataset_name(dataset_name)
    
    # Set default config file if not provided
    if config_file is None:
        script_dir = Path(__file__).resolve().parent
        config_file = script_dir.parent / 'acorn_configs' / 'gnn_stage_(2)' / 'gnn_train.yaml'
    else:
        config_file = Path(config_file)
        if not config_file.is_absolute():
            script_dir = Path(__file__).resolve().parent
            config_file = script_dir.parent / config_file
    
    # Load config
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    # Find checkpoint
    checkpoint_path = find_checkpoint(model_name)
    print(f"Found checkpoint: {checkpoint_path}")
    
    # Load checkpoint to get its hyperparameters (preserve architecture params)
    print("Loading model from checkpoint...")
    checkpoint = torch.load(str(checkpoint_path), map_location='cpu', weights_only=False)
    checkpoint_hparams = checkpoint.get('hyper_parameters', {})
    
    # Start with checkpoint hyperparameters (preserves architecture like 'hidden')
    # Then fill in any missing keys from config
    merged_hparams = {**checkpoint_hparams}
    for key, value in config.items():
        if key not in merged_hparams:
            merged_hparams[key] = value
    
    # Override paths and data-related params from current config
    script_dir = Path(__file__).resolve().parent
    input_dir = script_dir.parent / config['input_dir']
    merged_hparams['input_dir'] = str(input_dir)
    merged_hparams['stage_dir'] = str(input_dir)
    merged_hparams['data_split'] = config['data_split']
    merged_hparams['reprocess_classifier'] = True
    
    # Load model with merged hyperparameters
    model = gnn_class_for_ckpt(checkpoint_path).load_from_checkpoint(
        str(checkpoint_path),
        map_location='cuda' if torch.cuda.is_available() else 'cpu',
        **merged_hparams
    ) 
    
    # Load graph file directly from disk (sorted by filename) to match create_graph_visual.py
    dataset_dir = input_dir / dataset_name
    if not dataset_dir.exists():
        raise ValueError(f"Dataset directory not found: {dataset_dir}")
    
    # Get all .pyg files sorted (matching create_graph_visual.py behavior)
    graph_files = sorted([f for f in dataset_dir.glob('*.pyg')])
    if len(graph_files) == 0:
        raise ValueError(f"No graph files found in {dataset_dir}")
    
    if index < 1 or index > len(graph_files):
        raise ValueError(f"Index {index} out of range. Available: 1-{len(graph_files)}")
    
    # Select the graph file at the specified index (convert to 0-indexed)
    graph_path = graph_files[index - 1]
    print(f"Loading {dataset_name} graph #{index} from: {graph_path.name}")
    
    # Load the graph file
    graph = torch.load(graph_path, map_location='cpu', weights_only=False)
    
    # Setup model to get dataset for preprocessing
    if dataset_name == 'testset':
        model.setup(stage='test')
        loader = model.test_dataloader()
    else:
        model.setup(stage='fit')
        if dataset_name == 'trainset':
            loader = model.train_dataloader()
        else:  # valset
            loader = model.val_dataloader()
    
    dataset = loader.dataset
    
    # Apply preprocessing that the dataset would normally do
    # (matching GraphDataset.get() and preprocess_event() behavior)
    from acorn.utils.loading_utils import add_variable_name_prefix_in_pyg, infer_num_nodes
    
    if (not model.hparams.get("variable_with_prefix")) or model.hparams.get("add_variable_name_prefix_in_pyg"):
        graph = add_variable_name_prefix_in_pyg(graph)
    
    if dataset.preprocess:
        graph = dataset.preprocess_event(graph)
    
    print(f"Loaded {dataset_name} graph #{index}")
    print(f"  Nodes: {graph.num_nodes}")
    print(f"  Edges: {graph.edge_index.shape[1]}")
    print()
    
    # Run inference
    model.eval()
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model = model.to(device)
    
    print(f"Running inference on {device}...")
    with torch.no_grad():
        graph = graph.to(device)
        output = model(graph)
        scores = torch.sigmoid(output)
        predictions = (scores > edge_cut).cpu().numpy()
    
    # Get ground truth
    truth = graph.edge_y.cpu().numpy().astype(bool)
    predictions = predictions.astype(bool)
    
    # Classify edges into 4 categories
    true_positives = predictions & truth  # Correctly classified true edges
    true_negatives = ~predictions & ~truth  # Correctly classified false edges
    false_negatives = ~predictions & truth  # Missed true edges
    false_positives = predictions & ~truth  # False connections
    
    num_tp = true_positives.sum()
    num_tn = true_negatives.sum()
    num_fn = false_negatives.sum()
    num_fp = false_positives.sum()
    
    print(f"Classification results (threshold={edge_cut}):")
    print(f"  True Positives (TP):  {num_tp:>6} - Correctly found tracks (dark green)")
    print(f"  True Negatives (TN):  {num_tn:>6} - Correctly rejected false edges (light gray)")
    print(f"  False Negatives (FN): {num_fn:>6} - Missed tracks (red)")
    print(f"  False Positives (FP): {num_fp:>6} - False connections (orange)")
    print()
    
    accuracy = (num_tp + num_tn) / len(truth)
    precision = num_tp / (num_tp + num_fp) if (num_tp + num_fp) > 0 else 0
    recall = num_tp / (num_tp + num_fn) if (num_tp + num_fn) > 0 else 0
    
    print(f"Metrics:")
    print(f"  Accuracy:  {accuracy:.4f} ({accuracy*100:.2f}%)")
    print(f"  Precision: {precision:.4f} ({precision*100:.2f}%)")
    print(f"  Recall:    {recall:.4f} ({recall*100:.2f}%)")
    print()
    
    # Create visualization
    print("Creating visualization...")
    fig = create_classification_visualization(
        graph, 
        true_positives, 
        true_negatives, 
        false_negatives, 
        false_positives,
        model_name,
        dataset_name,
        index,
        edge_cut,
        hparams=model.hparams
    )
    
    # Determine output path
    if output_path is None:
        script_dir = Path(__file__).resolve().parent
        visuals_dir = script_dir.parent / 'data' / 'visuals' / 'evaluated' / dataset_name
        visuals_dir.mkdir(parents=True, exist_ok=True)
        output_path = visuals_dir / f"{model_name}_{dataset_name}{index:03d}_threshold{edge_cut}.html"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Writing interactive HTML to: {output_path}")
    fig.write_html(str(output_path), include_plotlyjs='cdn')
    
    print()
    print("=" * 70)
    print("✓ Visualization complete!")
    print("=" * 70)
    print(f"\nOpen in browser: {output_path.absolute()}")
    
    return output_path


def create_classification_visualization(graph, tp_mask, tn_mask, fn_mask, fp_mask, model_name, dataset_name, index, edge_cut, hparams=None):
    """
    Create interactive Plotly visualization with edges colored by classification.
    
    Args:
        graph: PyTorch Geometric Data object (may have scaled features)
        tp_mask: Boolean mask for true positives
        tn_mask: Boolean mask for true negatives
        fn_mask: Boolean mask for false negatives
        fp_mask: Boolean mask for false positives
        model_name: Name of the model
        dataset_name: Name of dataset ('trainset', 'valset', or 'testset')
        index: Index of graph in dataset
        edge_cut: Classification threshold
        hparams: Model hyperparameters (to check if features are scaled)
    """
    # Extract node positions - unscale if they were scaled during preprocessing
    # Check what feature names are actually in the graph (may be prefixed or not)
    if hasattr(graph, 'hit_r'):
        r = graph.hit_r.cpu().numpy()
        phi = graph.hit_phi.cpu().numpy()
        z = graph.hit_z.cpu().numpy()
    else:
        r = graph.r.cpu().numpy()
        phi = graph.phi.cpu().numpy()
        z = graph.z.cpu().numpy()
    
    # Unscale coordinates if they were scaled during preprocessing
    # The scale_features function scales based on node_features in config, which may or may not have prefix
    if hparams is not None and "node_scales" in hparams and "node_features" in hparams:
        node_features = hparams["node_features"]
        node_scales = hparams["node_scales"]
        
        # Find which features correspond to r, phi, z and unscale them
        for i, feature_name in enumerate(node_features):
            if i >= len(node_scales):
                continue
            
            scale = node_scales[i]
            
            # Check if this feature corresponds to r, phi, or z
            # Handle both prefixed (hit_r) and non-prefixed (r) names
            base_name = feature_name.replace("hit_", "").replace("_hit", "")
            
            if base_name == "r":
                r = r * scale
            elif base_name == "phi":
                phi = phi * scale
            elif base_name == "z":
                z = z * scale
    
    # Convert to Cartesian
    x, y, z_cart = cylindrical_to_cartesian(r, phi, z)
    num_nodes = len(x)
    
    # Build edges by category using visual_utils
    edge_index = graph.edge_index.cpu().numpy()
    edge_data = build_edge_coordinates(
        edge_index,
        x, y, z_cart,
        masks={
            'tn': tn_mask,
            'fp': fp_mask,
            'fn': fn_mask,
            'tp': tp_mask
        }
    )

    # Create figure
    fig = go.Figure()

    # Add edges for each category (order matters for layering)
    # 1. True Negatives (light gray, bottom layer, very transparent)
    if edge_data['tn'][3] > 0:
        edge_x, edge_y, edge_z, count = edge_data['tn']
        fig.add_trace(go.Scatter3d(
            x=edge_x, y=edge_y, z=edge_z,
            mode='lines',
            line=dict(color='#A9A9A9', width=0.5),
            opacity=0.05,
            name=f'True Negatives ({count})',
            hoverinfo='skip',
            showlegend=True
        ))

    # 2. False Positives (orange)
    if edge_data['fp'][3] > 0:
        edge_x, edge_y, edge_z, count = edge_data['fp']
        fig.add_trace(go.Scatter3d(
            x=edge_x, y=edge_y, z=edge_z,
            mode='lines',
            line=dict(color='orange', width=2),
            opacity=0.7,
            name=f'False Positives ({count})',
            hoverinfo='skip',
            showlegend=True
        ))

    # 3. False Negatives (red)
    if edge_data['fn'][3] > 0:
        edge_x, edge_y, edge_z, count = edge_data['fn']
        fig.add_trace(go.Scatter3d(
            x=edge_x, y=edge_y, z=edge_z,
            mode='lines',
            line=dict(color='red', width=2),
            opacity=0.7,
            name=f'False Negatives ({count})',
            hoverinfo='skip',
            showlegend=True
        ))

    # 4. True Positives (dark green, top layer)
    if edge_data['tp'][3] > 0:
        edge_x, edge_y, edge_z, count = edge_data['tp']
        fig.add_trace(go.Scatter3d(
            x=edge_x, y=edge_y, z=edge_z,
            mode='lines',
            line=dict(color='darkgreen', width=2.5),
            opacity=0.8,
            name=f'True Positives ({count})',
            hoverinfo='skip',
            showlegend=True
        ))
    
    # Add nodes with hover text
    hover_text = [
        create_node_hover_text(i, {
            'x': x[i], 'y': y[i], 'z_cart': z_cart[i],
            'r': r[i], 'phi': phi[i]
        })
        for i in range(num_nodes)
    ]
    
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z_cart,
        mode='markers',
        marker=dict(
            size=3,
            color='steelblue',
            opacity=0.6,
            line=dict(width=0.5, color='white')
        ),
        name=f'Nodes ({num_nodes})',
        text=hover_text,
        hoverinfo='text',
        showlegend=True
    ))
    
    # Calculate metrics for title
    total_edges = len(tp_mask)
    accuracy = (tp_mask.sum() + tn_mask.sum()) / total_edges
    precision = tp_mask.sum() / (tp_mask.sum() + fp_mask.sum()) if (tp_mask.sum() + fp_mask.sum()) > 0 else 0
    recall = tp_mask.sum() / (tp_mask.sum() + fn_mask.sum()) if (tp_mask.sum() + fn_mask.sum()) > 0 else 0
    
    # Update layout with standard scene configuration
    title_text = f"Model: {model_name} | {dataset_name.capitalize()} Graph #{index} | Threshold: {edge_cut}<br>"
    title_text += f"Accuracy: {accuracy:.3f} | Precision: {precision:.3f} | Recall: {recall:.3f}"

    fig.update_layout(**get_standard_scene_layout(title_text, margin_t=60))
    
    return fig


def main():
    parser = argparse.ArgumentParser(
        description='Visualize GNN edge classification results on a test graph',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python create_evaluated_visual.py epoch9_900 testset 1
    # Evaluate 1st test graph with model epoch9_900
  
  python create_evaluated_visual.py epoch9_900 valset 3
    # Evaluate 3rd validation graph
  
  python create_evaluated_visual.py epoch23_900 trainset 5 --edge-cut 0.3
    # Evaluate 5th training graph with threshold 0.3
  
  python create_evaluated_visual.py epoch9_900 testset 1 --output results/graph1.html
    # Save to custom location

Edge colors:
  Dark green  - True Positives (correctly found tracks)
  Light gray  - True Negatives (correctly rejected false edges)
  Red         - False Negatives (missed tracks)
  Orange      - False Positives (false connections)
        """
    )
    parser.add_argument(
        'model_name',
        type=str,
        help='Name of the model checkpoint (e.g., epoch9_900)'
    )
    parser.add_argument(
        'dataset',
        type=str,
        choices=['trainset', 'valset', 'testset'],
        help='Dataset to visualize: trainset, valset, or testset'
    )
    parser.add_argument(
        'index',
        type=int,
        help='Index of graph to visualize (1-based, e.g., 1 for first graph)'
    )
    parser.add_argument(
        '--edge-cut',
        type=float,
        default=0.5,
        help='Score threshold for edge classification (default: 0.5)'
    )
    parser.add_argument(
        '--config',
        type=str,
        default=None,
        help='Path to config file (default: ../acorn_configs/gnn_stage_(2)/gnn_train.yaml)'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        default=None,
        help='Output HTML file path (default: data/visuals/evaluated/testset/<model>_test<N>_threshold<T>.html)'
    )
    
    args = parser.parse_args()
    
    try:
        evaluate_and_visualize(
            args.model_name,
            args.dataset,
            args.index,
            args.edge_cut,
            args.config,
            args.output
        )
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
