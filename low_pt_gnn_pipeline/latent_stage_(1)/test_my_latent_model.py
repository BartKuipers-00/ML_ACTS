#!/usr/bin/env python3
"""
Evaluate metric learning (latent cluster learning) model on test set

Computes TP, FP, TN, FN and other metrics for graph construction performance.

Example:
    python test_my_latent_model.py saved_models/low_pt_latentmodel_100_mixed_f1=0.2752.ckpt --num-events 20 --dr-same-layer-cut 10 --knn 200 --r-max 0.4 --r-max-geometric 500 --primary-only
"""

#low_pt_latentmodel_val_loss_val_loss=0.0074.ckpt

import argparse
import sys
from pathlib import Path
import yaml
import torch
import numpy as np
from tqdm import tqdm
from torch_geometric.loader import DataLoader

# Add acorn to path
SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / 'acorn'))
sys.path.insert(0, str(PIPELINE_ROOT))

from acorn.stages.graph_construction.models.metric_learning import MetricLearning
from low_pt_custom_utils.graph_utils import build_edges, edge_truth_labels  # KNN + fast truth match
from acorn.stages.graph_construction.graph_construction_stage import EventDataset


def load_model(checkpoint_path):
    """Load trained metric learning model from checkpoint"""
    print(f"Loading model from: {checkpoint_path}")
    
    checkpoint = torch.load(checkpoint_path, map_location='cpu', weights_only=False)
    
    # Extract hyperparameters
    hparams = checkpoint['hyper_parameters']

    # node_scales must be saved in the checkpoint so test normalization matches
    # training exactly. No fallback: a missing value would silently mis-normalize.
    if "node_scales" not in hparams:
        raise KeyError(
            f"'node_scales' not found in checkpoint hyper_parameters of {checkpoint_path}. "
            "Cannot determine the normalization used at training time. "
            "Retrain with node_scales in the config, or test a checkpoint that stored it."
        )

    print(f"  Model architecture:")
    print(f"    Input features: {hparams['node_features']}")
    print(f"    Hidden layers: {hparams['nb_layer']} x {hparams['emb_hidden']}")
    print(f"    Embedding dim: {hparams['emb_dim']}D")
    print(f"    Activation: {hparams['activation']}")
    print()
    
    # Create model instance
    model = MetricLearning(hparams)
    model.load_state_dict(checkpoint['state_dict'])
    model.eval()
    
    return model, hparams


def build_particle_only_edges(batch):
    """Rebuild track_edges to span segment boundaries (all consecutive same-particle hits by time).
    Used when segmented=False so cross-segment edges count as true positives."""
    particle_ids = batch.hit_particle_id
    hit_t = batch.hit_t
    edges = []
    for pid in particle_ids.unique():
        if pid == 0:
            continue
        mask = particle_ids == pid
        idx = mask.nonzero(as_tuple=True)[0]
        times = hit_t[idx]
        order = torch.argsort(times)
        sorted_idx = idx[order]
        if len(sorted_idx) >= 2:
            src = sorted_idx[:-1]
            dst = sorted_idx[1:]
            edges.append(torch.stack([src, dst], dim=0))
    if not edges:
        return torch.zeros(2, 0, dtype=torch.long, device=particle_ids.device)
    return torch.cat(edges, dim=1)


def primary_hit_mask(batch):
    """Per-hit bool: True if the hit's particle is primary. Noise hits -> False."""
    pid = batch.hit_particle_id.long()
    tpid = batch.track_particle_id.long()
    tprim = batch.track_particle_primary.bool()
    order = torch.argsort(tpid)
    s_ids, s_prim = tpid[order], tprim[order]
    idx = torch.searchsorted(s_ids, pid).clamp(max=s_ids.numel() - 1)
    return s_prim[idx] & (s_ids[idx] == pid)  # exact-match guard


def evaluate_model(model, hparams, testset, knn_max=50, r_max=0.15, r_max_geometric=None, same_layer_cut="geometry", segmented=True, device='cpu', primary_only=False):
    """
    Evaluate metric learning model on test set
    
    Args:
        model: Trained MetricLearning model
        hparams: Model hyperparameters
        testset: Test dataset
        knn_max: Maximum number of neighbors for KNN graph construction
        r_max: Maximum radius for graph construction
        device: Device to run evaluation on
    
    Returns:
        Dictionary with evaluation metrics
    """
    model = model.to(device)
    model.eval()
    
    # Create data loader
    test_loader = DataLoader(testset, batch_size=1, num_workers=0, shuffle=False)
    
    # Accumulate metrics across all events
    total_tp = 0
    total_fp = 0
    total_tn = 0
    total_fn = 0
    total_pred_edges = 0
    total_true_edges = 0
    total_nodes = 0
    # Primary-only confusion (secondary true edges treated as don't-care)
    total_tp_p = 0
    total_fp_p = 0
    total_fn_p = 0
    total_true_edges_p = 0

    # Per-event metrics for analysis
    event_metrics = []
    
    print(f"Evaluating on {len(testset)} test events...")
    geo_str = f"{r_max_geometric} mm" if r_max_geometric is not None else "disabled"
    print(f"Graph construction parameters: r_max={r_max}, k_max={knn_max}, r_max_geometric={geo_str}, same_layer_cut={same_layer_cut}")
    print(f"segmented={segmented}  ({'cross-segment edges count as false' if segmented else 'cross-segment same-particle edges count as true'})")
    print()
    
    with torch.no_grad():
        for batch_idx, batch in enumerate(tqdm(test_loader, desc="Processing events")):
            batch = batch.to(device)
            
            # Get embeddings
            embedding = model.apply_embedding(batch)
            
            # Build graph using KNN in latent space
            pred_edges = build_edges(
                query=embedding,
                database=embedding,
                indices=None,
                r_max=r_max,
                k_max=knn_max,
                backend="FRNN",
            )

            # Apply physical geometric distance cut in 3D detector space
            # Note: EventDataset scales node_features by node_scales, so hit_r and hit_z
            # are NOT in mm. Use hit_x/hit_y (raw mm, not in node_features) and
            # unscale hit_z by node_scales[2] (typically 500).
            if r_max_geometric is not None:
                src, dst = pred_edges
                node_scales = hparams["node_scales"]
                hit_x = batch.hit_x.float()
                hit_y = batch.hit_y.float()
                hit_z_mm = batch.hit_z.float() * node_scales[2]
                dx = hit_x[src] - hit_x[dst]
                dy = hit_y[src] - hit_y[dst]
                dz = hit_z_mm[src] - hit_z_mm[dst]
                dist3d = torch.sqrt(dx**2 + dy**2 + dz**2)
                pred_edges = pred_edges[:, dist3d <= r_max_geometric]

            # Remove same-layer edges by detector identity: same layer iff same (volume_id,
            # layer_id). Geometry-exact (works for endcap disks); removes 0 true edges.
            if same_layer_cut == "geometry":
                if not hasattr(batch, "hit_volume_id"):
                    raise AttributeError(
                        "same_layer_cut='geometry' needs hit_volume_id/hit_layer_id on the graph. "
                        "Rebuild the feature store with these in hit_features (convert_csv_to_pyg_sets.yaml)."
                    )
                src, dst = pred_edges
                same_layer = (batch.hit_volume_id[src] == batch.hit_volume_id[dst]) & \
                             (batch.hit_layer_id[src] == batch.hit_layer_id[dst])
                pred_edges = pred_edges[:, ~same_layer]

            # Get truth edges from batch
            # If segmented=False, rebuild edges to span segment boundaries
            if not segmented:
                true_edges = build_particle_only_edges(batch)
            else:
                true_edges = batch.track_edges
            
            # Count unique true edges first (handle undirected case where edges might be duplicated)
            # True edges might be stored as undirected (each edge appears twice: (i,j) and (j,i))
            true_edges_unique = torch.unique(true_edges, dim=1)
            num_unique_true_edges = len(true_edges_unique[0])
            
            # Truth-label predicted edges. Fast int-key dedup+match (identical result to acorn's
            # graph_intersection, but ~20x faster on CPU — avoids its torch.unique(dim=1)).
            pred_edges, edge_y = edge_truth_labels(
                pred_edges,
                true_edges,
                num_nodes=embedding.shape[0],
                undirected=hparams.get("undirected", False),
            )
            
            # Convert to numpy for easier computation
            edge_y_np = edge_y.cpu().numpy().astype(bool)
            
            # For graph construction, we predict edges exist if they're in the graph
            # So: predicted = True for all edges in pred_edges
            predictions = np.ones(len(edge_y_np), dtype=bool)
            
            # Compute confusion matrix
            # TP: predicted edges that are true
            tp_raw = np.sum(predictions & edge_y_np)
            
            # Count how many unique true edges we found
            # If graph is undirected, pred_edges might contain duplicates of the same true edge
            # So we need to count unique true edges in our predictions
            true_pred_edges = pred_edges[:, edge_y_np]
            if true_pred_edges.shape[1] > 0:
                # Make edges canonical (smaller node first) to handle undirected duplicates
                true_pred_edges_canonical = torch.stack([
                    torch.min(true_pred_edges, dim=0)[0],
                    torch.max(true_pred_edges, dim=0)[0]
                ])
                true_pred_edges_unique = torch.unique(true_pred_edges_canonical, dim=1)
                tp = len(true_pred_edges_unique[0])
            else:
                tp = 0
            
            fp = np.sum(predictions & ~edge_y_np)
            
            # Calculate True Negatives: pairs of hits we correctly did NOT connect
            # Total possible edges in undirected graph: N*(N-1)/2
            num_nodes_val = batch.num_nodes
            if hasattr(num_nodes_val, 'item'):
                num_nodes_val = num_nodes_val.item()
            elif isinstance(num_nodes_val, (list, tuple)):
                num_nodes_val = num_nodes_val[0] if len(num_nodes_val) > 0 else 0
            
            total_possible_edges = num_nodes_val * (num_nodes_val - 1) // 2
            
            # Edges we created (pred_edges already deduped by edge_truth_labels)
            num_pred_edges_unique = pred_edges.shape[1]
            
            # False edges that exist (should not be connected)
            # Total false edges = total_possible - true_edges
            total_false_edges = total_possible_edges - num_unique_true_edges
            
            # TN = False edges we correctly did NOT create
            # = False edges that exist AND we didn't create
            # = total_false_edges - fp
            tn = max(0, total_false_edges - fp)  # Ensure non-negative
            
            fn = num_unique_true_edges - tp  # True edges not in our prediction

            # Primary-only: drop secondary true edges (don't-care, trained at weight 0)
            if primary_only:
                prim = primary_hit_mask(batch)
                true_edges_p = true_edges[:, prim[true_edges[0]]]
                num_true_p = torch.unique(true_edges_p, dim=1).shape[1]
                _, edge_y_p = edge_truth_labels(
                    pred_edges, true_edges_p, num_nodes=embedding.shape[0],
                    undirected=hparams.get("undirected", False),
                )
                tp_pred = pred_edges[:, edge_y_p.cpu().numpy().astype(bool)]
                if tp_pred.shape[1] > 0:
                    canon = torch.stack([tp_pred.min(0)[0], tp_pred.max(0)[0]])
                    tp_p = torch.unique(canon, dim=1).shape[1]
                else:
                    tp_p = 0
                fp_p = int(np.sum(~edge_y_np))  # secondary-true preds excluded (edge_y True)
                total_tp_p += tp_p
                total_fp_p += fp_p
                total_fn_p += num_true_p - tp_p
                total_true_edges_p += num_true_p

            # Accumulate
            total_tp += tp
            total_fp += fp
            total_tn += tn
            total_fn += fn
            total_pred_edges += len(pred_edges[0])
            total_true_edges += num_unique_true_edges  # Count unique true edges, not duplicates
            
            # Handle num_nodes - could be int or tensor
            num_nodes = batch.num_nodes
            if hasattr(num_nodes, 'item'):
                num_nodes = num_nodes.item()
            elif isinstance(num_nodes, (list, tuple)):
                num_nodes = num_nodes[0] if len(num_nodes) > 0 else 0
            total_nodes += num_nodes
            
            # Store per-event metrics
            # Handle num_nodes - could be int or tensor
            num_nodes = batch.num_nodes
            if hasattr(num_nodes, 'item'):
                num_nodes = num_nodes.item()
            elif isinstance(num_nodes, (list, tuple)):
                num_nodes = num_nodes[0] if len(num_nodes) > 0 else 0
            
            event_metrics.append({
                'event_id': batch.event_id[0] if hasattr(batch, 'event_id') and hasattr(batch.event_id, '__getitem__') else batch_idx,
                'num_nodes': num_nodes,
                'num_pred_edges': len(pred_edges[0]),
                'num_true_edges': len(true_edges[0]),
                'tp': int(tp),
                'fp': int(fp),
                'fn': int(fn),
                'precision': float(tp / (tp + fp)) if (tp + fp) > 0 else 0.0,
                'recall': float(tp / (tp + fn)) if (tp + fn) > 0 else 0.0,
            })
    
    # Compute aggregate metrics
    total_edges_evaluated = total_tp + total_fp + total_fn  # Note: TN=0 for graph construction
    
    precision = total_tp / (total_tp + total_fp) if (total_tp + total_fp) > 0 else 0.0
    recall = total_tp / (total_tp + total_fn) if (total_tp + total_fn) > 0 else 0.0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

    primary = None
    if primary_only:
        prec_p = total_tp_p / (total_tp_p + total_fp_p) if (total_tp_p + total_fp_p) > 0 else 0.0
        rec_p = total_tp_p / (total_tp_p + total_fn_p) if (total_tp_p + total_fn_p) > 0 else 0.0
        primary = {
            'tp': total_tp_p, 'fp': total_fp_p, 'fn': total_fn_p,
            'true_edges': total_true_edges_p, 'precision': prec_p, 'recall': rec_p,
        }

    return {
        'primary': primary,
        'total_tp': total_tp,
        'total_fp': total_fp,
        'total_tn': total_tn,
        'total_fn': total_fn,
        'total_pred_edges': total_pred_edges,
        'total_true_edges': total_true_edges,
        'total_nodes': total_nodes,
        'num_events': len(testset),
        'precision': precision,
        'recall': recall,
        'f1_score': f1_score,
        'event_metrics': event_metrics,
    }


def print_results(metrics, knn_max, r_max, r_max_geometric=None, same_layer_cut="geometry", segmented=True):
    """Print evaluation results in a formatted way"""
    print("\n" + "="*80)
    print("EVALUATION RESULTS")
    print("="*80)
    print(f"Graph construction parameters:")
    print(f"  r_max: {r_max}")
    print(f"  k_max: {knn_max}")
    geo_str = f"{r_max_geometric} mm" if r_max_geometric is not None else "disabled"
    print(f"  r_max_geometric: {geo_str}")
    print(f"  same_layer_cut: {same_layer_cut}")
    print(f"  segmented: {segmented}  ({'within-segment edges only' if segmented else 'cross-segment same-particle edges are true'})")
    print()
    # Calculate average connections per node
    avg_connections = metrics['total_pred_edges'] / metrics['total_nodes'] if metrics['total_nodes'] > 0 else 0

    print(f"Dataset:")
    print(f"  Number of events: {metrics['num_events']}")
    print(f"  Total nodes: {metrics['total_nodes']:,}")
    print(f"  Total true edges: {metrics['total_true_edges']:,}")
    print(f"  Total predicted edges: {metrics['total_pred_edges']:,}")
    print(f"  Avg connections/node: {avg_connections:.1f}")
    print()
    print("="*80)
    print("CONFUSION MATRIX")
    print("="*80)
    print(f"  True Positives (TP):  {metrics['total_tp']:>12,} ")
    print(f"  False Positives (FP): {metrics['total_fp']:>12,} ")
    print(f"  True Negatives (TN):  {metrics['total_tn']:>12,}  ")
    print(f"  False Negatives (FN): {metrics['total_fn']:>12,} ")
    print()
    print("="*80)
    print("PERFORMANCE METRICS")
    print("="*80)
    print(f"  Precision:  {metrics['precision']:.6f} ({metrics['precision']*100:.2f}%)")
    print(f"              TP / (TP + FP) = {metrics['total_tp']:,} / {metrics['total_tp'] + metrics['total_fp']:,}")
    print()
    print(f"  Recall:     {metrics['recall']:.6f} ({metrics['recall']*100:.2f}%)")
    print(f"              TP / (TP + FN) = {metrics['total_tp']:,} / {metrics['total_tp'] + metrics['total_fn']:,}")
    print()

    # Primary-only vs inclusive comparison (secondary true edges excluded)
    p = metrics.get('primary')
    if p is not None:
        n_sec = metrics['total_true_edges'] - p['true_edges']
        print("="*80)
        print("PRIMARY-ONLY ")
        print("="*80)
        print(f"  Secondary true edges excluded: {n_sec:,} of {metrics['total_true_edges']:,}")
        print(f"  Recall:    {p['recall']*100:.2f}%  (inclusive: {metrics['recall']*100:.2f}%, "
              f"Δ {(p['recall']-metrics['recall'])*100:+.2f})")
        print(f"  Precision: {p['precision']*100:.4f}%  (inclusive: {metrics['precision']*100:.4f}%)")
        print()
    print("="*80)
    
    # Per-event statistics
    if metrics['event_metrics']:
        precisions = [m['precision'] for m in metrics['event_metrics']]
        recalls = [m['recall'] for m in metrics['event_metrics']]
        print("PER-EVENT STATISTICS")
        print("="*80)
        print(f"  Precision: mean={np.mean(precisions):.4f}, std={np.std(precisions):.4f}")
        print(f"             min={np.min(precisions):.4f}, max={np.max(precisions):.4f}")
        print(f"  Recall:    mean={np.mean(recalls):.4f}, std={np.std(recalls):.4f}")
        print(f"             min={np.min(recalls):.4f}, max={np.max(recalls):.4f}")
        print("="*80)
    print()


def main():
    parser = argparse.ArgumentParser(
        description='Evaluate metric learning model on test set',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python test_my_latent_model.py saved_models/low_pt_latent_f1=0.0149.ckpt
  python test_my_latent_model.py saved_models/low_pt_latent_f1=0.0149.ckpt --knn 1000 --r-max 0.15
  python test_my_latent_model.py saved_models/low_pt_latent_f1=0.0149.ckpt --knn 2000 --r-max 0.15
        """
    )
    parser.add_argument(
        'checkpoint_path',
        type=str,
        help='Path to model checkpoint (.ckpt file)'
    )
    parser.add_argument(
        '--knn', '--k-max',
        type=int,
        default=None,
        help='Maximum number of neighbors for KNN graph construction (default: from graph_construction_latent.yaml)'
    )
    parser.add_argument(
        '--r-max',
        type=float,
        default=None,
        help='Maximum radius for graph construction (default: from graph_construction_latent.yaml)'
    )
    parser.add_argument(
        '--r-max-geometric',
        type=float,
        default=None,
        help='Maximum 3D Euclidean distance (mm) between connected hits; applied after KNN (default: from graph_construction_latent.yaml)'
    )
    parser.add_argument(
        '--same-layer-cut',
        type=str,
        default=None,
        choices=['geometry', 'none'],
        help="Same-layer edge removal: 'geometry' = drop edges sharing (volume_id, layer_id); 'none' = disabled (default: from graph_construction_latent.yaml)"
    )
    parser.add_argument(
        '--segmented',
        type=lambda x: x.lower() not in ('false', '0', 'no'),
        default=None,
        metavar='BOOL',
        help='Whether cross-segment same-particle edges are false (default: from graph_construction_latent.yaml)'
    )
    parser.add_argument(
        '--num-events',
        type=int,
        default=None,
        help='Number of test events to evaluate (default: all events in testset as defined by data_split in config)'
    )
    parser.add_argument(
        '--primary-only',
        action='store_true',
        help='Also report recall/precision with secondary true edges excluded (weight-0 in training)'
    )
    parser.add_argument(
        '--config',
        type=str,
        default=None,
        help='Path to training config file (default: acorn_configs/latent_stage_(1)/latent_cluster_learning_train.yaml)'
    )
    parser.add_argument(
        '--graph-config',
        type=str,
        default=None,
        help='Path to graph construction config file (default: acorn_configs/latent_stage_(1)/graph_construction_latent.yaml)'
    )
    parser.add_argument(
        '--device',
        type=str,
        default='auto',
        help='Device to use: auto (try GPU, fallback to CPU), cpu, or cuda (default: auto)'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("METRIC LEARNING MODEL EVALUATION")
    print("="*80)
    print()
    
    # Determine device - try GPU, fallback to CPU
    if args.device == 'auto':
        if torch.cuda.is_available():
            device = 'cuda'
            print(f"✓ GPU detected: {torch.cuda.get_device_name(0)}")
            print(f"  Using CUDA for evaluation")
        else:
            device = 'cpu'
            print("⚠ GPU not available, falling back to CPU")
    else:
        device = args.device
        if device == 'cuda' and not torch.cuda.is_available():
            print("⚠ WARNING: CUDA requested but not available, falling back to CPU")
            device = 'cpu'
        print(f"Using device: {device}")
    print()
    
    # Load checkpoint
    checkpoint_path = Path(args.checkpoint_path)
    if not checkpoint_path.is_absolute():
        checkpoint_path = PIPELINE_ROOT / checkpoint_path
    
    if not checkpoint_path.exists():
        print(f"ERROR: Checkpoint not found: {checkpoint_path}")
        sys.exit(1)
    
    model, hparams = load_model(checkpoint_path)
    
    # Load training config
    if args.config is None:
        config_path = PIPELINE_ROOT / 'acorn_configs' / 'latent_stage_(1)' / 'latent_cluster_learning_train.yaml'
    else:
        config_path = Path(args.config)

    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}")
        sys.exit(1)

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Load graph construction config
    if args.graph_config is None:
        graph_config_path = PIPELINE_ROOT / 'acorn_configs' / 'latent_stage_(1)' / 'graph_construction_latent.yaml'
    else:
        graph_config_path = Path(args.graph_config)

    if not graph_config_path.exists():
        print(f"ERROR: Graph construction config not found: {graph_config_path}")
        sys.exit(1)

    with open(graph_config_path, 'r') as f:
        graph_config = yaml.safe_load(f)

    # Resolve evaluation parameters: CLI args override graph construction config
    knn_max           = args.knn                if args.knn                is not None else graph_config.get('k_max', 1000)
    r_max             = args.r_max              if args.r_max              is not None else graph_config.get('r_max', 0.15)
    r_max_geo         = args.r_max_geometric    if args.r_max_geometric    is not None else graph_config.get('r_max_geometric', None)
    same_layer_cut    = args.same_layer_cut     if args.same_layer_cut     is not None else graph_config.get('same_layer_cut', 'geometry')
    segmented         = args.segmented          if args.segmented          is not None else graph_config.get('segmented', True)

    # Setup dataset - path is relative to pipeline root
    input_dir_config = config['input_dir']

    # Construct path relative to pipeline root (not script dir)
    input_dir = PIPELINE_ROOT / input_dir_config
    testset_dir = input_dir / 'testset'

    if not testset_dir.exists():
        print(f"ERROR: Test set directory not found:")
        print(f"  {testset_dir}")
        print(f"\nPlease check your config file 'input_dir' setting or directory structure.")
        sys.exit(1)
    
    print(f"Loading test set from: {testset_dir}")
    print(f"Using input_dir: {input_dir}")
    print()
    
    # Create dataset
    # Get number of test events from config
    data_split = config.get('data_split', [9000, 500, 500])
    num_test_events = args.num_events if args.num_events is not None else (data_split[2] if len(data_split) > 2 else None)
    
    testset = EventDataset(
        input_dir=str(input_dir),
        data_name='testset',
        num_events=num_test_events,
        hparams=hparams,
    )
    
    print(f"Loaded {len(testset)} test events")
    print()
    
    # Evaluate
    metrics = evaluate_model(
        model,
        hparams,
        testset,
        knn_max=knn_max,
        r_max=r_max,
        r_max_geometric=r_max_geo,
        same_layer_cut=same_layer_cut,
        segmented=segmented,
        device=device,
        primary_only=args.primary_only,
    )

    # Print results
    print_results(metrics, knn_max, r_max, r_max_geo, same_layer_cut=same_layer_cut, segmented=segmented)
    
    print("✓ Evaluation complete!")


if __name__ == '__main__':
    main()
