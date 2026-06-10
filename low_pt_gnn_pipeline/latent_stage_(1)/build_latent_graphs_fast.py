#!/usr/bin/env python3
"""
Graph construction using learned latent space embeddings


Usage:
    python build_latent_graphs_ultra_fast.py <model_name>

Example:
    python build_latent_graphs_ultra_fast.py latent_builder_10epochs
    python build_latent_graphs_ultra_fast.py last
"""

import sys
import os
from pathlib import Path
import yaml
import torch
import pandas as pd
import numpy as np
from tqdm import tqdm
from torch.utils.data import Dataset, DataLoader
from torch_geometric.data import Batch
import threading
import queue

# Add acorn to path
SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / 'acorn'))
sys.path.insert(0, str(PIPELINE_ROOT))

from acorn.stages.graph_construction.models.metric_learning import MetricLearning
from acorn.stages.graph_construction.models.utils import graph_intersection
from low_pt_custom_utils.graph_utils import build_edges, compute_edge_y  # proper KNN+radius (no FRNN needed)


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
    """Load trained metric learning model from checkpoint"""
    print(f"Loading model from: {checkpoint_path}")

    checkpoint = torch.load(checkpoint_path, map_location='cpu')
    hparams = checkpoint['hyper_parameters']

    print(f"  Model architecture:")
    print(f"    Input features: {hparams['node_features']}")
    print(f"    Hidden layers: {hparams['nb_layer']} x {hparams['emb_hidden']}")
    print(f"    Embedding dim: {hparams['emb_dim']}D")
    print(f"    Activation: {hparams['activation']}")
    print()

    model = MetricLearning(hparams)
    model.load_state_dict(checkpoint['state_dict'])
    model.eval()

    return model, hparams


class GraphDataset(Dataset):
    """Dataset that loads graphs and truth CSVs on-the-fly"""

    def __init__(self, event_files, input_dir, node_features):
        self.event_files = event_files
        self.input_dir = Path(input_dir)
        self.node_features = node_features

    def __len__(self):
        return len(self.event_files)

    def __getitem__(self, idx):
        event_file = self.event_files[idx]
        input_path = self.input_dir / event_file
        truth_csv_path = str(input_path).replace('-graph.pyg', '-truth.csv')

        # Load graph
        graph = torch.load(input_path)

        # Load truth CSV
        truth_df = pd.read_csv(truth_csv_path)
        hit_particle_ids = torch.tensor(truth_df['particle_id'].values, dtype=torch.long)

        # Extract Cartesian coordinates
        x = truth_df['x'].values
        y = truth_df['y'].values
        z_cart = truth_df['z'].values

        # Compute cylindrical coordinates
        r = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)

        # Store coordinates on graph
        graph.hit_r = torch.tensor(r, dtype=torch.float32)
        graph.hit_phi = torch.tensor(phi, dtype=torch.float32)
        graph.hit_z = torch.tensor(z_cart, dtype=torch.float32)
        graph.r = graph.hit_r.double()
        graph.phi = graph.hit_phi.double()
        graph.z = graph.hit_z.double()
        graph.x = torch.tensor(x, dtype=torch.float64)
        graph.hit_x = torch.tensor(x, dtype=torch.float32)
        graph.hit_y = torch.tensor(y, dtype=torch.float32)
        graph.particle_id = hit_particle_ids

        # Detector layer identity for the geometry-exact same-layer cut. Prefer baked pyg
        # attrs (rebuilt graphs); fall back to the truth_df already loaded here.
        if not hasattr(graph, 'hit_volume_id'):
            graph.hit_volume_id = torch.tensor(truth_df['volume_id'].values, dtype=torch.long)
        if not hasattr(graph, 'hit_layer_id'):
            graph.hit_layer_id = torch.tensor(truth_df['layer_id'].values, dtype=torch.long)

        # Segment ID might not exist in older datasets - default to None
        if not hasattr(graph, 'hit_segment_id'):
            graph.segment_id = None
        else:
            graph.segment_id = graph.hit_segment_id.long()

        # Store filename for saving later
        graph.filename = event_file

        return graph


def build_edges_latent_true_batch(model, graphs, node_features, k_max=500, r_max=0.15, r_max_geometric=None, same_layer_cut="geometry", segmented=True, one_in_one_out=False, node_scales=None, compute_truth=True):
    """
    Build edges using TRUE batching - process all graphs simultaneously on GPU

    Args:
        model: Trained MetricLearning model (already on correct device)
        graphs: List of PyG graphs
        node_features: List of feature names
        k_max: Maximum neighbors per hit
        r_max: Maximum radius in latent space
        compute_truth: If True (default), also compute edge_y and track_to_edge_map
            (needed for training/eval). Pass False at inference to skip them — they
            are unused downstream and track_to_edge_map costs an O(K*M) loop.

    Returns:
        graphs: List of graphs with edge_index and embeddings added; also edge_y
            and track_to_edge_map when compute_truth=True.
    """
    from torch_geometric.data import Batch

    device = next(model.parameters()).device

    # Extract features from all graphs
    all_features = []
    all_particle_ids = []
    all_segment_ids = []

    for graph in graphs:
        # Extract features for this graph
        feature_list = []
        for i, feat_name in enumerate(node_features):
            if hasattr(graph, feat_name):
                feat = getattr(graph, feat_name).float()
                if node_scales is not None:
                    feat = feat / node_scales[i]
                feature_list.append(feat)
            else:
                raise ValueError(f"Feature {feat_name} not found in graph")

        # Stack features [num_hits, num_features]
        x = torch.stack(feature_list, dim=1).float()
        all_features.append(x)
        all_particle_ids.append(graph.particle_id)
        all_segment_ids.append(graph.segment_id)

    # Concatenate all features into single batch
    batch_features = torch.cat(all_features, dim=0).to(device)

    # Create batch assignment (which graph each node belongs to)
    batch_assignment = torch.cat([
        torch.full((x.shape[0],), i, dtype=torch.long)
        for i, x in enumerate(all_features)
    ]).to(device)

    # Get embeddings for entire batch at once
    with torch.no_grad():
        all_embeddings = model(batch_features)  # [total_hits_in_batch, emb_dim]

    # Split embeddings back to individual graphs (keep on GPU for KNN matmul)
    node_offsets = [0]
    for x in all_features:
        node_offsets.append(node_offsets[-1] + x.shape[0])

    processed_graphs = []
    for i, graph in enumerate(graphs):
        start_idx = node_offsets[i]
        end_idx = node_offsets[i + 1]
        num_nodes = end_idx - start_idx

        # Extract embeddings for this graph (still on GPU)
        embeddings = all_embeddings[start_idx:end_idx]

        # Build edges — matmul is device-aware (CPU or GPU depending on config)
        graph_edges = build_edges(
            query=embeddings,
            database=embeddings,
            indices=None,
            r_max=r_max,
            k_max=k_max,
            backend="FRNN",
        )

        # Ensure edge list is on CPU before geometric cuts (hit tensors are on CPU)
        graph_edges = graph_edges.cpu()

        # Apply physical geometric distance cut in 3D detector space
        # Use raw Cartesian coordinates (mm) — hit_x/hit_y are unscaled raw CSV values
        if r_max_geometric is not None:
            src, dst = graph_edges
            hit_x = graph.hit_x
            hit_y = graph.hit_y
            hit_z = graph.hit_z
            dx = hit_x[src] - hit_x[dst]
            dy = hit_y[src] - hit_y[dst]
            dz = hit_z[src] - hit_z[dst]
            dist3d = torch.sqrt(dx**2 + dy**2 + dz**2)
            geo_mask = dist3d <= r_max_geometric
            graph_edges = graph_edges[:, geo_mask]

        # Remove same-layer edges by detector identity: two hits are on the same layer iff
        # they share (volume_id, layer_id). Geometry-exact — unlike a |Δr| proxy it works for
        # endcap disks (constant-z, span r) and removes 0 true edges (consecutive hits never
        # share a layer). See compare_same_layer_cut.py.
        if same_layer_cut == "geometry":
            src, dst = graph_edges
            same_layer = (graph.hit_volume_id[src] == graph.hit_volume_id[dst]) & \
                         (graph.hit_layer_id[src] == graph.hit_layer_id[dst])
            graph_edges = graph_edges[:, ~same_layer]
        elif same_layer_cut not in (None, "none", 0, 0.0, False):
            # Backward-compatible |Δr| cut (mm), for graphs without baked (volume_id, layer_id)
            # — e.g. the in-memory benchmark and older feature_store datasets. hit_r is raw mm.
            thr = float(same_layer_cut)
            src, dst = graph_edges
            dr = (graph.hit_r[src] - graph.hit_r[dst]).abs()
            graph_edges = graph_edges[:, dr >= thr]

        graph.edge_index = graph_edges
        graph.embeddings = embeddings.cpu()

        if not compute_truth:        # inference: skip unused edge_y + track_to_edge_map
            processed_graphs.append(graph)
            continue

        # Compute edge truth labels
        particle_id = all_particle_ids[i]
        segment_id = all_segment_ids[i] if segmented else None

        edge_y = compute_edge_y(
            graph_edges,
            particle_id,
            graph.track_edges,
            segment_id=segment_id,
            one_in_one_out=one_in_one_out,
        )

        # track_to_edge_map (ACORN standard): 1-D, length = num track_edges. Entry i =
        # index of true edge i in edge_index, or -1 if that true edge wasn't constructed.
        # undirected=True matches the GNN's undirected config and compute_edge_y's
        # undirected edge labelling. Consumed by the GNN stage's weighting (maps
        # track-level attrs like track_particle_primary onto edges).
        track_to_edge_map = graph_intersection(
            graph_edges,
            graph.track_edges,
            return_y_pred=False,
            return_truth_to_pred=True,
            undirected=True,
        )

        # Add truth attributes (edge_index/embeddings already set above)
        graph.edge_y = edge_y
        graph.track_to_edge_map = track_to_edge_map

        processed_graphs.append(graph)

    return processed_graphs


class AsyncSaver:
    """Background thread for saving graphs to disk"""

    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.queue = queue.Queue(maxsize=4)  # Buffer up to 4 batches
        self.thread = threading.Thread(target=self._worker, daemon=True)
        self.running = True
        self.thread.start()

    def _worker(self):
        """Worker thread that saves graphs"""
        while self.running:
            try:
                item = self.queue.get(timeout=0.1)
                if item is None:  # Poison pill
                    break

                graphs, output_dir = item
                for graph in graphs:
                    output_path = output_dir / graph.filename
                    filename = graph.filename
                    del graph.filename
                    torch.save(graph, output_path)

                self.queue.task_done()
            except queue.Empty:
                continue

    def save_batch(self, graphs, output_dir):
        """Queue a batch of graphs for saving"""
        self.queue.put((graphs, output_dir))

    def wait(self):
        """Wait for all pending saves to complete"""
        self.queue.join()

    def shutdown(self):
        """Shutdown the saver thread"""
        self.running = False
        self.queue.put(None)  # Poison pill
        self.thread.join()


def collate_fn(batch):
    """Custom collate that returns list of graphs"""
    return batch


def run_graph_construction(model, hparams, config):
    """Build edges using learned embeddings + KNN with TRUE batching"""

    input_dir = Path(config['input_dir'])
    output_dir = Path(config['output_dir'])
    k_max = config['k_max']
    r_max = config.get('r_max', 0.15)
    r_max_geometric = config.get('r_max_geometric', None)
    same_layer_cut = config.get('same_layer_cut', 'geometry')
    segmented = config.get('segmented', True)
    one_in_one_out = config.get('one_in_one_out', False)
    device = config.get('device', 'cpu')
    datasets_config = config.get('datasets', ['trainset', 'valset', 'testset'])
    # Support both list format (process all) and dict format {name: max_count}
    if isinstance(datasets_config, dict):
        datasets = list(datasets_config.keys())
        dataset_limits = datasets_config
    else:
        datasets = datasets_config
        dataset_limits = {}
    node_features = hparams['node_features']
    node_scales = hparams.get('node_scales', None)

    # Optimization parameters
    num_workers = config.get('num_workers', 8)
    prefetch_factor = config.get('prefetch_factor', 4)
    batch_size = config.get('process_batch_size', 16)

    print(f"Input:  {input_dir}")
    print(f"Output: {output_dir}")
    print(f"Max neighbors (KNN): {k_max}")
    print(f"Max radius (latent space): {r_max}")
    print(f"Max geometric distance (3D, mm): {r_max_geometric if r_max_geometric is not None else 'disabled'}")
    print(f"Same-layer cut: {same_layer_cut}")
    print(f"Device: {device}")
    print(f"Node features: {node_features}")
    print(f"Num workers: {num_workers}")
    print(f"Prefetch factor: {prefetch_factor}")
    print(f"Process batch size: {batch_size}")
    print(f"Using TRUE batched inference + async saving")
    print()

    # Move model to GPU once
    model = model.to(device)
    model.eval()

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

        # Create dataset and dataloader
        dataset = GraphDataset(event_files, input_dataset_dir, node_features)

        dataloader = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=False,
            num_workers=num_workers,
            collate_fn=collate_fn,
            prefetch_factor=prefetch_factor if num_workers > 0 else None,
            persistent_workers=True if num_workers > 0 else False,
            pin_memory=True if device == 'cuda' else False
        )

        # Create async saver
        saver = AsyncSaver(output_dataset_dir)

        total_nodes = 0
        total_edges = 0
        total_true_edges = 0
        num_graphs_processed = 0

        # Process batches
        pbar = tqdm(dataloader, desc=f"  {dataset_name}", total=len(dataloader))

        try:
            for batch_graphs in pbar:
                # Process batch on GPU with TRUE batching
                processed_graphs = build_edges_latent_true_batch(
                    model, batch_graphs, node_features, k_max, r_max, r_max_geometric, same_layer_cut,
                    segmented=segmented, one_in_one_out=one_in_one_out, node_scales=node_scales,
                )

                # Update statistics (do this before async save)
                for graph in processed_graphs:
                    total_nodes += graph.num_nodes
                    total_edges += graph.edge_index.shape[1]
                    # Count unique undirected true edges: canonicalize (i,j) and (j,i) to
                    # (min,max) so bidirectional true labels aren't double-counted
                    true_mask = graph.edge_y.bool()
                    true_ei = graph.edge_index[:, true_mask]
                    if true_ei.shape[1] > 0:
                        canonical = torch.stack([true_ei.min(dim=0)[0], true_ei.max(dim=0)[0]])
                        num_unique = canonical.T.unique(dim=0).shape[0]
                    else:
                        num_unique = 0
                    total_true_edges += num_unique
                    num_graphs_processed += 1

                # Save graphs asynchronously
                saver.save_batch(processed_graphs, output_dataset_dir)

                # Update progress bar
                if num_graphs_processed > 0:
                    pbar.set_postfix({
                        'avg_nodes': f'{total_nodes/num_graphs_processed:.0f}',
                        'avg_edges': f'{total_edges/num_graphs_processed:.0f}',
                        'queue': saver.queue.qsize()
                    })

            # Wait for all saves to complete
            saver.wait()

        finally:
            saver.shutdown()

        # Print statistics
        print(f"  ✓ Total nodes: {total_nodes}")
        print(f"  ✓ Total edges: {total_edges}")
        if total_edges > 0:
            print(f"  ✓ True edges: {total_true_edges} ({100*total_true_edges/total_edges:.1f}%)")
        if total_nodes > 0:
            print(f"  ✓ Avg edges/node: {total_edges/total_nodes:.1f}")

    print("\n✓ Graph construction complete!\n")


def main():
    """Run graph construction with learned embeddings"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Build graphs using learned latent space ",
        epilog="Example: python build_latent_graphs_fast.py latent_builder_10epochs"
    )
    parser.add_argument(
        "model_name",
        type=str,
        nargs='?',
        default=None,
        help="Model name (without .ckpt extension) from saved_models/ directory"
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=8,
        help="Number of data loading workers (default: 8)"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=16,
        help="Number of graphs to process at once (default: 16)"
    )

    args = parser.parse_args()

    print("\n" + "="*80)
    print("GNN TRAINING PIPELINE - LEARNED GRAPH CONSTRUCTION ")
    print("="*80 + "\n")

    # Load configuration
    config = load_config(SCRIPT_DIR / 'acorn_configs' / 'graph_construction_latent.yaml')

    # Override with command line args
    if args.num_workers is not None:
        config['num_workers'] = args.num_workers
    if args.batch_size is not None:
        config['process_batch_size'] = args.batch_size

    print("="*80)
    print("CONFIGURATION")
    print("="*80)
    print(yaml.dump(config))
    print("="*80 + "\n")

    # Get checkpoint name
    checkpoint_name = args.model_name or config.get('checkpoint')
    if not checkpoint_name:
        print("ERROR: No model specified.")
        print("Usage: python build_latent_graphs_fast.py <model_name>")
        print("\nExample: python build_latent_graphs_fast.py latent_builder_10epochs")
        sys.exit(1)

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
    print(f"Ready for training with: python train_myGNN.py")
    print(f"  (Update config to use input_dir: {config['output_dir']})\n")


if __name__ == "__main__":
    main()
