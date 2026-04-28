"""
Simple script to evaluate model accuracy on test data.
Uses the model's data loading pipeline to ensure proper preprocessing.

Example usage:
    python test_my_gnn.py Low_pt_gnn_300_MIXED_1i1o_best_val_loss=0.0014 --num-events 2 --edge-cut 0.8 --validation
"""



import argparse
import sys
from pathlib import Path
import yaml
import torch
from sklearn.metrics import roc_auc_score

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / 'acorn'))
sys.path.insert(0, str(PIPELINE_ROOT))

from acorn.stages.edge_classifier.models.interaction_gnn import InteractionGNN


def find_checkpoint(model_name, base_dir=None):
    """Find checkpoint file in saved_models directory."""
    if base_dir is None:
        base_dir = PIPELINE_ROOT / "saved_models"
    
    checkpoint_path = base_dir / f"{model_name}.ckpt"
    if checkpoint_path.exists():
        return checkpoint_path
    
    raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}")


def evaluate_accuracy(model_name, config_file='acorn_configs/gnn_train.yaml', edge_cut=0.5, use_validation=False, num_events=None):
    """
    Evaluate model accuracy on test or validation data.
    
    Args:
        model_name: Name of the model checkpoint (e.g., 'epoch9_900')
        config_file: Path to config file
        edge_cut: Score threshold for classifying edges as positive (default: 0.5)
        use_validation: If True, evaluate on validation set; if False, evaluate on test set
    """
    dataset_name = "validation" if use_validation else "test"
    print("=" * 70)
    print(f"EVALUATING MODEL ACCURACY: {model_name}")
    print(f"Dataset: {dataset_name.upper()} set")
    print("=" * 70)
    print()

    # Load config - resolve relative paths from PIPELINE_ROOT
    config_path = Path(config_file)
    if not config_path.is_absolute():
        config_path = PIPELINE_ROOT / config_file

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Find checkpoint
    checkpoint_path = find_checkpoint(model_name)

    # Load checkpoint to get its hyperparameters (preserve architecture params)
    print("Loading model from checkpoint...")
    checkpoint = torch.load(str(checkpoint_path), map_location='cpu', weights_only=False)
    checkpoint_hparams = checkpoint.get('hyper_parameters', {})
    
    # Start with checkpoint hyperparameters (preserves architecture like 'hidden')
    # Then fill in any missing keys from config
    merged_hparams = {**checkpoint_hparams}
    for key, value in config.items():
        if key not in merged_hparams:
            print(f"Setting missing hyperparameter from config: {key} = {value}")
            merged_hparams[key] = value
    
    # Override paths and data-related params from current config
    data_dir = PIPELINE_ROOT / "data" / "graph_constructed_latent"
    merged_hparams['input_dir'] = str(data_dir)
    merged_hparams['stage_dir'] = str(data_dir)
    
    # Count actual files in test/validation set directory
    if use_validation:
        eval_set_dir = data_dir / "valset"
    else:
        eval_set_dir = data_dir / "testset"
    
    actual_files = list(eval_set_dir.rglob("*.pyg"))
    num_files = len(actual_files)
    if num_events is not None:
        num_files = min(num_events, num_files)
    print(f"Found {num_files} files in {eval_set_dir.name}/")
    
    # Set data_split to use all available files in the evaluation set
    # Format: [trainset, valset, testset]
    if use_validation:
        merged_hparams['data_split'] = [0, num_files, 0]  # Only load validation
    else:
        merged_hparams['data_split'] = [0, 0, num_files]  # Only load test
    
    merged_hparams['reprocess_classifier'] = True
    
    # Load model with merged hyperparameters
    model = InteractionGNN.load_from_checkpoint(
        str(checkpoint_path),
        map_location='cuda' if torch.cuda.is_available() else 'cpu',
        **merged_hparams
    )
    
    # Setup data (this will preprocess it)
    if use_validation:
        model.setup(stage='fit')  # 'fit' stage loads train/val/test sets with preprocessing
        data_loader = model.val_dataloader()
    else:
        model.setup(stage='test')  # Now with reprocess_classifier=True, this will also preprocess
        data_loader = model.test_dataloader()
    
    num_samples = len(data_loader.dataset) if hasattr(data_loader.dataset, '__len__') else 'unknown'
    print(f"Loaded {num_samples} {dataset_name} graphs")
    print()
    
    model.eval()
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model = model.to(device)

    
    # Evaluate
    print("=" * 70)
    print("RUNNING EVALUATION")
    print("=" * 70)
    
    # Metrics for ALL edges
    total_edges = 0
    all_tp = 0
    all_fp = 0
    all_tn = 0
    all_fn = 0
    
    with torch.no_grad():
        for batch_idx, batch in enumerate(data_loader):
            batch = batch.to(device)
            
            # Forward pass
            output = model(batch)
            scores = torch.sigmoid(output)
            predictions = (scores > edge_cut).cpu()
            
            # Get truth labels
            all_truth = batch.edge_y.bool().cpu()
            
            # Compute metrics for ALL edges
            tp = (predictions & all_truth).sum().item()
            fp = (predictions & ~all_truth).sum().item()
            tn = (~predictions & ~all_truth).sum().item()
            fn = (~predictions & all_truth).sum().item()
            
            all_tp += tp
            all_fp += fp
            all_tn += tn
            all_fn += fn
            total_edges += len(all_truth)
            
            if (batch_idx + 1) % 10 == 0:
                print(f"Processed {batch_idx + 1}/{len(data_loader)} batches...")
    
    print()
    print("=" * 70)
    print("RESULTS")
    print("=" * 70)
    print(f"Edge classification threshold: {edge_cut}")
    print()
    print(f"Total edges evaluated: {total_edges:,}")
    print()
    print("Classification breakdown:")
    print(f"  True Positives (TP):  {all_tp:>10,} - Correctly identified true edges")
    print(f"  False Positives (FP): {all_fp:>10,} - False edges incorrectly called true")
    print(f"  True Negatives (TN):  {all_tn:>10,} - Correctly identified false edges")
    print(f"  False Negatives (FN): {all_fn:>10,} - True edges incorrectly called false")
    print()
    
    accuracy = (all_tp + all_tn) / total_edges if total_edges > 0 else 0
    precision = all_tp / (all_tp + all_fp) if (all_tp + all_fp) > 0 else 0
    recall = all_tp / (all_tp + all_fn) if (all_tp + all_fn) > 0 else 0
    f1_score = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    print("Performance metrics:")
    print(f"  Accuracy:  {accuracy:.4f} ({accuracy*100:.2f}%)")
    print(f"  Precision: {precision:.4f} ({precision*100:.2f}%) - TP / (TP + FP)")
    print(f"  Recall:    {recall:.4f} ({recall*100:.2f}%) - TP / (TP + FN)")
    print(f"  F1 Score:  {f1_score:.4f}")
    print()
    print("=" * 70)


def main():
    parser = argparse.ArgumentParser(description="Evaluate model accuracy on test data")
    parser.add_argument(
        'model_name',
        type=str,
        help='Name of the model checkpoint (e.g., epoch9_900)'
    )
    parser.add_argument(
        '--config',
        type=str,
        default='acorn_configs/gnn_stage_(2)/gnn_train.yaml',
        help='Path to config file'
    )
    parser.add_argument(
        '--edge-cut',
        type=float,
        default=0.5,
        help='Score threshold for edge classification (default: 0.5)'
    )
    parser.add_argument(
        '--validation',
        action='store_true',
        help='Evaluate on validation set instead of test set'
    )
    parser.add_argument(
        '--num-events',
        type=int,
        default=None,
        help='Number of events to evaluate (default: all)'
    )

    args = parser.parse_args()
    evaluate_accuracy(args.model_name, args.config, args.edge_cut, args.validation, args.num_events)


if __name__ == '__main__':
    main()
