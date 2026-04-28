"""Helper for saving plot data as JSON next to the plot file."""
import json
import numpy as np
from pathlib import Path


def _clean(obj):
    """Recursively convert numpy types and NaN → None for JSON serialization."""
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _clean(obj.tolist())
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        v = float(obj)
        return None if np.isnan(v) else v
    if isinstance(obj, float):
        return None if np.isnan(obj) else obj
    return obj


def save_plot_data_json(plot_path, plot_data):
    """Save plot_data dict as JSON next to plot_path (same stem, .json extension).

    NaN values are serialized as JSON null so notebooks can load the file with
    json.load() and reconstruct the exact same arrays.
    """
    json_path = Path(plot_path).with_suffix('.json')
    with open(json_path, 'w') as f:
        json.dump(_clean(plot_data), f, indent=2)
