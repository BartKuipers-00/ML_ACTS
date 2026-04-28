"""
Evaluation and plotting utilities for track building algorithms.

Contains general-purpose functions for:
- Loading reconstruction and particle data from graphs
- Evaluating track quality metrics (efficiency, fake rate, clone rate)
- Plotting efficiency, clone rate, purity, and completeness vs pT/eta
- Saving evaluation results to disk

These utilities can be used with any track building algorithm (Connected Components,
custom loop builder, etc.) as long as graphs contain the required attributes:
- hit_track_labels (reconstructed track assignments)
- hit_particle_id (ground truth particle IDs)
- track_priority (optional: track priority labels if available)
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
import scipy.stats
import torch
import matplotlib.pyplot as plt
from tqdm import tqdm

# Import ACORN utilities
import sys
SCRIPT_DIR = Path(__file__).resolve().parent.parent
WORKSPACE_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(WORKSPACE_ROOT / 'acorn'))

from acorn.stages.track_building.utils import get_matching_df, calculate_matching_fraction
from acorn.stages.track_building.track_building_stage import make_result_summary, GraphDataset

try:
    from atlasify import atlasify
    HAS_ATLASIFY = True
except ImportError:
    HAS_ATLASIFY = False

from low_pt_custom_utils.plot_data_utils import save_plot_data_json

# PDG code → merged species name (charge conjugates merged)
PDG_TO_SPECIES = {
    13: 'Muon', -13: 'Muon',
    211: 'Pion', -211: 'Pion',
    11: 'Electron', -11: 'Electron',
    2212: 'Proton', -2212: 'Proton',
    321: 'Kaon', -321: 'Kaon',
}

SPECIES_COLORS = {
    'Muon': '#1f77b4',      # blue
    'Pion': '#ff7f0e',      # orange
    'Electron': '#2ca02c',  # green
    'Proton': '#d62728',    # red
    'Kaon': '#9467bd',      # purple
}

SPECIES_MARKERS = {
    'Muon': 'o',
    'Pion': 's',
    'Electron': '^',
    'Proton': 'D',
    'Kaon': 'v',
}


# ─── Graph Data Loading ──────────────────────────────────────────────────────

def safe_load_reconstruction_df(graph):
    """
    Load reconstruction DataFrame from a graph with track labels.

    Handles tensor conversion and adds priority information if available.

    Args:
        graph: PyG Data object with hit_track_labels and hit_particle_id

    Returns:
        DataFrame with columns: hit_id, track_id, particle_id, priority
    """
    if hasattr(graph, "hit_id"):
        hit_id = graph.hit_id
        if isinstance(hit_id, torch.Tensor):
            hit_id = hit_id.cpu().numpy()
    else:
        hit_id = torch.arange(graph.num_nodes).numpy()

    track_labels = graph.hit_track_labels
    if isinstance(track_labels, torch.Tensor):
        track_labels = track_labels.cpu().numpy()

    reco_df = pd.DataFrame({"hit_id": hit_id, "track_id": track_labels})

    # Add priority information if available
    if hasattr(graph, 'track_priority'):
        track_priority = graph.track_priority
        if isinstance(track_priority, torch.Tensor):
            track_priority = track_priority.cpu().numpy()
        # Map track_id to priority
        priority_map = {i: p for i, p in enumerate(track_priority)}
        reco_df["priority"] = reco_df["track_id"].map(priority_map).fillna(-1).astype(int)
    else:
        reco_df["priority"] = -1  # Unknown priority

    if hasattr(graph, 'hit_particle_id'):
        hit_pid = graph.hit_particle_id
        if isinstance(hit_pid, torch.Tensor):
            reco_df["particle_id"] = hit_pid.cpu().numpy()
        else:
            reco_df["particle_id"] = hit_pid
    else:
        try:
            node_id = graph.track_edges.reshape(-1)
            pids = graph.track_particle_id.repeat(2)
            if len(node_id) != len(pids):
                raise ValueError(
                    f"track_edges and track_particle_id length mismatch: "
                    f"{len(node_id)} vs {len(pids)}"
                )
            pid_df = pd.DataFrame({"hit_id": node_id, "particle_id": pids})
            pid_df.drop_duplicates(subset=["hit_id", "particle_id"], inplace=True)
            reco_df = reco_df.merge(pid_df, on="hit_id", how="outer")
        except (ValueError, AttributeError) as e:
            print(f"Warning: Could not load particle_id from track_edges: {e}")
            reco_df["particle_id"] = 0

    reco_df.fillna({"track_id": -1, "particle_id": 0}, inplace=True)
    return reco_df


def safe_load_particles_df(graph, sel_conf):
    """
    Build particles DataFrame from hit_particle_id and available attributes.

    Args:
        graph: PyG Data object with hit_particle_id and particle properties
        sel_conf: Selection configuration dict

    Returns:
        DataFrame with columns: particle_id, pt, eta
    """
    if 'hit_particle_id' in graph:
        unique_pids = torch.unique(graph.hit_particle_id)
        unique_pids = unique_pids[unique_pids > 0].numpy()
    else:
        raise ValueError("Graph must have hit_particle_id attribute")

    particles_data = {"particle_id": unique_pids}

    # Map pt to particles
    if 'track_edges' in graph and 'pt' in graph:
        edge_hits = graph.track_edges.reshape(-1)
        edge_pids = graph.hit_particle_id[edge_hits].numpy()
        edge_pt = graph.pt.repeat(2).numpy()
        pid_pt_map = {}
        for pid, pval in zip(edge_pids, edge_pt):
            if pid > 0 and pid not in pid_pt_map:
                pid_pt_map[pid] = pval
        particles_data["pt"] = [pid_pt_map.get(pid, 0.0) for pid in unique_pids]
    elif 'track_particle_pt' in graph and 'track_edges' in graph:
        # Map particle IDs to pT using track_edges (each edge has a pT value)
        edge_pids = graph.hit_particle_id[graph.track_edges[0, :]].numpy()
        edge_pt = graph.track_particle_pt.numpy()
        pid_pt_map = {}
        for pid, pval in zip(edge_pids, edge_pt):
            if pid > 0 and pid not in pid_pt_map:
                pid_pt_map[pid] = pval
        particles_data["pt"] = [pid_pt_map.get(pid, 0.0) for pid in unique_pids]
    else:
        particles_data["pt"] = [0.0] * len(unique_pids)

    # Map eta to particles
    if 'track_edges' in graph and 'eta_particle' in graph:
        edge_hits = graph.track_edges.reshape(-1)
        edge_pids = graph.hit_particle_id[edge_hits].numpy()
        edge_eta = graph.eta_particle.repeat(2).numpy()
        pid_eta_map = {}
        for pid, eval in zip(edge_pids, edge_eta):
            if pid > 0 and pid not in pid_eta_map:
                pid_eta_map[pid] = eval
        particles_data["eta"] = [pid_eta_map.get(pid, 0.0) for pid in unique_pids]
    elif 'track_particle_eta' in graph and 'track_edges' in graph:
        # Map particle IDs to eta using track_edges (each edge has an eta value)
        edge_pids = graph.hit_particle_id[graph.track_edges[0, :]].numpy()
        edge_eta = graph.track_particle_eta.numpy()
        pid_eta_map = {}
        for pid, eval in zip(edge_pids, edge_eta):
            if pid > 0 and pid not in pid_eta_map:
                pid_eta_map[pid] = eval
        particles_data["eta"] = [pid_eta_map.get(pid, 0.0) for pid in unique_pids]

    # Map particle_type (PDG code) to particles
    if 'track_particle_type' in graph and 'track_edges' in graph:
        edge_pids = graph.hit_particle_id[graph.track_edges[0, :]].numpy()
        edge_ptype = graph.track_particle_type.numpy()
        pid_ptype_map = {}
        for pid, ptype in zip(edge_pids, edge_ptype):
            if pid > 0 and pid not in pid_ptype_map:
                pid_ptype_map[pid] = int(ptype)
        particles_data["particle_type"] = [
            PDG_TO_SPECIES.get(pid_ptype_map.get(pid, 0), 'Other')
            for pid in unique_pids
        ]

    particles_df = pd.DataFrame(particles_data)
    particles_df = particles_df.drop_duplicates(subset=["particle_id"])
    return particles_df


def safe_evaluate_labelled_graph(graph, sel_conf, matching_fraction, matching_style, min_track_length):
    """
    Evaluate a single graph with track labels against truth.

    Args:
        graph: PyG Data object with hit_track_labels and hit_particle_id
        sel_conf: Selection configuration dict (fiducial cuts)
        matching_fraction: Minimum fraction of shared hits for matching
        matching_style: "ATLAS", "one_way", or "two_way"
        min_track_length: Minimum number of hits per track

    Returns:
        DataFrame with matching information (purity, completeness, is_matched, etc.)
    """
    reco_df = safe_load_reconstruction_df(graph)

    particles_sel_conf = {k: v for k, v in sel_conf.items() if k != 'nhits'}
    particles_df = safe_load_particles_df(graph, particles_sel_conf)

    matching_sel_conf = {}
    for k, v in sel_conf.items():
        if k == 'nhits':
            matching_sel_conf['n_true_hits'] = v
        else:
            matching_sel_conf[k] = v

    matching_df = get_matching_df(reco_df, particles_df, matching_sel_conf, min_track_length=min_track_length)

    # Merge priority information back (get_matching_df doesn't preserve extra columns)
    if 'priority' in reco_df.columns:
        priority_per_track = reco_df[['track_id', 'priority']].drop_duplicates(subset=['track_id'])
        matching_df = matching_df.merge(priority_per_track, on='track_id', how='left')
        matching_df['priority'] = matching_df['priority'].fillna(-1).astype(int)

    event_id = graph.event_id
    while type(event_id) == list:
        event_id = event_id[0]
    matching_df["event_id"] = int(event_id)

    matching_df = calculate_matching_fraction(matching_df)

    if matching_style == "ATLAS":
        matching_df["is_matched"] = matching_df["is_reconstructed"] = (
            matching_df.purity_reco >= matching_fraction
        )
    elif matching_style == "one_way":
        matching_df["is_matched"] = matching_df.purity_reco >= matching_fraction
        matching_df["is_reconstructed"] = matching_df.eff_true >= matching_fraction
    elif matching_style == "two_way":
        matching_df["is_matched"] = matching_df["is_reconstructed"] = (
            matching_df.purity_reco >= matching_fraction
        ) & (matching_df.eff_true >= matching_fraction)

    return matching_df


# ─── Evaluation ──────────────────────────────────────────────────────────────

def run_evaluation(dataset_name, eval_config):
    """
    Evaluate all events in a dataset and return results.

    Args:
        dataset_name: "trainset", "valset", or "testset"
        eval_config: Config dict with input_dir, matching_fraction, etc.

    Returns:
        evaluated_events: concatenated matching DataFrame for all events
        summary: dict with summary statistics (efficiency, fake_rate, etc.)
        summary_text: formatted text summary string
    """
    # Determine input directory
    input_dir = Path(eval_config.get('input_dir', eval_config.get('stage_dir')))
    if not input_dir.is_absolute():
        # Assume relative to low_pt_gnn_pipeline directory
        input_dir = SCRIPT_DIR / input_dir

    dataset_dir = input_dir / dataset_name
    if not dataset_dir.exists():
        raise FileNotFoundError(f"Dataset directory not found: {dataset_dir}")

    graph_files = list(dataset_dir.glob('*.pyg'))
    # Respect data_split if provided — only evaluate the number of events that were built
    data_split = eval_config.get("data_split", None)
    if data_split is not None:
        split_map = {"trainset": data_split[0], "valset": data_split[1], "testset": data_split[2]}
        max_events = split_map.get(dataset_name, len(graph_files))
        if max_events > 0 and max_events < len(graph_files):
            graph_files.sort()
            graph_files = graph_files[:max_events]
    num_events = len(graph_files)
    if num_events == 0:
        raise ValueError(f"No graph files found in {dataset_dir}")

    print(f"Found {num_events} events in {dataset_name}")
    print(f"Matching fraction: {eval_config.get('matching_fraction', 0.5)}")
    print(f"Matching style: {eval_config.get('matching_style', 'ATLAS')}")
    print()

    dataset = GraphDataset(
        input_dir=str(input_dir),
        data_name=dataset_name,
        num_events=num_events,
        stage="test",
        hparams=eval_config
    )

    evaluated_events = []
    times = []

    for event in tqdm(dataset, desc=f"Evaluating {dataset_name}"):
        try:
            matching_df = safe_evaluate_labelled_graph(
                event,
                sel_conf=eval_config.get("target_tracks", {}),
                matching_fraction=eval_config.get("matching_fraction", 0.5),
                matching_style=eval_config.get("matching_style", "ATLAS"),
                min_track_length=eval_config.get("min_track_length", 1),
            )
            evaluated_events.append(matching_df)
            if hasattr(event, 'time_taken'):
                times.append(event.time_taken)
        except Exception as e:
            event_id = getattr(event, 'event_id', 'unknown')
            if isinstance(event_id, list):
                event_id = event_id[0] if len(event_id) > 0 else 'unknown'
            print(f"\nError evaluating event {event_id}: {e}")
            import traceback
            traceback.print_exc()
            raise

    evaluated_events = pd.concat(evaluated_events, ignore_index=True)

    if times:
        times = np.array(times)
        time_avg, time_std = np.mean(times), np.std(times)
    else:
        time_avg, time_std = 0.0, 0.0

    # Calculate metrics
    particles = evaluated_events[evaluated_events["is_reconstructable"]]
    reconstructed_particles = particles[particles["is_reconstructed"] & particles["is_matchable"]]
    tracks = evaluated_events[evaluated_events["is_matchable"]]
    matched_tracks = tracks[tracks["is_matched"]]

    n_particles = len(particles.drop_duplicates(subset=["event_id", "particle_id"]))
    n_reconstructed = len(reconstructed_particles.drop_duplicates(subset=["event_id", "particle_id"]))
    n_tracks = len(tracks.drop_duplicates(subset=["event_id", "track_id"]))
    n_matched = len(matched_tracks.drop_duplicates(subset=["event_id", "track_id"]))
    n_duplicates = len(reconstructed_particles) - n_reconstructed

    efficiency = n_reconstructed / n_particles if n_particles > 0 else 0.0
    fake_rate = 1 - (n_matched / n_tracks) if n_tracks > 0 else 0.0
    clone_rate = n_duplicates / n_reconstructed if n_reconstructed > 0 else 0.0

    # Calculate efficiency per priority
    priority_stats = {}
    if 'priority' in evaluated_events.columns and evaluated_events['priority'].max() >= 0:
        tracks_with_priority = evaluated_events[evaluated_events["is_matchable"]].copy()
        unique_tracks = tracks_with_priority.drop_duplicates(subset=["event_id", "track_id"])

        for priority in sorted(unique_tracks['priority'].unique()):
            if priority < 0:  # Skip invalid priorities
                continue

            priority_tracks = unique_tracks[unique_tracks['priority'] == priority]
            matched_priority_tracks = priority_tracks[priority_tracks['is_matched']]

            n_tracks_p = len(priority_tracks)
            n_matched_p = len(matched_priority_tracks)
            efficiency_p = n_matched_p / n_tracks_p if n_tracks_p > 0 else 0.0

            priority_name = {0: 'Simple', 1: 'P1(Straight)', 2: 'P2(Loop-Hi)', 3: 'P3(Loop-Lo)', 4: 'P4(Leftover)'}.get(priority, f'P{priority}')
            priority_stats[priority_name] = {
                'n_tracks': int(n_tracks_p),
                'n_matched': int(n_matched_p),
                'efficiency': float(efficiency_p),
            }

    summary = {
        'dataset': dataset_name,
        'n_events': num_events,
        'n_particles': int(n_particles),
        'n_reconstructed_particles': int(n_reconstructed),
        'n_tracks': int(n_tracks),
        'n_matched_tracks': int(n_matched),
        'n_duplicate_particles': int(n_duplicates),
        'efficiency': float(efficiency),
        'fake_rate': float(fake_rate),
        'clone_rate': float(clone_rate),
        'avg_purity': float(evaluated_events['purity_reco'].mean() if 'purity_reco' in evaluated_events.columns else 0.0),
        'avg_completeness': float(evaluated_events['eff_true'].mean() if 'eff_true' in evaluated_events.columns else 0.0),
        'time_avg_seconds': float(time_avg),
        'time_std_seconds': float(time_std),
        'matching_fraction': eval_config.get('matching_fraction', 0.5),
        'matching_style': eval_config.get('matching_style', 'ATLAS'),
        'priority_stats': priority_stats,  # Add priority statistics
    }

    summary_text = make_result_summary(
        n_reconstructed, n_particles, n_matched, n_tracks, n_duplicates,
        efficiency, fake_rate, clone_rate, time_avg, time_std,
    )

    # Append priority statistics to summary text
    if priority_stats:
        summary_text += "\n\nEfficiency by Priority:"
        summary_text += "\n" + "=" * 60
        for priority_name in ['Simple', 'P1(Straight)', 'P2(Loop-Hi)', 'P3(Loop-Lo)', 'P4(Leftover)']:
            if priority_name in priority_stats:
                stats = priority_stats[priority_name]
                pct = 100 * stats['efficiency']
                summary_text += f"\n  {priority_name:15s}  {stats['n_matched']:6d}/{stats['n_tracks']:<6d} tracks matched  ({pct:5.1f}%)"
        summary_text += "\n" + "=" * 60

    return evaluated_events, summary, summary_text


def save_evaluation_results(evaluated_events, summary, summary_text, dataset_name, output_dir, debug_files=True):
    """
    Save evaluation results to disk.

    Args:
        evaluated_events: DataFrame with matching information
        summary: dict with summary statistics
        summary_text: formatted text summary
        dataset_name: "trainset", "valset", or "testset"
        output_dir: Path to output directory
        debug_files: if False, skip large per-event DataFrames (matching_df, particles)
                     and only write the compact summary JSON and text file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Large per-event DataFrames — only written in debug mode
    if debug_files:
        evaluated_events.to_csv(output_dir / f"matching_df_{dataset_name}.csv", index=False)

        particles = evaluated_events[evaluated_events["is_reconstructable"]]
        particles_plot = particles.copy()
        particles_plot["is_reconstructed"] = particles_plot.groupby(
            ["event_id", "particle_id"]
        )["is_reconstructed"].transform("any")
        particles_plot = particles_plot.drop_duplicates(subset=["event_id", "particle_id"])
        particles_plot.to_csv(output_dir / f"particles_{dataset_name}.csv", index=False)

    # Summary JSON — always written (tiny, contains the key metrics)
    with open(output_dir / f"summary_{dataset_name}.json", 'w') as f:
        json.dump(summary, f, indent=2)

    # Summary text — always written
    with open(output_dir / f"results_summary_{summary.get('matching_style', 'ATLAS')}_{dataset_name}.txt", 'w') as f:
        f.write(summary_text)

    print(f"Evaluation results saved to: {output_dir}")


# ─── Plotting ────────────────────────────────────────────────────────────────

def _clopper_pearson_errors(passed_arr, total_arr, level=0.68):
    """Clopper-Pearson confidence interval errors for arrays of (passed, total).
    Returns (lo_err, hi_err) arrays — distances from the central value."""
    alpha = (1 - level) / 2
    lo = np.where(passed_arr > 0,
                  scipy.stats.beta.ppf(alpha, passed_arr, total_arr - passed_arr + 1), 0.0)
    hi = np.where(passed_arr < total_arr,
                  scipy.stats.beta.ppf(1 - alpha, passed_arr + 1, total_arr - passed_arr), 1.0)
    eff = np.where(total_arr > 0, passed_arr / total_arr, np.nan)
    return np.maximum(0, eff - lo), np.maximum(0, hi - eff)


def _get_bins(var, varconf):
    """Determine histogram bins from variable config."""
    n_bins = varconf.get('n_bins', 20)
    if 'x_bins' in varconf:
        return np.array(varconf['x_bins'])
    elif 'x_lim' in varconf:
        return np.linspace(varconf['x_lim'][0], varconf['x_lim'][1], n_bins)
    else:
        if var == 'pt':
            return np.linspace(0.05, 1.0, n_bins)
        else:
            return np.linspace(-4, 4, n_bins)


def _get_species_list(df):
    """Get sorted list of species present in the DataFrame."""
    if 'particle_type' not in df.columns:
        return []
    return sorted(df['particle_type'].dropna().unique())


def plot_efficiency_vs_variable(particles_df, var, varconf, output_path, summary=None, font_sizes=None):
    """Plot efficiency vs pT or eta, with one curve per particle species."""
    if font_sizes is None:
        font_sizes = {'axis_label': 18, 'tick_label': 14, 'legend': 14, 'title': 18}

    x_all = particles_df[var].values
    if 'x_scale' in varconf:
        x_all = x_all * float(varconf['x_scale'])

    reconstructable = particles_df['is_reconstructable']
    reconstructed = particles_df['is_reconstructable'] & particles_df['is_reconstructed']

    if reconstructable.sum() == 0:
        print(f"Warning: No reconstructable particles with {var} data. Skipping.")
        return

    x_bins = _get_bins(var, varconf)
    xvals = (x_bins[1:] + x_bins[:-1]) / 2
    xerrs = (x_bins[1:] - x_bins[:-1]) / 2

    fig, ax = plt.subplots(figsize=(8, 6))

    species_list = _get_species_list(particles_df)
    if not species_list:
        species_list = ['All']

    title = f"Track Efficiency vs {var.upper()}"
    if summary and 'efficiency' in summary:
        title += f" (Overall: {summary['efficiency']:.3f})"

    plot_data = {
        "plot_type": "efficiency",
        "variable": var,
        "xlabel": varconf.get('x_label', var),
        "ylabel": "Efficiency",
        "ylim": varconf.get('y_lim', [0, 1.1]),
        "title": title,
        "series": [],
    }

    for species in species_list:
        if species == 'All':
            mask = np.ones(len(particles_df), dtype=bool)
        else:
            mask = (particles_df['particle_type'] == species).values

        true_x = x_all[mask & reconstructable.values]
        reco_x = x_all[mask & reconstructed.values]

        if len(true_x) == 0:
            continue

        true_vals, _ = np.histogram(true_x, bins=x_bins)
        reco_vals, _ = np.histogram(reco_x, bins=x_bins)

        with np.errstate(divide='ignore', invalid='ignore'):
            eff = np.where(true_vals > 0, reco_vals / true_vals, np.nan)
        lo_err, hi_err = _clopper_pearson_errors(reco_vals, true_vals)
        err = [lo_err, hi_err]

        n_reco = int(reco_vals.sum())
        n_true = int(true_vals.sum())
        label = f'{species} ({n_reco}/{n_true})'
        color = SPECIES_COLORS.get(species, 'black')
        marker = SPECIES_MARKERS.get(species, 'o')

        ax.errorbar(xvals, eff, xerr=xerrs, yerr=err, fmt=marker, color=color,
                    label=label, capsize=3, capthick=1.5, markersize=5)
        plot_data["series"].append({
            "label": label, "color": color, "marker": marker,
            "x": xvals, "xerr": xerrs,
            "y": eff, "yerr_lo": lo_err, "yerr_hi": hi_err,
        })

    ax.set_xlabel(varconf.get('x_label', var), fontsize=font_sizes['axis_label'])
    ax.set_ylabel('Efficiency', fontsize=font_sizes['axis_label'])
    ax.set_ylim(varconf.get('y_lim', [0, 1.1]))
    ax.tick_params(axis='both', labelsize=font_sizes['tick_label'])
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=font_sizes['legend'])
    if summary:
        ax.set_title(title, fontsize=font_sizes['title'])
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    save_plot_data_json(output_path, plot_data)
    plt.close()
    print(f"Saved: {output_path}")


def plot_clone_rate_vs_variable(matching_df, var, varconf, output_path, summary=None, font_sizes=None):
    """Plot clone rate vs pT or eta, with one curve per particle species."""
    if font_sizes is None:
        font_sizes = {'axis_label': 18, 'tick_label': 14, 'legend': 14, 'title': 18}

    particles = matching_df[matching_df['is_reconstructable'] & matching_df['is_reconstructed']].copy()
    particle_track_counts = particles.groupby(['event_id', 'particle_id']).size().reset_index(name='n_tracks')
    particle_track_counts['is_cloned'] = particle_track_counts['n_tracks'] > 1

    merge_cols = ['event_id', 'particle_id', var]
    if 'particle_type' in particles.columns:
        merge_cols.append('particle_type')
    particles_unique = particles.drop_duplicates(subset=['event_id', 'particle_id'])
    particle_track_counts = particle_track_counts.merge(
        particles_unique[merge_cols], on=['event_id', 'particle_id']
    )

    x_all = particle_track_counts[var].values
    if 'x_scale' in varconf:
        x_all = x_all * float(varconf['x_scale'])

    x_bins = _get_bins(var, varconf)
    xvals = (x_bins[1:] + x_bins[:-1]) / 2
    xerrs = (x_bins[1:] - x_bins[:-1]) / 2

    fig, ax = plt.subplots(figsize=(8, 6))

    species_list = _get_species_list(particle_track_counts)
    if not species_list:
        species_list = ['All']

    title = f"Clone Rate vs {var.upper()}"
    if summary and 'clone_rate' in summary:
        title += f" (Overall: {summary['clone_rate']:.3f})"

    plot_data = {
        "plot_type": "clone_rate",
        "variable": var,
        "xlabel": varconf.get('x_label', var),
        "ylabel": "Clone Rate",
        "ylim": [0, 1.1],
        "title": title,
        "series": [],
    }

    for species in species_list:
        if species == 'All':
            mask = np.ones(len(particle_track_counts), dtype=bool)
        else:
            mask = (particle_track_counts['particle_type'] == species).values

        subset = particle_track_counts[mask]
        sx = x_all[mask]
        cloned_x = sx[subset['is_cloned'].values]

        if len(sx) == 0:
            continue

        all_vals, _ = np.histogram(sx, bins=x_bins)
        cloned_vals, _ = np.histogram(cloned_x, bins=x_bins)

        with np.errstate(divide='ignore', invalid='ignore'):
            clone_rate = np.where(all_vals > 0, cloned_vals / all_vals, np.nan)
        lo_err, hi_err = _clopper_pearson_errors(cloned_vals, all_vals)
        err = [lo_err, hi_err]

        color = SPECIES_COLORS.get(species, 'orange')
        marker = SPECIES_MARKERS.get(species, 'o')
        ax.errorbar(xvals, clone_rate, xerr=xerrs, yerr=err, fmt=marker, color=color,
                    label=species, capsize=3, capthick=1.5, markersize=5)
        plot_data["series"].append({
            "label": species, "color": color, "marker": marker,
            "x": xvals, "xerr": xerrs,
            "y": clone_rate, "yerr_lo": lo_err, "yerr_hi": hi_err,
        })

    ax.set_xlabel(varconf.get('x_label', var), fontsize=font_sizes['axis_label'])
    ax.set_ylabel('Clone Rate', fontsize=font_sizes['axis_label'])
    ax.set_ylim([0, 1.1])
    ax.tick_params(axis='both', labelsize=font_sizes['tick_label'])
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=font_sizes['legend'])
    ax.set_title(title, fontsize=font_sizes['title'])
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    save_plot_data_json(output_path, plot_data)
    plt.close()
    print(f"Saved: {output_path}")


def _plot_binned_mean_vs_variable(matching_df, var, varconf, value_col, output_path,
                                   y_label, title_prefix, summary_key, summary=None, font_sizes=None):
    """Plot mean of a value column vs pT or eta, with one curve per particle species.

    Used for purity and completeness plots.
    """
    if font_sizes is None:
        font_sizes = {'axis_label': 18, 'tick_label': 14, 'legend': 14, 'title': 18}

    matched = matching_df[matching_df['is_matched'] & matching_df['is_matchable']].copy()
    if len(matched) == 0:
        print(f"Warning: No matched tracks for {title_prefix.lower()} plot. Skipping.")
        return

    x_all = matched[var].values
    if 'x_scale' in varconf:
        x_all = x_all * float(varconf['x_scale'])
    values = matched[value_col].values

    x_bins = _get_bins(var, varconf)
    xvals = (x_bins[1:] + x_bins[:-1]) / 2
    xerrs = (x_bins[1:] - x_bins[:-1]) / 2
    n_bins = len(x_bins) - 1

    fig, ax = plt.subplots(figsize=(8, 6))

    species_list = _get_species_list(matched)
    if not species_list:
        species_list = ['All']

    title = f"{title_prefix} vs {var.upper()}"
    if summary and summary_key in summary:
        title += f" (Overall: {summary[summary_key]:.3f})"

    plot_data = {
        "plot_type": title_prefix.lower().replace(' ', '_'),
        "variable": var,
        "xlabel": varconf.get('x_label', var),
        "ylabel": y_label,
        "ylim": [0, 1.1],
        "title": title,
        "series": [],
    }

    for species in species_list:
        if species == 'All':
            mask = np.ones(len(matched), dtype=bool)
        else:
            mask = (matched['particle_type'] == species).values

        sx = x_all[mask]
        sv = values[mask]
        if len(sx) == 0:
            continue

        bin_indices = np.digitize(sx, x_bins) - 1
        mean_vals = np.full(n_bins, np.nan)
        err_vals = np.full(n_bins, np.nan)

        for i in range(n_bins):
            bmask = bin_indices == i
            if np.sum(bmask) > 0:
                bv = sv[bmask]
                mean_vals[i] = np.mean(bv)
                err_vals[i] = np.std(bv) / np.sqrt(len(bv))

        color = SPECIES_COLORS.get(species, 'blue')
        marker = SPECIES_MARKERS.get(species, 'o')
        ax.errorbar(xvals, mean_vals, xerr=xerrs, yerr=err_vals, fmt=marker, color=color,
                    label=species, capsize=3, capthick=1.5, markersize=5)
        plot_data["series"].append({
            "label": species, "color": color, "marker": marker,
            "x": xvals, "xerr": xerrs,
            "y": mean_vals, "yerr_lo": err_vals, "yerr_hi": err_vals,
        })

    ax.set_xlabel(varconf.get('x_label', var), fontsize=font_sizes['axis_label'])
    ax.set_ylabel(y_label, fontsize=font_sizes['axis_label'])
    ax.set_ylim([0, 1.1])
    ax.tick_params(axis='both', labelsize=font_sizes['tick_label'])
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=font_sizes['legend'])
    ax.set_title(title, fontsize=font_sizes['title'])
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    save_plot_data_json(output_path, plot_data)
    plt.close()
    print(f"Saved: {output_path}")


def plot_purity_vs_variable(matching_df, var, varconf, output_path, summary=None, font_sizes=None):
    """Plot mean purity of matched tracks vs pT or eta, per species."""
    _plot_binned_mean_vs_variable(
        matching_df, var, varconf, 'purity_reco', output_path,
        y_label='Purity (fraction of correct hits)',
        title_prefix='Track Purity', summary_key='avg_purity',
        summary=summary, font_sizes=font_sizes,
    )


def plot_completeness_vs_variable(matching_df, var, varconf, output_path, summary=None, font_sizes=None):
    """Plot mean completeness of matched tracks vs pT or eta, per species."""
    _plot_binned_mean_vs_variable(
        matching_df, var, varconf, 'eff_true', output_path,
        y_label='Completeness (fraction of true hits found)',
        title_prefix='Track Completeness', summary_key='avg_completeness',
        summary=summary, font_sizes=font_sizes,
    )


def run_plotting(evaluated_events, summary, dataset_name, output_dir, plot_config=None):
    """
    Create all metric plots.

    Args:
        evaluated_events: DataFrame with matching information
        summary: dict with summary statistics
        dataset_name: "trainset", "valset", or "testset"
        output_dir: Path to output directory
        plot_config: Optional dict with plot configuration (variables, limits, etc.)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Default plot config if not provided
    if plot_config is None:
        plot_config = {}

    # Load particles for efficiency plots
    particles = evaluated_events[evaluated_events["is_reconstructable"]].copy()
    # A particle is reconstructed if ANY of its matching rows has is_reconstructed=True
    particles["is_reconstructed"] = particles.groupby(
        ["event_id", "particle_id"]
    )["is_reconstructed"].transform("any")
    particles_df = particles.drop_duplicates(subset=["event_id", "particle_id"])

    # Default plot config (fallback if not in config file)
    default_plot_config = {
        'pt': {'x_label': '$p_T$ [GeV]', 'x_lim': [0.05, 1.0], 'y_lim': [0, 1.1]},
        'eta': {'x_label': '$\\eta$', 'x_lim': [-4, 4], 'y_lim': [0, 1.1]},
    }

    # Extract plot config from eval config (under 'tracking_efficiency' -> 'variables')
    config_vars = plot_config.get('tracking_efficiency', {}).get('variables', {})

    # Extract font sizes from config (under 'tracking_efficiency' -> 'font_sizes')
    font_sizes = plot_config.get('tracking_efficiency', {}).get('font_sizes', {
        'axis_label': 18, 'tick_label': 14, 'legend': 14, 'title': 18
    })

    for var in ['pt', 'eta']:
        # Use config from file if available, otherwise use default
        varconf = config_vars.get(var, default_plot_config[var])

        eff_path = output_dir / f"efficiency_vs_{var}_{dataset_name}.png"
        plot_efficiency_vs_variable(particles_df, var, varconf, eff_path, summary, font_sizes)

        clone_path = output_dir / f"clone_rate_vs_{var}_{dataset_name}.png"
        plot_clone_rate_vs_variable(evaluated_events, var, varconf, clone_path, summary, font_sizes)

        purity_path = output_dir / f"purity_vs_{var}_{dataset_name}.png"
        plot_purity_vs_variable(evaluated_events, var, varconf, purity_path, summary, font_sizes)

        completeness_path = output_dir / f"completeness_vs_{var}_{dataset_name}.png"
        plot_completeness_vs_variable(evaluated_events, var, varconf, completeness_path, summary, font_sizes)

    print(f"Plots saved to: {output_dir}")
