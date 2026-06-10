#!/usr/bin/env python3
"""
create_visual_evaluated_tracks.py - Visualize evaluated track building results

Loads a built-track graph from data/track_building/<dataset>/, evaluates it
against ground truth, and produces an interactive 3D visualization with:

  Green:  Correctly matched tracks — legend shows particle species, ID,
          purity and completeness for each track.
  Red:    Fake tracks — legend shows majority-particle purity and completeness
          even though the match is below the matching threshold.
  Orange: Missed particles — ground-truth spacepoints for particles that were
          not reconstructed (may spatially overlap with other tracks).

Default matching: two_way at 70 % threshold.

Usage:
  python create_visual_evaluated_tracks.py <dataset> <index>
  python create_visual_evaluated_tracks.py testset 1
  python create_visual_evaluated_tracks.py valset 3 --matching-fraction 0.5
  python create_visual_evaluated_tracks.py testset 2 --output my_plot.html
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import torch

try:
    import plotly.graph_objects as go
except ImportError:
    print("Error: plotly is required. Install with: pip install plotly")
    sys.exit(1)

# ─── Path Setup ───────────────────────────────────────────────────────────────

SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / "acorn"))
sys.path.insert(0, str(PIPELINE_ROOT))

from low_pt_custom_utils.track_evaluation_utils import safe_evaluate_labelled_graph
from visual_utils import (
    validate_dataset_name,
    extract_node_coordinates,
    get_standard_scene_layout,
    build_hit_particle_type_map,
    pdg_to_particle_name,
)

# ─── Constants ────────────────────────────────────────────────────────────────

# Evaluation colours
COLOR_MATCHED = "green"
COLOR_FAKE = "red"
COLOR_MISSED = "orange"
COLOR_NOISE = "rgba(150,150,150,0.25)"


# ─── Graph Loading ────────────────────────────────────────────────────────────

def load_track_building_graph(dataset_name, index):
    """Load a graph from data/track_building/<dataset_name>/."""
    validate_dataset_name(dataset_name)

    dataset_dir = PIPELINE_ROOT / "data" / "track_building" / dataset_name
    if not dataset_dir.exists():
        raise FileNotFoundError(f"Dataset directory not found: {dataset_dir}")

    graph_files = sorted(dataset_dir.glob("*.pyg"))
    if not graph_files:
        raise ValueError(f"No graph files found in {dataset_dir}")

    if index < 1 or index > len(graph_files):
        raise ValueError(f"Index {index} out of range. Available: 1-{len(graph_files)}")

    graph_path = graph_files[index - 1]
    print(f"Loading {dataset_name} graph #{index}: {graph_path.name}")

    graph = torch.load(graph_path, map_location="cpu", weights_only=False)

    if not hasattr(graph, "hit_track_labels"):
        raise ValueError(
            f"Graph missing 'hit_track_labels'. Re-run track building:\n"
            f"  python helix_segmentmatcher_walkthrough.py {dataset_name}"
        )
    if not hasattr(graph, "hit_particle_id"):
        raise ValueError("Graph missing 'hit_particle_id'.")

    return graph, graph_path


# ─── Evaluation ───────────────────────────────────────────────────────────────

def evaluate_graph(graph, matching_fraction, matching_style,
                   min_track_length, min_particle_hits):
    """Run per-event evaluation and return matching DataFrame."""
    # nhits: minimum number of true hits for a particle to be reconstructable.
    # Without this, 1-hit noise particles inflate the 'missed' count.
    sel_conf = {"nhits": min_particle_hits} if min_particle_hits > 0 else {}
    matching_df = safe_evaluate_labelled_graph(
        graph,
        sel_conf=sel_conf,
        matching_fraction=matching_fraction,
        matching_style=matching_style,
        min_track_length=min_track_length,
    )
    return matching_df


# ─── Classification Helpers ───────────────────────────────────────────────────

def classify_tracks_from_matching(matching_df, hit_track_labels):
    """
    Return dicts with per-track and per-missed-particle info.

    matched_tracks: {track_id -> {particle_id, purity, completeness}}
    fake_tracks:    {track_id -> {particle_id, purity, completeness}}
    missed_pids:    set of particle_ids that are reconstructable but not reconstructed
    """
    matched_tracks = {}
    fake_tracks = {}

    # Ground-truth of which track_ids are matched comes from matching_df rows
    # where is_matched=True.  We do NOT filter by is_matchable here because
    # with an empty sel_conf it can be False for all rows, hiding all fakes.
    matched_in_df = set(
        matching_df[matching_df["is_matched"] == True]["track_id"].unique()
    )

    # Build per-track best-particle lookup from matching_df
    track_best = {}  # track_id -> {particle_id, purity, completeness}
    reco_rows = matching_df[matching_df["track_id"] >= 0].copy()
    for track_id, grp in reco_rows.groupby("track_id"):
        best_idx = grp["purity_reco"].idxmax()
        best = grp.loc[best_idx]
        track_best[int(track_id)] = {
            "particle_id": int(best["particle_id"]),
            "purity": float(best["purity_reco"]),
            "completeness": float(best["eff_true"]),
        }

    # Only consider tracks that pass is_matchable (i.e. meet min_track_length).
    # Exclude track_id == -1 (isolated hits, not a real cluster).
    matchable_track_ids = set(
        matching_df[(matching_df["is_matchable"] == True) & (matching_df["track_id"] >= 0)]["track_id"].unique()
    )
    for track_id in matchable_track_ids:
        info = track_best.get(track_id, {"particle_id": 0, "purity": 0.0, "completeness": 0.0})
        if track_id in matched_in_df:
            matched_tracks[track_id] = info
        else:
            fake_tracks[track_id] = info

    # Missed particles: reconstructable but not reconstructed
    reconstructable = matching_df[matching_df["is_reconstructable"] == True].copy()
    reconstructed_pids = set(
        reconstructable[reconstructable["is_reconstructed"] == True]["particle_id"].unique()
    )
    all_reconstructable_pids = set(reconstructable["particle_id"].unique())
    missed_pids = all_reconstructable_pids - reconstructed_pids

    return matched_tracks, fake_tracks, missed_pids


# ─── Visualization ────────────────────────────────────────────────────────────

def create_evaluated_visualization(graph, matching_df, dataset_name, index):
    """
    Build Plotly Figure with green/red/orange colour coding.

    Returns:
        plotly.graph_objects.Figure
    """
    # Node coordinates
    coords = extract_node_coordinates(graph, return_cylindrical=True, return_cartesian=True)
    x, y, z = coords["x"], coords["y"], coords["z_cart"]
    r, phi, z_cyl = coords["r"], coords["phi"], coords["z"]
    num_nodes = len(x)

    # Particle type lookup
    particle_type_map = build_hit_particle_type_map(graph)

    # Hit arrays
    hit_track_labels = np.array(graph.hit_track_labels.tolist())
    hit_particle_ids = np.array(graph.hit_particle_id.tolist())

    # Classify
    matched_tracks, fake_tracks, missed_pids = classify_tracks_from_matching(
        matching_df, hit_track_labels
    )

    # Hit-level time for GT trajectory ordering
    hit_t = np.array(graph.hit_t.tolist()) if hasattr(graph, "hit_t") else None

    fig = go.Figure()

    # ── Matched tracks (green) ────────────────────────────────────────────────
    for track_id, info in matched_tracks.items():
        mask = hit_track_labels == track_id
        if not mask.any():
            continue

        pid = info["particle_id"]
        pdg_code = particle_type_map.get(pid, 0)
        species = pdg_to_particle_name(pdg_code) if pdg_code else "?"
        purity = info["purity"]
        completeness = info["completeness"]

        legend_label = (
            f"[OK] Track {track_id} | {species} PID {pid} | "
            f"pur={purity:.2f} cmp={completeness:.2f}"
        )

        node_indices = np.where(mask)[0]
        hover = [
            f"Track {track_id} (matched)<br>"
            f"Particle: {pid} ({species})<br>"
            f"Purity: {purity:.3f}<br>"
            f"Completeness: {completeness:.3f}<br>"
            f"r: {r[i]:.1f} mm  phi: {phi[i]:.3f}  z: {z_cyl[i]:.1f} mm"
            for i in node_indices
        ]

        fig.add_trace(go.Scatter3d(
            x=x[mask], y=y[mask], z=z[mask],
            mode="markers",
            marker=dict(size=4, color=COLOR_MATCHED, opacity=0.85,
                        line=dict(width=0.5, color="darkgreen")),
            name=legend_label,
            legendgroup=legend_label,
            text=hover,
            hoverinfo="text",
            showlegend=True,
        ))

        # GT trajectory edges — same legendgroup so they toggle together
        _add_gt_trajectory_edges(
            fig, pid, hit_particle_ids, hit_t, x, y, z,
            color=COLOR_MATCHED, opacity=0.5, legendgroup=legend_label,
        )

    # ── Fake tracks (red) ─────────────────────────────────────────────────────
    for track_id, info in fake_tracks.items():
        mask = hit_track_labels == track_id
        if not mask.any():
            continue

        pid = info["particle_id"]
        pdg_code = particle_type_map.get(pid, 0)
        species = pdg_to_particle_name(pdg_code) if pdg_code else "?"
        purity = info["purity"]
        completeness = info["completeness"]

        legend_label = (
            f"[FAKE] Track {track_id} | best PID {pid} ({species}) | "
            f"pur={purity:.2f} cmp={completeness:.2f}"
        )

        node_indices = np.where(mask)[0]
        hover = [
            f"Track {track_id} (FAKE)<br>"
            f"Best particle: {pid} ({species})<br>"
            f"Best purity: {purity:.3f}<br>"
            f"Best completeness: {completeness:.3f}<br>"
            f"r: {r[i]:.1f} mm  phi: {phi[i]:.3f}  z: {z_cyl[i]:.1f} mm"
            for i in node_indices
        ]

        fig.add_trace(go.Scatter3d(
            x=x[mask], y=y[mask], z=z[mask],
            mode="markers",
            marker=dict(size=4, color=COLOR_FAKE, opacity=0.75,
                        line=dict(width=0.5, color="darkred")),
            name=legend_label,
            legendgroup=legend_label,
            text=hover,
            hoverinfo="text",
            showlegend=True,
        ))

    # ── Missed particles (orange) ─────────────────────────────────────────────
    # Look up pt per particle from matching_df for hover
    pid_pt_map = {}
    if "pt" in matching_df.columns:
        for pid_val, grp in matching_df[matching_df["is_reconstructable"]].groupby("particle_id"):
            pid_pt_map[int(pid_val)] = float(grp["pt"].iloc[0])

    for pid in sorted(missed_pids):
        mask = hit_particle_ids == pid
        if not mask.any():
            continue

        pdg_code = particle_type_map.get(pid, 0)
        species = pdg_to_particle_name(pdg_code) if pdg_code else "?"
        pt_str = f"{pid_pt_map[pid]*1000:.0f} MeV" if pid in pid_pt_map else "?"

        legend_label = f"[MISSED] PID {pid} ({species}) pT={pt_str}"

        node_indices = np.where(mask)[0]
        hover = [
            f"Missed particle {pid} ({species})<br>"
            f"pT: {pt_str}<br>"
            f"r: {r[i]:.1f} mm  phi: {phi[i]:.3f}  z: {z_cyl[i]:.1f} mm"
            for i in node_indices
        ]

        fig.add_trace(go.Scatter3d(
            x=x[mask], y=y[mask], z=z[mask],
            mode="markers",
            marker=dict(size=5, color=COLOR_MISSED, opacity=0.9,
                        symbol="diamond",
                        line=dict(width=0.5, color="darkorange")),
            name=legend_label,
            legendgroup=legend_label,
            text=hover,
            hoverinfo="text",
            showlegend=True,
        ))

        # GT trajectory edges — same legendgroup so they toggle together
        _add_gt_trajectory_edges(
            fig, pid, hit_particle_ids, hit_t, x, y, z,
            color=COLOR_MISSED, opacity=0.7, legendgroup=legend_label,
        )

    # ── Noise / unassigned hits (gray, no legend) ─────────────────────────────
    # Hits that are neither in a matched/fake track nor a missed particle
    assigned_mask = np.zeros(num_nodes, dtype=bool)
    for track_id in {**matched_tracks, **fake_tracks}:
        assigned_mask |= (hit_track_labels == track_id)
    for pid in missed_pids:
        assigned_mask |= (hit_particle_ids == pid)

    noise_mask = ~assigned_mask
    if noise_mask.any():
        fig.add_trace(go.Scatter3d(
            x=x[noise_mask], y=y[noise_mask], z=z[noise_mask],
            mode="markers",
            marker=dict(size=2, color="gray", opacity=0.2),
            name=f"Unassigned ({noise_mask.sum()} hits)",
            hoverinfo="skip",
            showlegend=True,
        ))

    # ── Summary statistics ────────────────────────────────────────────────────
    n_matched = len(matched_tracks)
    n_fake = len(fake_tracks)
    n_missed = len(missed_pids)
    n_total_reco = n_matched + n_fake
    eff = n_matched / max(n_matched + n_missed, 1)
    fake_rate = n_fake / max(n_total_reco, 1)

    event_id = graph.event_id
    if isinstance(event_id, list):
        event_id = event_id[0]

    title = (
        f"Evaluated Tracks: {dataset_name} graph #{index}  (event {event_id})<br>"
        f"Matched (green): {n_matched}  |  Fake (red): {n_fake}  |  "
        f"Missed (orange): {n_missed}  |  "
        f"Efficiency: {eff:.1%}  Fake rate: {fake_rate:.1%}"
    )

    layout = get_standard_scene_layout(title, margin_t=80, margin_l=0)
    # Override legend: spread it along the full left side of the screen
    layout["legend"] = dict(
        x=0,
        xanchor="left",
        y=1.0,
        yanchor="top",
        bgcolor="rgba(255,255,255,0.92)",
        bordercolor="rgba(0,0,0,0.15)",
        borderwidth=1,
        font=dict(size=11, family="monospace"),
        itemsizing="constant",
        tracegroupgap=4,      # tight vertical spacing so all items fit
    )
    # Push the 3D scene right to leave room for the left-side legend.
    # scene.domain shifts the 3D canvas; x=[0.22, 1] leaves ~22 % for the legend.
    layout["scene"]["domain"] = dict(x=[0.22, 1.0], y=[0.0, 1.0])
    layout["margin"] = dict(l=10, r=10, b=10, t=80)

    fig.update_layout(**layout)
    return fig


def _add_gt_trajectory_edges(fig, pid, hit_particle_ids, hit_t, x, y, z,
                               color, opacity, legendgroup):
    """Add time-ordered GT trajectory edges for one particle as a line trace.

    Shares legendgroup with the spacepoints trace so clicking the legend
    item hides/shows both together.
    """
    mask = hit_particle_ids == pid
    node_indices = np.where(mask)[0]
    if len(node_indices) < 2:
        return

    if hit_t is not None:
        order = np.argsort(hit_t[mask])
        node_indices = node_indices[order]

    edge_x, edge_y, edge_z = [], [], []
    for i in range(len(node_indices) - 1):
        a, b = node_indices[i], node_indices[i + 1]
        edge_x.extend([x[a], x[b], None])
        edge_y.extend([y[a], y[b], None])
        edge_z.extend([z[a], z[b], None])

    fig.add_trace(go.Scatter3d(
        x=edge_x, y=edge_y, z=edge_z,
        mode="lines",
        line=dict(color=color, width=2),
        opacity=opacity,
        legendgroup=legendgroup,
        name=legendgroup,          # must be set for legendgroup to work
        hoverinfo="skip",
        showlegend=False,          # only the spacepoints trace shows in legend
    ))


# ─── Main Entry ───────────────────────────────────────────────────────────────

def visualize_evaluated_tracks(dataset_name, index,
                                matching_fraction=0.7,
                                matching_style="two_way",
                                min_track_length=3,
                                min_particle_hits=3,
                                output_path=None):
    print("=" * 70)
    print(f"EVALUATED TRACK VISUALIZATION: {dataset_name} graph #{index}")
    print(f"  Matching: {matching_style}  threshold={matching_fraction}  min_particle_hits={min_particle_hits}")
    print("=" * 70)
    print()

    graph, _ = load_track_building_graph(dataset_name, index)

    print(f"  Nodes: {graph.num_nodes}")
    print()

    print("Running evaluation...")
    matching_df = evaluate_graph(graph, matching_fraction, matching_style,
                                 min_track_length, min_particle_hits)

    matched_tracks, fake_tracks, missed_pids = classify_tracks_from_matching(
        matching_df, np.array(graph.hit_track_labels.tolist())
    )
    print(f"  Matched tracks : {len(matched_tracks)}")
    print(f"  Fake tracks    : {len(fake_tracks)}")
    print(f"  Missed particles: {len(missed_pids)}")
    print()

    print("Creating visualization...")
    fig = create_evaluated_visualization(graph, matching_df, dataset_name, index)

    # Output path
    if output_path is None:
        out_dir = PIPELINE_ROOT / "data" / "visuals" / "evaluated_tracks" / dataset_name
        out_dir.mkdir(parents=True, exist_ok=True)
        output_path = out_dir / f"{dataset_name}{index:03d}_evaluated.html"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.write_html(str(output_path), include_plotlyjs="cdn")
    print(f"Saved: {output_path.absolute()}")
    print()
    print("=" * 70)

    return output_path


def main():
    parser = argparse.ArgumentParser(
        description="Visualize evaluated track building results (green/red/orange)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Color coding:
  Green  : Correctly matched tracks (species, particle ID, purity, completeness)
  Red    : Fake tracks (majority particle stats shown even below threshold)
  Orange : Missed particles (GT spacepoints for unrecovered particles)

Examples:
  python create_visual_evaluated_tracks.py testset 1
  python create_visual_evaluated_tracks.py valset 3 --matching-fraction 0.5
  python create_visual_evaluated_tracks.py testset 2 --matching-style ATLAS
        """,
    )
    parser.add_argument(
        "dataset",
        choices=["trainset", "valset", "testset"],
        help="Dataset to visualize",
    )
    parser.add_argument(
        "index",
        type=int,
        help="Graph index (1-based)",
    )
    parser.add_argument(
        "--matching-fraction",
        type=float,
        default=0.7,
        help="Minimum shared-hit fraction for matching (default: 0.7)",
    )
    parser.add_argument(
        "--matching-style",
        choices=["ATLAS", "one_way", "two_way"],
        default="two_way",
        help="Matching style (default: two_way)",
    )
    parser.add_argument(
        "--min-track-length",
        type=int,
        default=3,
        help="Minimum hits per reconstructed track (default: 3)",
    )
    parser.add_argument(
        "--min-particle-hits",
        type=int,
        default=3,
        help="Minimum true hits for a particle to be reconstructable (default: 3)",
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default=None,
        help="Output HTML path (default: data/visuals/evaluated_tracks/<dataset>/...)",
    )

    args = parser.parse_args()

    try:
        visualize_evaluated_tracks(
            args.dataset,
            args.index,
            matching_fraction=args.matching_fraction,
            matching_style=args.matching_style,
            min_track_length=args.min_track_length,
            min_particle_hits=args.min_particle_hits,
            output_path=args.output,
        )
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
