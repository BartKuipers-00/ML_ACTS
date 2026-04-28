"""
Segment-level evaluation utilities.

Evaluates the quality of proposed CC/Wrangler segments against ground truth
segments, before the matching stage. This diagnoses whether poor track-level
performance is due to bad segment extraction or bad matching.

For each true segment (particle_id, segment_id) passing fiducial cuts, we find
the proposed segment with the most shared hits and report:
  - completeness = shared_hits / n_true_hits
  - purity       = shared_hits / n_proposed_hits

No matching threshold is applied: if no proposed segment shares any hits with a
true segment, both metrics are 0.
"""

from collections import defaultdict
from pathlib import Path
from typing import List, Dict, Any

import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
import torch

from low_pt_custom_utils.mini_gnn_segment_embedding import get_segment_particle_id
from low_pt_custom_utils.plot_data_utils import save_plot_data_json

# PDG code → species name (same mapping as track_evaluation_utils)
PDG_TO_SPECIES = {
    13: 'Muon', -13: 'Muon',
    211: 'Pion', -211: 'Pion',
    11: 'Electron', -11: 'Electron',
    2212: 'Proton', -2212: 'Proton',
    321: 'Kaon', -321: 'Kaon',
}

SPECIES_COLORS = {
    'Muon':     '#1f77b4',
    'Pion':     '#ff7f0e',
    'Electron': '#2ca02c',
    'Proton':   '#d62728',
    'Kaon':     '#9467bd',
}

SPECIES_MARKERS = {
    'Muon':     'o',
    'Pion':     's',
    'Electron': '^',
    'Proton':   'D',
    'Kaon':     'v',
}


# ─── Per-Event Helpers ────────────────────────────────────────────────────────


def _build_pid_maps(graph) -> tuple:
    """
    Build particle_id → pT and particle_id → species maps from graph attributes.

    track_particle_pt / track_particle_type have one value per ground truth edge.
    We index the first hit of each edge to get its particle_id, then deduplicate.

    Returns: (pid_to_pt, pid_to_species)  — both keyed by raw ACTS particle_id.
    """
    pid_to_pt      = {}
    pid_to_species = {}

    if not (hasattr(graph, "track_particle_pt") and hasattr(graph, "track_edges")):
        return pid_to_pt, pid_to_species

    edge_pids = graph.hit_particle_id[graph.track_edges[0, :]].cpu().numpy()
    edge_pt   = graph.track_particle_pt.cpu().numpy()

    for pid, pt in zip(edge_pids, edge_pt):
        if pid > 0 and pid not in pid_to_pt:
            pid_to_pt[int(pid)] = float(pt)

    if hasattr(graph, "track_particle_type"):
        edge_ptype = graph.track_particle_type.cpu().numpy()
        for pid, ptype in zip(edge_pids, edge_ptype):
            if pid > 0 and pid not in pid_to_species:
                pid_to_species[int(pid)] = PDG_TO_SPECIES.get(int(ptype), "Other")

    return pid_to_pt, pid_to_species


# ─── Per-Event Evaluation ─────────────────────────────────────────────────────


def evaluate_segments_for_event(
    graph,
    proposed_segments,
    fiducial_config: Dict[str, Any],
) -> List[Dict[str, Any]]:
    """
    For each true segment passing fiducial cuts, find the best-matching proposed
    segment (most shared hits) and compute purity and completeness.

    Args:
        graph:             PyG Data object with hit_particle_id, hit_segment_id,
                           track_particle_pt, track_particle_type, track_edges.
        proposed_segments: List[SegmentInfo] from CC/Wrangler extraction.
        fiducial_config:   Dict with keys:
                             'pt'    : [pt_min, pt_max]   (GeV)
                             'nhits' : [nhits_min, nhits_max]

    Returns:
        List of dicts, one per true segment passing fiducial cuts:
          { 'pt', 'particle_type', 'nhits_true',
            'purity', 'completeness', 'has_match' }
    """
    if not hasattr(graph, "hit_segment_id"):
        return []

    hit_particle_id = graph.hit_particle_id.cpu().numpy().astype(np.int64)
    hit_segment_id  = graph.hit_segment_id.cpu().numpy().astype(np.int64)

    pt_range    = fiducial_config.get("pt",    [0.1, np.inf])
    nhits_range = fiducial_config.get("nhits", [3,   np.inf])
    pt_min,    pt_max    = float(pt_range[0]),  float(pt_range[1])
    nhits_min = int(nhits_range[0])
    nhits_max = int(nhits_range[1]) if nhits_range[1] != float("inf") else 10**9

    pid_to_pt, pid_to_species = _build_pid_maps(graph)

    # Build eta map (same pattern as evaluate_matching_efficiency_for_event)
    pid_to_eta: dict = {}
    if hasattr(graph, "track_particle_eta") and hasattr(graph, "track_edges"):
        edge_pids = graph.hit_particle_id[graph.track_edges[0, :]].cpu().numpy()
        edge_eta  = graph.track_particle_eta.cpu().numpy()
        for pid, eta in zip(edge_pids, edge_eta):
            if int(pid) > 0 and int(pid) not in pid_to_eta:
                pid_to_eta[int(pid)] = float(eta)

    # ── Group signal hits by (particle_id, segment_id) ──────────────────────
    # Raw ACTS particle IDs are ~10^15 — cannot pack into a single int64 key.
    signal_mask    = hit_particle_id > 0
    signal_pids    = hit_particle_id[signal_mask]
    signal_sids    = hit_segment_id[signal_mask]
    signal_indices = np.where(signal_mask)[0]

    seg_to_hits: dict = defaultdict(list)
    for idx, (pid, sid) in enumerate(zip(signal_pids.tolist(), signal_sids.tolist())):
        seg_to_hits[(pid, sid)].append(signal_indices[idx])

    # Pre-build proposed segment hit sets for fast intersection
    proposed_hit_sets = [set(seg.hits) for seg in proposed_segments]

    results = []
    for (pid, sid), hit_list in seg_to_hits.items():
        true_hit_indices = set(hit_list)
        nhits_true       = len(true_hit_indices)

        # Fiducial: hit count
        if nhits_true < nhits_min or nhits_true > nhits_max:
            continue

        # Fiducial: pT
        pt = pid_to_pt.get(pid, None)
        if pt is None:
            continue
        if pt < pt_min or pt > pt_max:
            continue

        # Find best proposed segment by shared hit count
        best_shared         = 0
        best_proposed_nhits = 0

        for proposed_set in proposed_hit_sets:
            shared = len(true_hit_indices & proposed_set)
            if shared > best_shared:
                best_shared         = shared
                best_proposed_nhits = len(proposed_set)

        completeness = best_shared / nhits_true         if nhits_true > 0         else 0.0
        purity       = best_shared / best_proposed_nhits if best_proposed_nhits > 0 else 0.0

        results.append({
            "pt":            pt,
            "eta":           pid_to_eta.get(pid, float("nan")),
            "particle_type": pid_to_species.get(pid, "Other"),
            "nhits_true":    nhits_true,
            "purity":        purity,
            "completeness":  completeness,
            "has_match":     best_shared > 0,
        })

    return results


# ─── Summary Text ─────────────────────────────────────────────────────────────


def make_segment_evaluation_summary(
    all_results: List[Dict[str, Any]],
    n_events: int,
) -> str:
    """Return a human-readable summary string for segment-level evaluation."""
    if not all_results:
        return "No true segments found (check hit_segment_id availability)."

    df = pd.DataFrame(all_results)
    n_true  = len(df)
    n_match = int(df["has_match"].sum())

    n_perfect_purity       = int((df["purity"]       == 1.0).sum())
    n_perfect_completeness = int((df["completeness"] == 1.0).sum())

    lines = [
        f"Segment Evaluation  ({n_events} events)",
        f"  True segments (fiducial):    {n_true}",
        f"  With ≥1 shared proposed hit: {n_match}  ({100*n_match/n_true:.1f}%)",
        f"",
        f"  Mean completeness:   {df['completeness'].mean():.4f}  "
        f"(std {df['completeness'].std():.4f})",
        f"  Perfect completeness:{n_perfect_completeness:6d} / {n_true}  "
        f"({100*n_perfect_completeness/n_true:.1f}%)",
        f"",
        f"  Mean purity:         {df['purity'].mean():.4f}  "
        f"(std {df['purity'].std():.4f})",
        f"  Perfect purity:      {n_perfect_purity:6d} / {n_true}  "
        f"({100*n_perfect_purity/n_true:.1f}%)",
    ]
    return "\n".join(lines)


# ─── Plotting helpers ─────────────────────────────────────────────────────────


def _clopper_pearson_errors(passed_arr, total_arr, level=0.68):
    """Clopper-Pearson confidence interval errors (same as track_evaluation_utils)."""
    alpha = (1 - level) / 2
    lo = np.where(passed_arr > 0,
                  scipy.stats.beta.ppf(alpha, passed_arr, total_arr - passed_arr + 1), 0.0)
    hi = np.where(passed_arr < total_arr,
                  scipy.stats.beta.ppf(1 - alpha, passed_arr + 1, total_arr - passed_arr), 1.0)
    eff = np.where(total_arr > 0, passed_arr / total_arr, np.nan)
    return np.maximum(0, eff - lo), np.maximum(0, hi - eff)


def _get_species_list(df):
    if "particle_type" not in df.columns:
        return []
    return sorted(df["particle_type"].dropna().unique())


# ─── Plotting ─────────────────────────────────────────────────────────────────


def plot_segment_evaluation(
    all_results: List[Dict[str, Any]],
    output_dir: Path,
    dataset_name: str,
    plot_config: Dict[str, Any],
) -> None:
    """
    Save three separate figures (efficiency, completeness, purity vs pT),
    one curve per particle species — same style as track_evaluation_utils plots.

    Args:
        all_results:  Flat list of per-segment dicts from evaluate_segments_for_event.
        output_dir:   Directory to save figures.
        dataset_name: Used in titles and filenames.
        plot_config:  Dict read from plots.segment_evaluation in the config YAML.
                      Expected sub-keys: title, font_sizes, variables.pt
    """
    if not all_results:
        print("  No segment evaluation data to plot.")
        return

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(all_results)

    fs       = plot_config.get("font_sizes", {})
    fs_ax    = fs.get("axis_label", 18)
    fs_tk    = fs.get("tick_label", 14)
    fs_ti    = fs.get("title",      18)
    fs_lg    = fs.get("legend",     14)

    pt_var   = plot_config.get("variables", {}).get("pt", {})
    x_lim    = pt_var.get("x_lim",   [0.1, 0.5])
    y_lim    = pt_var.get("y_lim",   [0.0, 1.05])
    n_bins   = pt_var.get("n_bins",  30)
    x_label  = pt_var.get("x_label", "$p_T$ [GeV]")

    x_bins = np.linspace(x_lim[0], x_lim[1], n_bins + 1)
    xvals  = (x_bins[1:] + x_bins[:-1]) / 2
    xerrs  = (x_bins[1:] - x_bins[:-1]) / 2

    species_list = _get_species_list(df) or ["All"]

    # ── Overall scores for titles ─────────────────────────────────────────────
    overall_efficiency          = float(df["has_match"].mean())
    overall_completeness        = float(df["completeness"].mean())
    overall_purity              = float(df["purity"].mean())
    overall_perfect_purity      = float((df["purity"]       == 1.0).mean())
    overall_perfect_completeness = float((df["completeness"] == 1.0).mean())

    # ── Helper ───────────────────────────────────────────────────────────────

    def _make_fig(ylabel, title):
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.set_xlabel(x_label,  fontsize=fs_ax)
        ax.set_ylabel(ylabel,   fontsize=fs_ax)
        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim)
        ax.tick_params(axis="both", labelsize=fs_tk)
        ax.grid(True, alpha=0.3)
        ax.set_title(title, fontsize=fs_ti)
        return fig, ax

    def _save(fig, ax, suffix, plot_data):
        ax.legend(fontsize=fs_lg)
        plt.tight_layout()
        path = output_dir / f"segment_{suffix}_{dataset_name}.png"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
        save_plot_data_json(path, plot_data)
        plt.close(fig)
        print(f"  Saved: {path}")

    # ── 1. Segment Finding Efficiency ────────────────────────────────────────
    fig_eff, ax_eff = _make_fig(
        "Segment Finding Efficiency",
        f"Segment Finding Efficiency (Overall: {overall_efficiency:.3f})",
    )
    pd_eff = {"plot_type": "segment_efficiency", "xlabel": x_label,
              "ylabel": "Segment Finding Efficiency", "xlim": x_lim, "ylim": y_lim,
              "title": f"Segment Finding Efficiency (Overall: {overall_efficiency:.3f})",
              "series": []}
    for species in species_list:
        mask = np.ones(len(df), dtype=bool) if species == "All" \
               else (df["particle_type"] == species).values

        sub       = df[mask]
        pt_all    = sub["pt"].values
        pt_match  = sub.loc[sub["has_match"], "pt"].values

        total_vals, _ = np.histogram(pt_all,   bins=x_bins)
        match_vals, _ = np.histogram(pt_match, bins=x_bins)

        with np.errstate(divide="ignore", invalid="ignore"):
            eff = np.where(total_vals > 0, match_vals / total_vals, np.nan)
        lo_err, hi_err = _clopper_pearson_errors(match_vals, total_vals)

        n_tot = int(total_vals.sum())
        n_mat = int(match_vals.sum())
        label  = f"{species} ({n_mat}/{n_tot})"
        color  = SPECIES_COLORS.get(species, "black")
        marker = SPECIES_MARKERS.get(species, "o")

        ax_eff.errorbar(xvals, eff, xerr=xerrs, yerr=[lo_err, hi_err],
                        fmt=marker, color=color, label=label,
                        capsize=3, capthick=1.5, markersize=5)
        pd_eff["series"].append({"label": label, "color": color, "marker": marker,
                                  "x": xvals, "xerr": xerrs,
                                  "y": eff, "yerr_lo": lo_err, "yerr_hi": hi_err})

    _save(fig_eff, ax_eff, "efficiency", pd_eff)

    # ── 2. Mean Completeness ─────────────────────────────────────────────────
    fig_com, ax_com = _make_fig(
        "Completeness (fraction of true hits found)",
        f"Segment Completeness (Overall: {overall_completeness:.3f})",
    )
    pd_com = {"plot_type": "segment_completeness", "xlabel": x_label,
              "ylabel": "Completeness (fraction of true hits found)",
              "xlim": x_lim, "ylim": y_lim,
              "title": f"Segment Completeness (Overall: {overall_completeness:.3f})",
              "series": []}
    for species in species_list:
        mask = np.ones(len(df), dtype=bool) if species == "All" \
               else (df["particle_type"] == species).values

        sub = df[mask]
        sx  = sub["pt"].values
        sv  = sub["completeness"].values

        if len(sx) == 0:
            continue

        bin_idx   = np.digitize(sx, x_bins) - 1
        mean_vals = np.full(n_bins, np.nan)
        err_vals  = np.full(n_bins, np.nan)

        for i in range(n_bins):
            bm = bin_idx == i
            if bm.sum() > 0:
                bv = sv[bm]
                mean_vals[i] = np.mean(bv)
                err_vals[i]  = np.std(bv) / np.sqrt(len(bv)) if len(bv) > 1 else 0.0

        color  = SPECIES_COLORS.get(species, "black")
        marker = SPECIES_MARKERS.get(species, "o")
        ax_com.errorbar(xvals, mean_vals, xerr=xerrs, yerr=err_vals,
                        fmt=marker, color=color, label=species,
                        capsize=3, capthick=1.5, markersize=5)
        pd_com["series"].append({"label": species, "color": color, "marker": marker,
                                  "x": xvals, "xerr": xerrs,
                                  "y": mean_vals, "yerr_lo": err_vals, "yerr_hi": err_vals})

    _save(fig_com, ax_com, "completeness", pd_com)

    # ── 3. Mean Purity ───────────────────────────────────────────────────────
    fig_pur, ax_pur = _make_fig(
        "Purity (fraction of correct hits in proposed segment)",
        f"Segment Purity (Overall: {overall_purity:.3f})",
    )
    pd_pur = {"plot_type": "segment_purity", "xlabel": x_label,
              "ylabel": "Purity (fraction of correct hits in proposed segment)",
              "xlim": x_lim, "ylim": y_lim,
              "title": f"Segment Purity (Overall: {overall_purity:.3f})",
              "series": []}
    for species in species_list:
        mask = np.ones(len(df), dtype=bool) if species == "All" \
               else (df["particle_type"] == species).values

        sub = df[mask]
        sx  = sub["pt"].values
        sv  = sub["purity"].values

        if len(sx) == 0:
            continue

        bin_idx   = np.digitize(sx, x_bins) - 1
        mean_vals = np.full(n_bins, np.nan)
        err_vals  = np.full(n_bins, np.nan)

        for i in range(n_bins):
            bm = bin_idx == i
            if bm.sum() > 0:
                bv = sv[bm]
                mean_vals[i] = np.mean(bv)
                err_vals[i]  = np.std(bv) / np.sqrt(len(bv)) if len(bv) > 1 else 0.0

        color  = SPECIES_COLORS.get(species, "black")
        marker = SPECIES_MARKERS.get(species, "o")
        ax_pur.errorbar(xvals, mean_vals, xerr=xerrs, yerr=err_vals,
                        fmt=marker, color=color, label=species,
                        capsize=3, capthick=1.5, markersize=5)
        pd_pur["series"].append({"label": species, "color": color, "marker": marker,
                                  "x": xvals, "xerr": xerrs,
                                  "y": mean_vals, "yerr_lo": err_vals, "yerr_hi": err_vals})

    _save(fig_pur, ax_pur, "purity", pd_pur)

    # ── 4. Fraction with perfect purity (purity == 1.0) ──────────────────────
    fig_pp, ax_pp = _make_fig(
        "Fraction with Perfect Purity",
        f"Segment Perfect Purity (Overall: {overall_perfect_purity:.3f})",
    )
    pd_pp = {"plot_type": "segment_perfect_purity", "xlabel": x_label,
             "ylabel": "Fraction with Perfect Purity",
             "xlim": x_lim, "ylim": y_lim,
             "title": f"Segment Perfect Purity (Overall: {overall_perfect_purity:.3f})",
             "series": []}
    for species in species_list:
        mask = np.ones(len(df), dtype=bool) if species == "All" \
               else (df["particle_type"] == species).values

        sub          = df[mask]
        pt_all       = sub["pt"].values
        pt_perfect   = sub.loc[sub["purity"] == 1.0, "pt"].values

        total_vals, _ = np.histogram(pt_all,     bins=x_bins)
        perf_vals,  _ = np.histogram(pt_perfect, bins=x_bins)

        with np.errstate(divide="ignore", invalid="ignore"):
            frac = np.where(total_vals > 0, perf_vals / total_vals, np.nan)
        lo_err, hi_err = _clopper_pearson_errors(perf_vals, total_vals)

        n_tot  = int(total_vals.sum())
        n_perf = int(perf_vals.sum())
        label  = f"{species} ({n_perf}/{n_tot})"
        color  = SPECIES_COLORS.get(species, "black")
        marker = SPECIES_MARKERS.get(species, "o")

        ax_pp.errorbar(xvals, frac, xerr=xerrs, yerr=[lo_err, hi_err],
                       fmt=marker, color=color, label=label,
                       capsize=3, capthick=1.5, markersize=5)
        pd_pp["series"].append({"label": label, "color": color, "marker": marker,
                                 "x": xvals, "xerr": xerrs,
                                 "y": frac, "yerr_lo": lo_err, "yerr_hi": hi_err})

    _save(fig_pp, ax_pp, "perfect_purity", pd_pp)

    # ── 5. Fraction with perfect completeness (completeness == 1.0) ──────────
    fig_pc, ax_pc = _make_fig(
        "Fraction with Perfect Completeness",
        f"Segment Perfect Completeness (Overall: {overall_perfect_completeness:.3f})",
    )
    pd_pc = {"plot_type": "segment_perfect_completeness", "xlabel": x_label,
             "ylabel": "Fraction with Perfect Completeness",
             "xlim": x_lim, "ylim": y_lim,
             "title": f"Segment Perfect Completeness (Overall: {overall_perfect_completeness:.3f})",
             "series": []}
    for species in species_list:
        mask = np.ones(len(df), dtype=bool) if species == "All" \
               else (df["particle_type"] == species).values

        sub        = df[mask]
        pt_all     = sub["pt"].values
        pt_perfect = sub.loc[sub["completeness"] == 1.0, "pt"].values

        total_vals, _ = np.histogram(pt_all,     bins=x_bins)
        perf_vals,  _ = np.histogram(pt_perfect, bins=x_bins)

        with np.errstate(divide="ignore", invalid="ignore"):
            frac = np.where(total_vals > 0, perf_vals / total_vals, np.nan)
        lo_err, hi_err = _clopper_pearson_errors(perf_vals, total_vals)

        n_tot  = int(total_vals.sum())
        n_perf = int(perf_vals.sum())
        label  = f"{species} ({n_perf}/{n_tot})"
        color  = SPECIES_COLORS.get(species, "black")
        marker = SPECIES_MARKERS.get(species, "o")

        ax_pc.errorbar(xvals, frac, xerr=xerrs, yerr=[lo_err, hi_err],
                       fmt=marker, color=color, label=label,
                       capsize=3, capthick=1.5, markersize=5)
        pd_pc["series"].append({"label": label, "color": color, "marker": marker,
                                 "x": xvals, "xerr": xerrs,
                                 "y": frac, "yerr_lo": lo_err, "yerr_hi": hi_err})

    _save(fig_pc, ax_pc, "perfect_completeness", pd_pc)

    # ── Eta plots (if eta data available) ─────────────────────────────────────
    if "eta" not in df.columns or df["eta"].isna().all():
        return

    eta_var   = plot_config.get("variables", {}).get("eta", {})
    y_lim_e   = eta_var.get("y_lim",   [0.0, 1.05])
    x_label_e = eta_var.get("x_label", r"$\eta$")
    n_bins_e  = eta_var.get("n_bins",  30)

    eta_vals_all = df["eta"].dropna().values
    eta_lo = float(np.floor(eta_vals_all.min() * 10) / 10)
    eta_hi = float(np.ceil( eta_vals_all.max() * 10) / 10)
    x_lim_e = [eta_lo, eta_hi]

    x_bins_e = np.linspace(x_lim_e[0], x_lim_e[1], n_bins_e + 1)
    xvals_e  = (x_bins_e[1:] + x_bins_e[:-1]) / 2
    xerrs_e  = (x_bins_e[1:] - x_bins_e[:-1]) / 2

    def _make_fig_eta(ylabel, title):
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.set_xlabel(x_label_e, fontsize=fs_ax)
        ax.set_ylabel(ylabel,    fontsize=fs_ax)
        ax.set_xlim(x_lim_e)
        ax.set_ylim(y_lim_e)
        ax.tick_params(axis="both", labelsize=fs_tk)
        ax.grid(True, alpha=0.3)
        ax.set_title(title, fontsize=fs_ti)
        return fig, ax

    # 6. Segment Finding Efficiency vs eta
    fig, ax = _make_fig_eta(
        "Segment Finding Efficiency",
        f"Segment Finding Efficiency vs η (Overall: {overall_efficiency:.3f})",
    )
    pd_e = {"plot_type": "segment_efficiency_eta", "xlabel": x_label_e,
            "ylabel": "Segment Finding Efficiency", "xlim": x_lim_e, "ylim": y_lim_e,
            "title": f"Segment Finding Efficiency vs η (Overall: {overall_efficiency:.3f})",
            "series": []}
    for species in species_list:
        mask = np.ones(len(df), dtype=bool) if species == "All" \
               else (df["particle_type"] == species).values
        sub      = df[mask].dropna(subset=["eta"])
        eta_all  = sub["eta"].values
        eta_match = sub.loc[sub["has_match"], "eta"].values
        total_vals, _ = np.histogram(eta_all,   bins=x_bins_e)
        match_vals, _ = np.histogram(eta_match, bins=x_bins_e)
        with np.errstate(divide="ignore", invalid="ignore"):
            eff = np.where(total_vals > 0, match_vals / total_vals, np.nan)
        lo_err, hi_err = _clopper_pearson_errors(match_vals, total_vals)
        n_tot = int(total_vals.sum()); n_mat = int(match_vals.sum())
        label = f"{species} ({n_mat}/{n_tot})"
        color = SPECIES_COLORS.get(species, "black"); marker = SPECIES_MARKERS.get(species, "o")
        ax.errorbar(xvals_e, eff, xerr=xerrs_e, yerr=[lo_err, hi_err],
                    fmt=marker, color=color, label=label, capsize=3, capthick=1.5, markersize=5)
        pd_e["series"].append({"label": label, "color": color, "marker": marker,
                                "x": xvals_e, "xerr": xerrs_e,
                                "y": eff, "yerr_lo": lo_err, "yerr_hi": hi_err})
    _save(fig, ax, "efficiency_eta", pd_e)

    # 7. Completeness vs eta
    fig, ax = _make_fig_eta(
        "Completeness (fraction of true hits found)",
        f"Segment Completeness vs η (Overall: {overall_completeness:.3f})",
    )
    pd_e = {"plot_type": "segment_completeness_eta", "xlabel": x_label_e,
            "ylabel": "Completeness (fraction of true hits found)",
            "xlim": x_lim_e, "ylim": y_lim_e,
            "title": f"Segment Completeness vs η (Overall: {overall_completeness:.3f})",
            "series": []}
    for species in species_list:
        mask = np.ones(len(df), dtype=bool) if species == "All" \
               else (df["particle_type"] == species).values
        sub = df[mask].dropna(subset=["eta"])
        sx = sub["eta"].values; sv = sub["completeness"].values
        if len(sx) == 0:
            continue
        bin_idx   = np.digitize(sx, x_bins_e) - 1
        mean_vals = np.full(n_bins_e, np.nan); err_vals = np.full(n_bins_e, np.nan)
        for i in range(n_bins_e):
            bm = bin_idx == i
            if bm.sum() > 0:
                bv = sv[bm]; mean_vals[i] = np.mean(bv)
                err_vals[i] = np.std(bv) / np.sqrt(len(bv)) if len(bv) > 1 else 0.0
        color = SPECIES_COLORS.get(species, "black"); marker = SPECIES_MARKERS.get(species, "o")
        ax.errorbar(xvals_e, mean_vals, xerr=xerrs_e, yerr=err_vals,
                    fmt=marker, color=color, label=species, capsize=3, capthick=1.5, markersize=5)
        pd_e["series"].append({"label": species, "color": color, "marker": marker,
                                "x": xvals_e, "xerr": xerrs_e,
                                "y": mean_vals, "yerr_lo": err_vals, "yerr_hi": err_vals})
    _save(fig, ax, "completeness_eta", pd_e)

    # 8. Purity vs eta
    fig, ax = _make_fig_eta(
        "Purity (fraction of correct hits in proposed segment)",
        f"Segment Purity vs η (Overall: {overall_purity:.3f})",
    )
    pd_e = {"plot_type": "segment_purity_eta", "xlabel": x_label_e,
            "ylabel": "Purity (fraction of correct hits in proposed segment)",
            "xlim": x_lim_e, "ylim": y_lim_e,
            "title": f"Segment Purity vs η (Overall: {overall_purity:.3f})",
            "series": []}
    for species in species_list:
        mask = np.ones(len(df), dtype=bool) if species == "All" \
               else (df["particle_type"] == species).values
        sub = df[mask].dropna(subset=["eta"])
        sx = sub["eta"].values; sv = sub["purity"].values
        if len(sx) == 0:
            continue
        bin_idx   = np.digitize(sx, x_bins_e) - 1
        mean_vals = np.full(n_bins_e, np.nan); err_vals = np.full(n_bins_e, np.nan)
        for i in range(n_bins_e):
            bm = bin_idx == i
            if bm.sum() > 0:
                bv = sv[bm]; mean_vals[i] = np.mean(bv)
                err_vals[i] = np.std(bv) / np.sqrt(len(bv)) if len(bv) > 1 else 0.0
        color = SPECIES_COLORS.get(species, "black"); marker = SPECIES_MARKERS.get(species, "o")
        ax.errorbar(xvals_e, mean_vals, xerr=xerrs_e, yerr=err_vals,
                    fmt=marker, color=color, label=species, capsize=3, capthick=1.5, markersize=5)
        pd_e["series"].append({"label": species, "color": color, "marker": marker,
                                "x": xvals_e, "xerr": xerrs_e,
                                "y": mean_vals, "yerr_lo": err_vals, "yerr_hi": err_vals})
    _save(fig, ax, "purity_eta", pd_e)

    # 9. Perfect Purity vs eta
    fig, ax = _make_fig_eta(
        "Fraction with Perfect Purity",
        f"Segment Perfect Purity vs η (Overall: {overall_perfect_purity:.3f})",
    )
    pd_e = {"plot_type": "segment_perfect_purity_eta", "xlabel": x_label_e,
            "ylabel": "Fraction with Perfect Purity",
            "xlim": x_lim_e, "ylim": y_lim_e,
            "title": f"Segment Perfect Purity vs η (Overall: {overall_perfect_purity:.3f})",
            "series": []}
    for species in species_list:
        mask = np.ones(len(df), dtype=bool) if species == "All" \
               else (df["particle_type"] == species).values
        sub = df[mask].dropna(subset=["eta"])
        eta_all    = sub["eta"].values
        eta_perfect = sub.loc[sub["purity"] == 1.0, "eta"].values
        total_vals, _ = np.histogram(eta_all,     bins=x_bins_e)
        perf_vals,  _ = np.histogram(eta_perfect, bins=x_bins_e)
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = np.where(total_vals > 0, perf_vals / total_vals, np.nan)
        lo_err, hi_err = _clopper_pearson_errors(perf_vals, total_vals)
        n_tot = int(total_vals.sum()); n_perf = int(perf_vals.sum())
        label = f"{species} ({n_perf}/{n_tot})"
        color = SPECIES_COLORS.get(species, "black"); marker = SPECIES_MARKERS.get(species, "o")
        ax.errorbar(xvals_e, frac, xerr=xerrs_e, yerr=[lo_err, hi_err],
                    fmt=marker, color=color, label=label, capsize=3, capthick=1.5, markersize=5)
        pd_e["series"].append({"label": label, "color": color, "marker": marker,
                                "x": xvals_e, "xerr": xerrs_e,
                                "y": frac, "yerr_lo": lo_err, "yerr_hi": hi_err})
    _save(fig, ax, "perfect_purity_eta", pd_e)

    # 10. Perfect Completeness vs eta
    fig, ax = _make_fig_eta(
        "Fraction with Perfect Completeness",
        f"Segment Perfect Completeness vs η (Overall: {overall_perfect_completeness:.3f})",
    )
    pd_e = {"plot_type": "segment_perfect_completeness_eta", "xlabel": x_label_e,
            "ylabel": "Fraction with Perfect Completeness",
            "xlim": x_lim_e, "ylim": y_lim_e,
            "title": f"Segment Perfect Completeness vs η (Overall: {overall_perfect_completeness:.3f})",
            "series": []}
    for species in species_list:
        mask = np.ones(len(df), dtype=bool) if species == "All" \
               else (df["particle_type"] == species).values
        sub = df[mask].dropna(subset=["eta"])
        eta_all    = sub["eta"].values
        eta_perfect = sub.loc[sub["completeness"] == 1.0, "eta"].values
        total_vals, _ = np.histogram(eta_all,     bins=x_bins_e)
        perf_vals,  _ = np.histogram(eta_perfect, bins=x_bins_e)
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = np.where(total_vals > 0, perf_vals / total_vals, np.nan)
        lo_err, hi_err = _clopper_pearson_errors(perf_vals, total_vals)
        n_tot = int(total_vals.sum()); n_perf = int(perf_vals.sum())
        label = f"{species} ({n_perf}/{n_tot})"
        color = SPECIES_COLORS.get(species, "black"); marker = SPECIES_MARKERS.get(species, "o")
        ax.errorbar(xvals_e, frac, xerr=xerrs_e, yerr=[lo_err, hi_err],
                    fmt=marker, color=color, label=label, capsize=3, capthick=1.5, markersize=5)
        pd_e["series"].append({"label": label, "color": color, "marker": marker,
                                "x": xvals_e, "xerr": xerrs_e,
                                "y": frac, "yerr_lo": lo_err, "yerr_hi": hi_err})
    _save(fig, ax, "perfect_completeness_eta", pd_e)


# ─── Segment Matching Efficiency ──────────────────────────────────────────────


def evaluate_matching_efficiency_for_event(
    graph,
    segments: list,
    matched_tracks: list,
    unmatched: list,
    fiducial_config: Dict[str, Any],
) -> List[Dict[str, Any]]:
    """
    For each looping particle (≥2 proposed segments), check whether all its
    segments were correctly grouped into a single matched track.

    A segment's dominant particle is determined by majority-vote on
    hit_particle_id (noise hits with pid=0 are excluded).

    'Correctly matched' for a particle with N segments: all N segments
    appear together in a single entry of matched_tracks.

    Args:
        graph:          PyG Data object with hit_particle_id, track_particle_pt,
                        track_particle_type, track_edges.
        segments:       Full list of SegmentInfo (all segments in this event).
        matched_tracks: List of lists of SegmentInfo (each list = one matched track).
        unmatched:      List of unmatched SegmentInfo.
        fiducial_config: Dict with keys 'pt' and 'nhits'.

    Returns:
        List of dicts, one per looping particle:
          { 'pt', 'eta', 'particle_type', 'n_segments', 'correctly_matched' }
    """
    if not hasattr(graph, "hit_particle_id"):
        return []

    pt_range         = fiducial_config.get("pt", [0.1, float("inf")])
    min_segment_hits = int(fiducial_config.get("min_segment_hits", 3))
    pt_min, pt_max   = float(pt_range[0]), float(pt_range[1])

    pid_to_pt, pid_to_species = _build_pid_maps(graph)

    # Build eta map
    pid_to_eta: Dict[int, float] = {}
    if hasattr(graph, "track_particle_eta") and hasattr(graph, "track_edges"):
        edge_pids = graph.hit_particle_id[graph.track_edges[0, :]].cpu().numpy()
        edge_eta  = graph.track_particle_eta.cpu().numpy()
        for pid, eta in zip(edge_pids, edge_eta):
            if int(pid) > 0 and int(pid) not in pid_to_eta:
                pid_to_eta[int(pid)] = float(eta)

    # Build ground-truth looper set: particles that truly have ≥2 distinct segment IDs.
    # This gates the denominator on physics, not on the GNN output — a non-looping
    # particle whose track is broken into 2 CC clusters must not appear here.
    gt_looper_pids: set = set()
    if hasattr(graph, "hit_segment_id"):
        hit_pids = graph.hit_particle_id.cpu().numpy().astype(np.int64)
        hit_sids = graph.hit_segment_id.cpu().numpy().astype(np.int64)
        pid_to_gt_segs: Dict[int, set] = defaultdict(set)
        for pid, sid in zip(hit_pids, hit_sids):
            if pid > 0:
                pid_to_gt_segs[int(pid)].add(int(sid))
        gt_looper_pids = {pid for pid, sids in pid_to_gt_segs.items() if len(sids) >= 2}

    # Label each segment with its dominant particle (reuse existing helper)
    seg_dominant: List[int] = [
        get_segment_particle_id(seg, graph) for seg in segments
    ]

    # Build pid → list of segment indices (all segments, including unmatched)
    pid_to_seg_indices: Dict[int, List[int]] = defaultdict(list)
    for seg_idx, pid in enumerate(seg_dominant):
        if pid > 0:
            pid_to_seg_indices[pid].append(seg_idx)

    # Build a lookup: segment object → matched_track index (or -1 if unmatched)
    # Use id(seg) as key to avoid equality issues with dataclasses
    seg_id_to_track: Dict[int, int] = {}
    for track_idx, track_segs in enumerate(matched_tracks):
        for seg in track_segs:
            seg_id_to_track[id(seg)] = track_idx
    # unmatched → no shared track
    for seg in unmatched:
        seg_id_to_track[id(seg)] = -1 - id(seg)  # unique negative sentinel

    results = []
    for pid, seg_indices in pid_to_seg_indices.items():
        if len(seg_indices) < 2:
            continue  # Only evaluate particles with ≥2 proposed segments

        # Gate on ground-truth: skip particles that don't truly loop.
        # A non-looping track broken by the GNN into 2 CC clusters would
        # otherwise pollute the denominator (e.g. high-pT muons above 0.35 GeV).
        if gt_looper_pids and pid not in gt_looper_pids:
            continue

        pt = pid_to_pt.get(pid, None)
        if pt is None or pt < pt_min or pt > pt_max:
            continue

        # Fiducial: each segment must be individually reconstructable
        if any(len(segments[i].hits) < min_segment_hits for i in seg_indices):
            continue

        # Check if all segments landed in the same matched_track entry
        track_ids = [seg_id_to_track[id(segments[i])] for i in seg_indices]
        # Correctly matched: all segments share the same positive track index
        correctly_matched = (len(set(track_ids)) == 1) and (track_ids[0] >= 0)

        results.append({
            "pt":                pt,
            "eta":               pid_to_eta.get(pid, float("nan")),
            "particle_type":     pid_to_species.get(pid, "Other"),
            "n_segments":        len(seg_indices),
            "correctly_matched": correctly_matched,
        })

    return results


def make_matching_efficiency_summary(
    all_results: List[Dict[str, Any]],
    n_events: int = None,
) -> str:
    """Return a human-readable summary of segment matching efficiency."""
    if not all_results:
        return "No looping particles (≥2 segments) found in this dataset."

    df = pd.DataFrame(all_results)
    n_total   = len(df)
    n_correct = int(df["correctly_matched"].sum())
    eff       = n_correct / n_total if n_total > 0 else 0.0

    header = f"Segment Matching Efficiency  ({n_events} events)" if n_events is not None \
             else "Segment Matching Efficiency"
    lines = [
        header,
        f"  Looping particles (≥2 segs): {n_total}",
        f"  Correctly matched:           {n_correct}  ({100*eff:.1f}%)",
    ]

    species_list = sorted(df["particle_type"].dropna().unique())
    if len(species_list) > 1:
        lines.append("")
        for sp in species_list:
            sub = df[df["particle_type"] == sp]
            n_sp  = len(sub)
            n_ok  = int(sub["correctly_matched"].sum())
            sp_eff = n_ok / n_sp if n_sp > 0 else 0.0
            lines.append(f"  {sp:<12s}  {n_ok:5d}/{n_sp:<5d}  ({100*sp_eff:.1f}%)")

    return "\n".join(lines)


def plot_matching_efficiency(
    all_results: List[Dict[str, Any]],
    output_dir: Path,
    dataset_name: str,
    plot_config: Dict[str, Any],
) -> None:
    """
    Plot segment matching efficiency vs pT (and optionally η) per particle species.

    Saves:
        matching_efficiency_pt_<dataset>.png
        matching_efficiency_eta_<dataset>.png  (if eta data available)
    """
    if not all_results:
        print("  No segment matching efficiency data to plot.")
        return

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(all_results)

    fs      = plot_config.get("font_sizes", {})
    fs_ax   = fs.get("axis_label", 18)
    fs_tk   = fs.get("tick_label", 14)
    fs_ti   = fs.get("title",      18)
    fs_lg   = fs.get("legend",     14)

    overall_eff = float(df["correctly_matched"].mean())
    species_list = _get_species_list(df) or ["All"]

    min_bin_fraction = plot_config.get("min_bin_fraction", 0.05)

    # ── Helper ───────────────────────────────────────────────────────────────
    def _plot_eff_vs_x(x_vals_all, x_matched_all, x_lim, x_label, n_bins,
                       ax, species, color, marker):
        x_bins = np.linspace(x_lim[0], x_lim[1], n_bins + 1)
        xvals  = (x_bins[1:] + x_bins[:-1]) / 2
        xerrs  = (x_bins[1:] - x_bins[:-1]) / 2

        total_raw, _ = np.histogram(x_vals_all,   bins=x_bins)
        match_raw, _ = np.histogram(x_matched_all, bins=x_bins)

        # Suppress bins with too few entries relative to the peak bin
        peak = total_raw.max() if total_raw.max() > 0 else 1
        sparse = total_raw < min_bin_fraction * peak
        total_vals = total_raw.astype(float)
        total_vals[sparse] = np.nan
        match_vals = match_raw.astype(float)
        match_vals[sparse] = np.nan

        with np.errstate(divide="ignore", invalid="ignore"):
            eff = np.where(~sparse & (total_vals > 0), match_vals / total_vals, np.nan)
        lo_err, hi_err = _clopper_pearson_errors(
            np.nan_to_num(match_vals).astype(int),
            np.nan_to_num(total_vals).astype(int),
        )
        lo_err[sparse] = np.nan
        hi_err[sparse] = np.nan

        n_tot = int(np.nansum(total_vals))
        n_mat = int(np.nansum(match_vals))
        label = f"{species} ({n_mat}/{n_tot})"
        ax.errorbar(xvals, eff, xerr=xerrs, yerr=[lo_err, hi_err],
                    fmt=marker, color=color, label=label,
                    capsize=3, capthick=1.5, markersize=5)
        return {"label": label, "color": color, "marker": marker,
                "x": xvals, "xerr": xerrs,
                "y": eff, "yerr_lo": lo_err, "yerr_hi": hi_err,
                "n": total_raw}

    # ── 1. Efficiency vs pT ───────────────────────────────────────────────────
    pt_var  = plot_config.get("variables", {}).get("pt", {})
    x_lim   = pt_var.get("x_lim",   [0.1, 0.5])
    n_bins  = pt_var.get("n_bins",  30)
    x_label = pt_var.get("x_label", "$p_T$ [GeV]")
    y_lim   = pt_var.get("y_lim",   [0.0, 1.05])

    pt_title = f"Segment Matching Efficiency (Overall: {overall_eff:.3f})"
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xlabel(x_label, fontsize=fs_ax)
    ax.set_ylabel("Segment Matching Efficiency", fontsize=fs_ax)
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.tick_params(axis="both", labelsize=fs_tk)
    ax.grid(True, alpha=0.3)
    ax.set_title(pt_title, fontsize=fs_ti)

    pd_pt = {"plot_type": "matching_efficiency", "xlabel": x_label,
             "ylabel": "Segment Matching Efficiency",
             "xlim": x_lim, "ylim": y_lim, "title": pt_title, "series": []}
    for species in species_list:
        mask = np.ones(len(df), dtype=bool) if species == "All" \
               else (df["particle_type"] == species).values
        sub      = df[mask]
        pt_all   = sub["pt"].values
        pt_match = sub.loc[sub["correctly_matched"], "pt"].values
        color    = SPECIES_COLORS.get(species, "black")
        marker   = SPECIES_MARKERS.get(species, "o")
        series = _plot_eff_vs_x(pt_all, pt_match, x_lim, x_label, n_bins, ax, species, color, marker)
        pd_pt["series"].append(series)

    ax.legend(fontsize=fs_lg)
    plt.tight_layout()
    path = output_dir / f"matching_efficiency_pt_{dataset_name}.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    save_plot_data_json(path, pd_pt)
    plt.close(fig)
    print(f"  Saved: {path}")

    # ── 2. Efficiency vs η ────────────────────────────────────────────────────
    if "eta" not in df.columns or df["eta"].isna().all():
        return

    eta_var   = plot_config.get("variables", {}).get("eta", {})
    y_lim_e   = eta_var.get("y_lim",   [0.0, 1.05])
    x_label_e = eta_var.get("x_label", r"$\eta$")
    n_bins_e  = eta_var.get("n_bins",  30)

    eta_vals_all = df["eta"].dropna().values
    eta_lo = float(np.floor(eta_vals_all.min() * 10) / 10)
    eta_hi = float(np.ceil( eta_vals_all.max() * 10) / 10)
    x_lim_e = [eta_lo, eta_hi]

    eta_title = f"Segment Matching Efficiency vs η (Overall: {overall_eff:.3f})"
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xlabel(x_label_e, fontsize=fs_ax)
    ax.set_ylabel("Segment Matching Efficiency", fontsize=fs_ax)
    ax.set_xlim(x_lim_e)
    ax.set_ylim(y_lim_e)
    ax.tick_params(axis="both", labelsize=fs_tk)
    ax.grid(True, alpha=0.3)
    ax.set_title(eta_title, fontsize=fs_ti)

    pd_eta = {"plot_type": "matching_efficiency", "xlabel": x_label_e,
              "ylabel": "Segment Matching Efficiency",
              "xlim": x_lim_e, "ylim": y_lim_e, "title": eta_title, "series": []}
    for species in species_list:
        mask = np.ones(len(df), dtype=bool) if species == "All" \
               else (df["particle_type"] == species).values
        sub       = df[mask].dropna(subset=["eta"])
        eta_all   = sub["eta"].values
        eta_match = sub.loc[sub["correctly_matched"], "eta"].values
        color     = SPECIES_COLORS.get(species, "black")
        marker    = SPECIES_MARKERS.get(species, "o")
        series = _plot_eff_vs_x(eta_all, eta_match, x_lim_e, x_label_e, n_bins_e, ax, species, color, marker)
        pd_eta["series"].append(series)

    ax.legend(fontsize=fs_lg)
    plt.tight_layout()
    path = output_dir / f"matching_efficiency_eta_{dataset_name}.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    save_plot_data_json(path, pd_eta)
    plt.close(fig)
    print(f"  Saved: {path}")
