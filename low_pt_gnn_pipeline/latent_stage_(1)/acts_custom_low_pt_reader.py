"""
Custom ACTS Reader for Low-pT Data with Time-Based Trajectory Ordering

This reader extends ActsReader to create trajectory-ordered sequential edges
instead of the layer-based ordering (not suitable for loops).

Expects hit_segment_id to be pre-computed in the CSV by
simulation_(0)/clean_loops_and_attribute_segments.py, which assigns radial
segments (outgoing=1, incoming=2, outgoing=3, ...) and enforces loop fraction.

Edges are built per (particle, segment) — no edges cross segment boundaries.

Usage:
    In convert_csv_to_pyg_sets.yaml, set:
        model: ActsCustomLowPTReader
"""

import sys
from pathlib import Path

# Add acorn to path
SCRIPT_DIR = Path(__file__).resolve().parent
PIPELINE_ROOT = SCRIPT_DIR.parent
WORKSPACE_ROOT = PIPELINE_ROOT.parent
sys.path.insert(0, str(WORKSPACE_ROOT / 'acorn'))
sys.path.insert(0, str(PIPELINE_ROOT))

from acorn.stages.data_reading.models.acts_reader import ActsReader
from low_pt_custom_utils.region_label_utils import resolve_region_labels
import logging
import os
from itertools import product
import numpy as np
import pandas as pd
import torch
from torch.utils.data import random_split


class ActsCustomLowPTReader(ActsReader):
    """
    Custom reader for low-pT data that uses time-based trajectory ordering.

    For looping particles, we need to order hits by their actual trajectory order
    (using time), not by layer/module.
    """

    def __init__(self, config):
        """
        Override to support event_range parameter for parallel processing.

        If config contains 'event_range': [start, end], only process events
        in that range from the raw CSV files.
        """
        # Select the detector region map (config['detector'] = 'generic' | 'odd').
        # No auto-detection: prompts interactively / errors in batch if unset.
        resolve_region_labels(config, config.get("detector"))
        super().__init__(config)
        # Apply event_range if specified (for parallel chunk processing)
        if 'event_range' in config:
            start, end = config['event_range']
            print(f"Applying event_range filter: [{start}, {end})")
            print(f"Before: {len(self.raw_events)} events available")
            self.raw_events = self.raw_events[start:end]                       # Slice raw_events to only include the requested range
            print(f"After: {len(self.raw_events)} events in this chunk")
            num_events = sum(self.config["data_split"])
            assert num_events <= len(self.raw_events), f"Requested {num_events} events but only {len(self.raw_events)} available in chunk"

            self.trainset = self.raw_events[:num_events]
            self.valset = []
            self.testset = []

            print(f"Processing {len(self.trainset)} events sequentially (no split)")
            print()

    def _process_particles(self, particles, hits):
        """
        Override to add a `primary` flag matching ACTS exactly.

        ACTS defines secondaries in ActsFatras::Particle::isSecondary()
        (Fatras/include/ActsFatras/EventData/Particle.hpp):

            secondary  <=>  vertexSecondary != 0
                            OR generation     != 0
                            OR subParticle    != 0

        i.e. a particle is PRIMARY iff all three barcode sub-fields are zero.
        In the ACTS CSV these are particle_id_sv / particle_id_gen /
        particle_id_subpart. `primary` is the logical negation (1 = primary,
        0 = secondary), stored as int so it survives as a numeric graph feature.

        Note ACTS keeps secondaries by default (ParticleSelector.removeSecondaries
        defaults to false); this only labels them so downstream stages can
        weight/cut on `track_particle_primary` if desired.
        """
        particles = super()._process_particles(particles, hits)

        barcode_subfields = ["particle_id_sv", "particle_id_gen", "particle_id_subpart"]
        if all(col in particles.columns for col in barcode_subfields):
            is_secondary = (
                (particles["particle_id_sv"] != 0)
                | (particles["particle_id_gen"] != 0)
                | (particles["particle_id_subpart"] != 0)
            )
            particles["primary"] = (~is_secondary).astype(int)
        else:
            # Barcode sub-fields unavailable (e.g. non-ACTS input); label all primary.
            print(
                "WARNING: barcode sub-fields not found; setting primary=1 for all "
                "particles (cannot distinguish secondaries)."
            )
            particles["primary"] = 1

        return particles

    def _process_measurements(self, measurements, simhits, simhit_map):
        """
        Override to propagate segment_id and time from simhits to measurements.

        The parent only joins geometry, hit_id, and particle_id. We additionally
        map segment_id and tt (renamed to t) so the rest of the pipeline sees the
        same column names as the truth-hits path.

        simhit_map.hit_id is the SimHitContainer row position (original pandas
        row index of hits.csv before any rows were removed).  clean_loops
        preserves this as the 'simhit_id' column; fall back to pandas row index
        for backwards-compatibility with raw ACTS output.

        If config['geometry_from_csv'] is set (e.g. custom detectors like ALICE 3
        with no detectors.csv), geometry columns (volume_id/layer_id/module_id) and
        global positions (global_x/y/z) are taken straight from the measurements CSV
        instead of being merged/recomputed from self.detector.
        """
        if self.config.get("geometry_from_csv", False):
            result = self._process_measurements_from_csv(measurements, simhits, simhit_map)
        else:
            result = super()._process_measurements(measurements, simhits, simhit_map)

        # Build simhit_id -> tt / segment_id lookups using the stable SimHitContainer
        # position stored in the 'simhit_id' column (added by clean_loops).
        if "simhit_id" in simhits.columns:
            key = simhits["simhit_id"]
        else:
            key = simhits.index  # raw ACTS output: row index == container position

        result["t"] = result["hit_id"].map(dict(zip(key, simhits["tt"])))
        result["segment_id"] = result["hit_id"].map(dict(zip(key, simhits["segment_id"])))

        return result

    def _process_measurements_from_csv(self, measurements, simhits, simhit_map):
        """Geometry + global positions already in the measurements CSV (no detectors.csv).

        Mirrors the parent's particle/hit-id join but skips the detector merge and the
        local->global transform, using the CSV columns directly.
        """
        m = measurements.rename(
            columns={"global_x": "x", "global_y": "y", "global_z": "z"}
        ).copy()
        m["hit_id"] = m["measurement_id"].map(
            dict(zip(simhit_map.measurement_id, simhit_map.hit_id))
        )
        if "simhit_id" in simhits.columns:
            key = simhits["simhit_id"]
        else:
            key = simhits.index
        m["particle_id"] = m["hit_id"].map(dict(zip(key, simhits.particle_id)))
        return m

    def _build_true_tracks(self, hits):
        """
        Override to create TRAJECTORY-ORDERED, LAYER-TO-LAYER EDGES using
        pre-computed segments.

        This reproduces the stock ACORN ground-truth convention (connect every
        hit in a layer to every hit in the next layer, no in-layer edges) but
        replaces the "radially out from production" ordering proxy — which
        zigzags for looping low-pT particles — with TIME ordering within each
        radial segment.

        Per (particle, segment):
          1. Order hits by time (= trajectory order).
          2. Group consecutive hits into detector layers; because hits are
             time-sorted and grouped with sort=False, layers appear in
             trajectory order.
          3. Connect every hit in layer k to every hit in layer k+1
             (Cartesian product). Hits in the SAME layer are NOT connected.

        No edges cross segment boundaries.

        Requires hit_segment_id (from clean_loops_and_attribute_segments.py) and
        layer-identifying columns (hit_volume_id + hit_layer_id by default) to be
        present in the hits DataFrame.
        """
        # Verify required columns are present
        required_cols = ["hit_particle_id", "hit_id"]
        assert all(col in hits.columns for col in required_cols), \
            f"Missing required columns. Need: {required_cols}"

        time_col = "hit_t"
        assert time_col in hits.columns, \
            f"Time column '{time_col}' not found (check diffrent naming? tt, t? )"

        assert "hit_segment_id" in hits.columns, \
            "hit_segment_id not found. Run clean_loops_and_attribute_segments.py first."

        # Columns that uniquely identify a detector layer. Override via the
        # 'layer_columns' config key if needed. Default to (volume, layer); fall
        # back to (region, layer_disk) if the raw geometry ids are unavailable.
        layer_columns = self.config.get("layer_columns")
        if layer_columns is None:
            if {"hit_volume_id", "hit_layer_id"}.issubset(hits.columns):
                layer_columns = ["hit_volume_id", "hit_layer_id"]
            elif {"hit_region", "hit_layer_disk"}.issubset(hits.columns):
                layer_columns = ["hit_region", "hit_layer_disk"]
            else:
                raise AssertionError(
                    "No layer-identifying columns found. Expected hit_volume_id+"
                    "hit_layer_id or hit_region+hit_layer_disk, or set "
                    "'layer_columns' in the config."
                )
        assert all(c in hits.columns for c in layer_columns), \
            f"layer_columns {layer_columns} not all present in hits."

        # Filter signal hits and sort by trajectory order (particle_id, then time)
        signal = hits[(hits.hit_particle_id != 0)].copy()
        signal = signal.sort_values(["hit_particle_id", time_col]).reset_index(drop=False)

        # --- Build layer-to-layer edges per (particle, segment) ---
        # First group by (particle, segment, layer): hits sharing a layer. Because
        # the frame is time-sorted and sort=False, layer groups appear in
        # trajectory order. Then regroup per (particle, segment) into an ordered
        # list of layers. No edges cross segment boundaries.
        layer_index_list = (
            signal.groupby(
                ["hit_particle_id", "hit_segment_id"] + layer_columns, sort=False
            )["index"]
            .agg(list)                                  # hits within one layer
            .groupby(level=[0, 1], sort=False)          # regroup per (particle, segment)
            .agg(list)                                  # ordered list of layers
        )

        track_index_edges = []
        for layers in layer_index_list.values:
            # Connect every hit in layer k to every hit in layer k+1 (all-to-all).
            # Hits within the same layer are intentionally NOT connected.
            for layer_a, layer_b in zip(layers[:-1], layers[1:]):
                track_index_edges.extend(product(layer_a, layer_b))

        if len(track_index_edges) == 0:
            return np.array([]), np.array([]), np.array([])

        # Convert to numpy array format [2, num_edges]
        track_index_edges = np.array(track_index_edges).T

        track_edges = hits.hit_id.values[track_index_edges]

        track_features = self._get_track_features(hits, track_index_edges, track_edges)

        # Remap
        track_edges, track_features, hits = self.remap_edges(
            track_edges, track_features, hits
        )

        return track_edges, track_features, hits
