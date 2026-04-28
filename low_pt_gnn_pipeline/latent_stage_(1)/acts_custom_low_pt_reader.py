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
import logging
import os
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
        """
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

    def _build_true_tracks(self, hits):
        """
        Override to create TRAJECTORY-ORDERED SEQUENTIAL EDGES using pre-computed segments.

        Requires hit_segment_id to already be present in the hits DataFrame,
        assigned by clean_loops_and_attribute_segments.py in the simulation stage.

        Edges are built per (particle, segment) — no edges cross segment boundaries.
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

        # Filter signal hits and sort by trajectory order (particle_id, then time)
        signal = hits[(hits.hit_particle_id != 0)].copy()
        signal = signal.sort_values(["hit_particle_id", time_col]).reset_index(drop=False)

        # --- Build edges per (particle, segment) — no edges cross segment boundaries ---
        signal_index_list = (
            signal.groupby(["hit_particle_id", "hit_segment_id"], sort=False)["index"]
            .agg(lambda x: list(x))
        )

        track_index_edges = []
        for segment_hits in signal_index_list.values:
            # Create sequential chain: hit[0]→hit[1], hit[1]→hit[2], ..., hit[n-2]→hit[n-1]
            if len(segment_hits) >= 2:
                for i in range(len(segment_hits) - 1):
                    track_index_edges.append((segment_hits[i], segment_hits[i + 1]))

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
