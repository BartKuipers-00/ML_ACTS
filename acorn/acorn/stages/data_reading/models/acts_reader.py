import os
import logging
import pandas as pd
import numpy as np
import torch
from torch.utils.data import random_split

from ..data_reading_stage import EventReader
from . import trackml_utils


class ActsReader(EventReader):
    def __init__(self, config):
        super().__init__(config)
        """
        Here we initialize and load any attributes that are needed for the _build_single_csv function.
        """

        # ACTS allows to configure this, thus it can in vary
        if "simhit_stem" in config:
            self.simhit_stem = config["simhit_stem"]
        else:
            self.simhit_stem = "truth"

        # Determine if we try to add cell information
        cell_features = [
            "cell_count",
            "cell_val",
            "leta",
            "lphi",
            "lx",
            "ly",
            "lz",
            "geta",
            "gphi",
        ]
        self.use_cell_information = any(
            [
                feat in self.config["feature_sets"]["hit_features"]
                for feat in cell_features
            ]
        )

        # Import files
        input_dir = self.config["input_dir"]
        files = [
            self.simhit_stem,
            "measurements",
            "particles",
            "measurement-simhit-map",
        ]

        if self.use_cell_information:
            files.append("cells")

        self.raw_events = self.get_file_names(input_dir, filename_terms=files)

        assert type(self.config["data_split"]) == list
        assert len(self.config["data_split"]) == 3

        num_events = sum(self.config["data_split"])
        assert num_events <= len(self.raw_events)

        # Set random seed for reproducible splits (default: 42)
        split_seed = self.config.get("split_seed", 42)
        torch.manual_seed(split_seed)
        
        self.trainset, self.valset, self.testset = random_split(
            self.raw_events[:num_events], self.config["data_split"]
        )

        # Let's get the detector
        self.detector, self.detector_dims = trackml_utils.load_detector(
            self.config["detector_path"]
        )

        # Extract translation information additionally
        self.detector_dims["translations"] = np.zeros(
            (*trackml_utils.determine_array_size(self.detector), 3)
        )

        for i, r in self.detector.iterrows():
            v, l, m = tuple(map(int, (r.volume_id, r.layer_id, r.module_id)))
            self.detector_dims["translations"][v, l, m] = np.array(
                [r.cx.item(), r.cy.item(), r.cz.item()]
            )

    def _build_single_csv(self, event, output_dir=None):
        # Check if file already exists
        event_id = event["event_id"]
        particles_file_exists = os.path.exists(
            os.path.join(output_dir, "event{:09}-particles.csv".format(int(event_id)))
        )
        truth_file_exists = os.path.exists(
            os.path.join(
                output_dir, "event{:09}-{}.csv".format(int(event_id), self.simhit_stem)
            )
        )
        if particles_file_exists and truth_file_exists:
            logging.info(f"File {event_id} already exists, skipping...")
            return

        # Load CSV files
        simhits = pd.read_csv(event[self.simhit_stem])
        particles = pd.read_csv(event["particles"])
        measurements = pd.read_csv(event["measurements"])
        measurement_simhit_map = pd.read_csv(event["measurement-simhit-map"])

        # Process hits
        if self.config["use_truth_hits"]:
            truth = self._process_true_hits(simhits)
        else:
            truth = self._process_measurements(
                measurements, simhits, measurement_simhit_map
            )

        # Adjust column names so we can use trackml utils
        if self.use_cell_information:
            cells = pd.read_csv(event["cells"]).rename(
                {"channel0": "ch0", "channel1": "ch1"}, axis=1
            )
            assert len(cells) > 0

            truth = truth[truth.measurement_id.isin(np.unique(cells.hit_id))]

            print(len(truth), len(np.unique(cells.hit_id)))

            assert len(np.unique(cells.hit_id)) == len(truth)
            truth = trackml_utils.add_cell_info(truth, cells, self.detector_dims)

        # Ensure DataFrame has proper index before adding region labels
        # Reset index to ensure clean RangeIndex (fixes pandas assignment issues)
        if len(truth) > 0:
            truth = truth.reset_index(drop=True)
        else:
            logging.warning(f"Empty truth DataFrame for event {event_id}, skipping")
            return

        # Add module_id and region_id
        truth = trackml_utils.add_region_labels(truth, self.config["region_labels"])
        
        # Add required hit geometry columns for acorn
        # barrel volumes: 8, 13, 17 -> 0
        # endcap volumes: 7, 9, 12, 14, 16, 18 -> 2
        barrel_volumes = [8, 13, 17]
        truth['hit_barrel_endcap'] = truth['volume_id'].apply(
            lambda v: 0 if v in barrel_volumes else 2
        )
        
        # Add hit_layer_disk (use layer_id from detector geometry)
        truth['hit_layer_disk'] = truth['layer_id']
        
        # Add hit_eta_module and hit_phi_module (use module_id as placeholder)
        truth['hit_eta_module'] = truth['module_id'] 
        truth['hit_phi_module'] = truth['module_id']

        # Process particles
        particles = self._process_particles(particles, truth)

        # Save to CSV
        truth.to_csv(
            os.path.join(output_dir, "event{:09}-truth.csv".format(int(event_id))),
            index=False,
        )
        particles.to_csv(
            os.path.join(output_dir, "event{:09}-particles.csv".format(int(event_id))),
            index=False,
        )

    def _process_measurements(self, measurements, truth, simhit_map):
        """
        Create the hit information from the digitized measurements provided by ACTS
        """
        # Add decoded geometry information
        measurements = measurements.merge(
            self.detector[["geometry_id", "volume_id", "layer_id", "module_id"]],
            on="geometry_id",
        )

        # Add global positions
        measurements = self._measurements_add_global_pos(measurements)

        # Add hit and particle id.
        # simhit_map.hit_id is the SimHitContainer row position (== pandas row
        # index of the original, unmodified hits.csv).  clean_loops preserves
        # this as the 'simhit_id' column; fall back to pandas row index for
        # backwards-compatibility with raw ACTS output that has not been through
        # clean_loops yet.
        measurements["hit_id"] = measurements["measurement_id"].map(
            dict(zip(simhit_map.measurement_id, simhit_map.hit_id))
        )
        if "simhit_id" in truth.columns:
            simhit_id_to_particle = dict(zip(truth["simhit_id"], truth.particle_id))
        else:
            simhit_id_to_particle = dict(zip(truth.index, truth.particle_id))
        measurements["particle_id"] = measurements["hit_id"].map(simhit_id_to_particle)

        return measurements

    def _measurements_add_global_pos(self, measurements):
        """
        As the measruements are given in 2D surface-local coordinates, transform them to the global detector frame
        """
        global_pos = np.zeros((len(measurements.index), 3))

        for i, meas in enumerate(measurements.itertuples()):
            v, l, m = meas.volume_id, meas.layer_id, meas.module_id
            trans = self.detector_dims["translations"][v, l, m]
            rot = self.detector_dims["rotations"][v, l, m]

            x = np.array([meas.local0, meas.local1, 0.0])
            global_pos[i] = trans + (x @ rot.T)

        measurements[["x", "y", "z"]] = global_pos
        return measurements

    def _process_true_hits(self, truth):
        """
        Creates the hit information from the true simulated hits by ACTS
        """
        # Add decoded geometry information
        truth = truth.merge(
            self.detector[["geometry_id", "volume_id", "layer_id", "module_id"]],
            on="geometry_id",
        )

        # Drop index axis as this has a different meaning in this context and confuses pandas
        truth = truth.drop("index", axis=1)

        # Rename to hits scheme
        truth = truth.rename(columns={"tx": "x", "ty": "y", "tz": "z", "tt": "t"})

        # Assign a hit id
        truth["hit_id"] = np.arange(len(truth))
        truth["measurement_id"] = truth["hit_id"]

        return truth

    def _process_particles(self, particles, hits):
        """
        Calculate addtional information per particle
        """

        # Calculate the radius of the particle
        particles["radius"] = np.sqrt(particles["vx"] ** 2 + particles["vy"] ** 2)

        # Calculate the pT of the particle
        particles["pt"] = np.sqrt(particles["px"] ** 2 + particles["py"] ** 2)

        # Calculate the eta of the particle
        particles["particle_eta"] = np.arctanh(
            particles["pz"] / np.sqrt(particles["pt"] ** 2 + particles["pz"] ** 2)
        )

        # Add nhits information
        particles["nhits"] = (
            particles["particle_id"]
            .map(hits.groupby("particle_id").size())
            .fillna(0)
            .astype(int)
        )

        return particles
