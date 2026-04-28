#!/usr/bin/env python3
"""
Geant4 simulation for GNN training data (multi-species particle gun).

Config: acorn_configs/simulation_(0)/geant4_simulation.yaml



Downstream notes
----------------
- CSV format is identical (same ACTS writers) so clean_loops and
  convert_csv_to_pyg_sets.py work without changes.
- ODD has different volume_id values than GenericDetector, so you may want to
  update region_labels and detector_path in
  acorn_configs/latent_stage_(1)/convert_csv_to_pyg_sets.yaml when running
  on data produced by this script.
"""

from pathlib import Path
import sys
import yaml

# ---------------------------------------------------------------------------
# Ensure the acts_new build is on sys.path.
#
# This script needs acts_new (not the old acts):
#   - acts_new is built with Geant4 + DD4HEP/ODD support
#   - acts (used by the rest of the pipeline) lacks Geant4/ODD bindings
#
# How to run:
#   source /data/alice/bkuipers/acts_new/this_acts_withdeps.sh
#   python geant4_event_generator.py
#
# The source command sets PYTHONPATH so that "import acts" picks up acts_new.
# The block below inserts acts_new/build/python at the front of sys.path as a
# fallback, so the script is also self-contained when run directly.
# ---------------------------------------------------------------------------
_acts_new_python = (
    Path(__file__).resolve().parent  # simulation_(0)/
    .parent                          # low_pt_gnn_pipeline/
    .parent                          # /data/alice/bkuipers/
    / "acts_new" / "build" / "python"
)
if _acts_new_python.is_dir():
    sys.path.insert(0, str(_acts_new_python))

import acts
import acts.examples
import acts.examples.hepmc3
from acts import UnitConstants as u
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory
from acts.examples.simulation import (
    addGeant4,
    addDigitization,
)


# ---------------------------------------------------------------------------
# Configuration helpers
# ---------------------------------------------------------------------------

def load_config(config_path=None):
    """Load simulation config from YAML (default: geant4_simulation.yaml)."""
    if config_path is None:
        script_dir = Path(__file__).resolve().parent
        pipeline_root = script_dir.parent
        config_path = (
            pipeline_root
            / "acorn_configs"
            / "simulation_(0)"
            / "geant4_simulation.yaml"
        )

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # Make output base_dir absolute relative to pipeline root
    pipeline_root = Path(__file__).resolve().parent.parent
    if "output" in config and "base_dir" in config["output"]:
        if not Path(config["output"]["base_dir"]).is_absolute():
            config["output"]["base_dir"] = str(
                pipeline_root / config["output"]["base_dir"]
            )

    return config


def get_particle_type(particle_name):
    """Map particle name string to ACTS PDG particle type."""
    particle_map = {
        "muon":     acts.PdgParticle.eMuon,
        "mu":       acts.PdgParticle.eMuon,
        "pion":     acts.PdgParticle.ePionPlus,
        "pi":       acts.PdgParticle.ePionPlus,
        "electron": acts.PdgParticle.eElectron,
        "e":        acts.PdgParticle.eElectron,
        "proton":   acts.PdgParticle.eProton,
        "p":        acts.PdgParticle.eProton,
        "kaon":     acts.PdgParticle.eKaonPlus,
        "k":        acts.PdgParticle.eKaonPlus,
        "photon":   acts.PdgParticle.eGamma,
        "gamma":    acts.PdgParticle.eGamma,
    }
    key = particle_name.lower()
    if key not in particle_map:
        raise ValueError(
            f"Unknown particle type: {particle_name!r}. "
            f"Available: {', '.join(particle_map)}"
        )
    return particle_map[key]


def get_log_level(level_name):
    """Map log-level string to ACTS logging level."""
    return {
        "VERBOSE": acts.logging.VERBOSE,
        "DEBUG":   acts.logging.DEBUG,
        "INFO":    acts.logging.INFO,
        "WARNING": acts.logging.WARNING,
        "ERROR":   acts.logging.ERROR,
        "FATAL":   acts.logging.FATAL,
    }.get(level_name.upper(), acts.logging.INFO)


# ---------------------------------------------------------------------------
# Main simulation function
# ---------------------------------------------------------------------------

def run_geant4_simulation(config=None):
    if config is None:
        config = load_config()

    # ---- top-level config values ----
    num_events          = config["num_events"]
    particles_per_vertex = config["particles_per_vertex"]
    multiplicity        = config["multiplicity"]
    random_seed         = config["random_seed"]

    # ---- particle species ----
    gun_config = config["particle_gun"]
    if "particle_types" in gun_config:
        species_list = gun_config["particle_types"]
        if not isinstance(species_list, list):
            raise ValueError("particle_types must be a list of particle name strings")
    elif "particle_type" in gun_config:
        species_list = [gun_config["particle_type"]]
    else:
        raise ValueError(
            "Config must have either 'particle_type' (single) or "
            "'particle_types' (multiple) under particle_gun"
        )
    particle_pdg_list = [get_particle_type(s) for s in species_list]

    # ---- output directory ----
    output_dir = Path(config["output"]["base_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_dir = output_dir / "csv"
    csv_dir.mkdir(exist_ok=True)

    # ---- simulation settings ----
    sim_config      = config["simulation"]
    log_level       = get_log_level(sim_config["log_level"])
    kill_after_time = sim_config.get("kill_after_time_ns", 40.0)  # ns

    # ---- kinematic ranges ----
    mom_config = gun_config["momentum"]
    eta_config = gun_config["eta"]
    phi_config = gun_config["phi"]

    # ---- print banner ----
    print("=" * 80)
    print("ACTS DATA GENERATION — Geant4 + OpenDataDetector")
    print("=" * 80)
    print(f"Config:              acorn_configs/simulation_(0)/geant4_simulation.yaml")
    print(f"Events:              {num_events}")
    print(f"Particle species:    {', '.join(species_list)}")
    print(f"Particles/vertex:    {particles_per_vertex} per species")
    print(f"Multiplicity:        {multiplicity} vertex/event")
    total = particles_per_vertex * len(species_list) * multiplicity
    print(f"Total/event:         {total}  ({particles_per_vertex}×{len(species_list)}×{multiplicity})")
    print(f"pT range:            {mom_config['min']}–{mom_config['max']} GeV")
    print(f"η range:             {eta_config['min']} to {eta_config['max']}")
    print(f"killAfterTime:       {kill_after_time} ns")
    print(f"Output:              {output_dir}")
    print("=" * 80 + "\n")

    # ---- detector (ODD — required for Geant4) ----
    geo_dir = getOpenDataDetectorDirectory()

    # Material map: use configured path or fall back to ODD default
    mat_map_cfg = config.get("detector", {}).get("material_map", None)
    if mat_map_cfg:
        mat_map_path = Path(mat_map_cfg)
    else:
        mat_map_path = geo_dir / "data" / "odd-material-maps.root"

    if mat_map_path.exists():
        material_deco = acts.IMaterialDecorator.fromFile(mat_map_path)
    else:
        print(f"Warning: material map not found at {mat_map_path}, running without material effects")
        material_deco = None

    detector_kwargs = {}
    if material_deco is not None:
        detector_kwargs["materialDecorator"] = material_deco

    detector         = getOpenDataDetector(odd_dir=geo_dir, **detector_kwargs)
    tracking_geometry = detector.trackingGeometry()

    # ---- magnetic field ----
    field_strength = config.get("magnetic_field", {}).get("strength", 2.0)
    field = acts.ConstantBField(acts.Vector3(0, 0, field_strength * u.T))

    # ---- sequencer ----
    s = acts.examples.Sequencer(
        events=num_events,
        numThreads=1,       # Geant4 requires single-threaded
        logLevel=log_level,
    )

    rnd = acts.examples.RandomNumbers(seed=random_seed)

    # ---- multi-species particle gun via EventGenerator ----
    # One Generator per species, each producing `particles_per_vertex` particles
    # at `multiplicity` vertices.  This mirrors the FATRAS event_generator exactly.
    print(f"Setting up particle gun for species: {', '.join(species_list)}")

    evGen = acts.examples.EventGenerator(
        level=log_level,
        generators=[
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=multiplicity),
                vertex=acts.examples.GaussianVertexGenerator(
                    mean=acts.Vector4(0, 0, 0, 0),
                    stddev=acts.Vector4(0, 0, 0, 0),
                ),
                particles=acts.examples.ParametricParticleGenerator(
                    p=(mom_config["min"] * u.GeV, mom_config["max"] * u.GeV),
                    pTransverse=mom_config["transverse"],
                    eta=(eta_config["min"], eta_config["max"]),
                    phi=(phi_config["min"] * u.degree, phi_config["max"] * u.degree),
                    etaUniform=eta_config.get("uniform", True),
                    numParticles=particles_per_vertex,
                    pdg=pdg,
                    randomizeCharge=gun_config["randomize_charge"],
                ),
            )
            for pdg in particle_pdg_list
        ],
        randomNumbers=rnd,
        outputEvent="particle_gun_event",
    )
    s.addReader(evGen)

    # Convert HepMC3 event to particle collection on the whiteboard
    hepmc3_converter = acts.examples.hepmc3.HepMC3InputConverter(
        level=log_level,
        inputEvent=evGen.config.outputEvent,
        outputParticles="particles_generated",
        outputVertices="vertices_generated",
        mergePrimaries=False,
    )
    s.addAlgorithm(hepmc3_converter)

    # addGeant4 reads inputParticles="particles_generated_selected" by default.
    # addParticleGun also sets "particles_selected", so we mirror both.
    s.addWhiteboardAlias("particles_generated_selected", hepmc3_converter.config.outputParticles)
    s.addWhiteboardAlias("particles_selected",           hepmc3_converter.config.outputParticles)
    s.addWhiteboardAlias("vertices_truth",               hepmc3_converter.config.outputVertices)
    print(f"  ✓ {len(species_list)} species configured\n")

    # ---- Geant4 simulation ----
    print("Adding Geant4 simulation...")
    addGeant4(
        s,
        detector,
        tracking_geometry,
        field,
        rnd=rnd,
        outputDirCsv=csv_dir,
        outputDirRoot=None,   # set to outputDir if you also want ROOT output
        outputDirObj=None,
        killVolume=tracking_geometry.highestTrackingVolume,
        killAfterTime=kill_after_time * u.ns,
    )

    # ---- digitization ----
    print("Adding digitization...")
    digi_config = config.get("digitization", {})
    digi_config_file_rel = digi_config.get(
        "config_file", "Examples/Configs/odd-digi-smearing-config.json"
    )

    # Try to locate digi config relative to acts_new (sibling of pipeline root)
    acts_new_dir = Path(__file__).resolve().parent.parent.parent / "acts_new"
    digi_config_path = acts_new_dir / digi_config_file_rel

    if not digi_config_path.exists():
        # Fall back to looking relative to the ODD geo directory
        alt = geo_dir.parent / digi_config_file_rel
        if alt.exists():
            digi_config_path = alt
        else:
            print(f"Warning: digitization config not found at {digi_config_path}")
            digi_config_path = None

    addDigitization(
        s,
        tracking_geometry,
        field,
        digiConfigFile=digi_config_path,
        outputDirCsv=csv_dir,
        rnd=rnd,
    )

    # ---- detector geometry CSV (for downstream ACORN reader) ----
    if config.get("output", {}).get("write_geometry", True):
        print("Writing detector geometry CSV...")
        s.addWriter(
            acts.examples.CsvTrackingGeometryWriter(
                level=log_level,
                trackingGeometry=tracking_geometry,
                outputDir=str(csv_dir),
                writePerEvent=False,
            )
        )

    # ---- run ----
    print("Running Geant4 simulation...\n")
    s.run()

    print(f"\n✓ Simulation complete!")
    print(f"Output files in: {csv_dir}")
    print("\nGenerated files:")
    for f in sorted(csv_dir.glob("*.csv"))[:20]:
        print(f"  - {f.name}")
    csv_files = list(csv_dir.glob("*.csv"))
    if len(csv_files) > 20:
        print(f"  ... and {len(csv_files) - 20} more")

    return csv_dir


if __name__ == "__main__":
    config = load_config()
    run_geant4_simulation(config)
