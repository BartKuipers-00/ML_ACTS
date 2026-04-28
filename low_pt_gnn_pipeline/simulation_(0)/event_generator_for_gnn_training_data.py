#!/usr/bin/env python3
"""
Script to generate training data for GNN using ACTS EventGenerator.
This creates CSV files compatible with acorn's ActsReader.

Configuration is loaded from acorn_configs/simulation_(0)/event_generator_simulation.yaml

"""

from pathlib import Path
import yaml
import acts
import acts.examples
import acts.examples.hepmc3  
from acts import UnitConstants as u
from acts.examples import GenericDetector

# Import simulation functions
from acts.examples.simulation import (
    addFatras,
    addDigitization,
)


def load_simulation_config(config_path=None):
    """Load ACTS simulation configuration from YAML file"""
    if config_path is None:
        script_dir = Path(__file__).resolve().parent
        pipeline_root = script_dir.parent
        config_path = pipeline_root / "acorn_configs" / "simulation_(0)" / "event_generator_simulation.yaml"

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    pipeline_root = Path(__file__).resolve().parent.parent
    if 'output' in config and 'base_dir' in config['output']:
        if not Path(config['output']['base_dir']).is_absolute():
            config['output']['base_dir'] = str(pipeline_root / config['output']['base_dir'])

    return config


def get_particle_type(particle_name):
    """Map particle name to ACTS PDG particle type"""
    particle_map = {
        'muon': acts.PdgParticle.eMuon,
        'mu': acts.PdgParticle.eMuon,
        'pion': acts.PdgParticle.ePionPlus,
        'pi': acts.PdgParticle.ePionPlus,
        'electron': acts.PdgParticle.eElectron,
        'e': acts.PdgParticle.eElectron,
        'proton': acts.PdgParticle.eProton,
        'p': acts.PdgParticle.eProton,
        'kaon': acts.PdgParticle.eKaonPlus,
        'k': acts.PdgParticle.eKaonPlus,
        'photon': acts.PdgParticle.eGamma,
        'gamma': acts.PdgParticle.eGamma,
    }

    particle_name_lower = particle_name.lower()
    if particle_name_lower not in particle_map:
        raise ValueError(
            f"Unknown particle type: {particle_name}. "
            f"Available types: {', '.join(particle_map.keys())}"
        )

    return particle_map[particle_name_lower]


def get_log_level(level_name):
    """Map log level name to ACTS logging level"""
    level_map = {
        'VERBOSE': acts.logging.VERBOSE,
        'DEBUG': acts.logging.DEBUG,
        'INFO': acts.logging.INFO,
        'WARNING': acts.logging.WARNING,
        'ERROR': acts.logging.ERROR,
        'FATAL': acts.logging.FATAL,
    }
    return level_map.get(level_name.upper(), acts.logging.INFO)


def generate_minimal_training_data(config=None):

    # Load configuration
    if config is None:
        config = load_simulation_config()

    # Extract configuration values
    num_events = config['num_events']
    particles_per_vertex = config['particles_per_vertex']
    multiplicity = config['multiplicity']
    random_seed = config['random_seed']

    # Extract and parse particle species (NEW: supports list or single)
    gun_config = config['particle_gun']

    if 'particle_types' in gun_config:
        species_list = gun_config['particle_types']
        if not isinstance(species_list, list):
            raise ValueError("particle_types must be a list of particle names")
    elif 'particle_type' in gun_config:
        # Backward compatible: treat single species as list of one
        species_list = [gun_config['particle_type']]
    else:
        raise ValueError(
            "Config must specify either 'particle_type' (single species) "
            "or 'particle_types' (multiple species)"
        )

    # Convert species names to PDG particle types
    particle_pdg_list = [get_particle_type(species) for species in species_list]

    # Setup output directory
    output_dir = Path(config['output']['base_dir'])
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_dir = output_dir / "csv"
    csv_dir.mkdir(exist_ok=True)

    # Enhanced logging showing species information
    print("="*80)
    print("ACTS DATA GENERATION (Multi-Species EventGenerator)")
    print("="*80)
    print(f"Configuration: acorn_configs/simulation_(0)/event_generator_simulation.yaml")
    print(f"Events: {num_events}")
    print(f"Particle species: {', '.join(species_list)} ({len(species_list)} species)")
    print(f"Particles per species per vertex: {particles_per_vertex}")
    print(f"Multiplicity (vertices per event): {multiplicity}")
    total_particles = particles_per_vertex * len(species_list) * multiplicity
    print(f"Total particles per event: {total_particles} ({particles_per_vertex}×{len(species_list)}×{multiplicity})")
    print(f"pT range: {gun_config['momentum']['min']}-{gun_config['momentum']['max']} GeV")
    print(f"Eta range: {gun_config['eta']['min']} to {gun_config['eta']['max']}")
    print(f"Output directory: {output_dir}")
    print("="*80 + "\n")

    # Setup detector
    detector = GenericDetector()
    trackingGeometry = detector.trackingGeometry()

    # Setup magnetic field (solenoid field along z-axis)
    field_strength = config['magnetic_field']['strength']
    field = acts.ConstantBField(acts.Vector3(0, 0, field_strength * u.T))

    # Create sequencer
    sim_config = config['simulation']
    log_level = get_log_level(sim_config['log_level'])
    acts.logging.setLevel(log_level)

    # Mask FLTUND (underflow) FPEs from Bethe-Bloch energy loss for nearly-stopped particles.
    # Empty string for file matches all files (ends_with("") is always true).
    fpe_masks = [
        acts.examples.Sequencer.FpeMask(
            "", (0, 999999), acts.FpeType.FLTUND, 1000
        )
    ]

    s = acts.examples.Sequencer(
        events=num_events,
        numThreads=sim_config['num_threads'],
        logLevel=log_level,
        fpeMasks=fpe_masks,
    )

    # Add random number generator
    rnd = acts.examples.RandomNumbers(seed=random_seed)

    # Add EventGenerator 
    print("Adding EventGenerator...")
    print(f"  Species: {', '.join(species_list)}")

    mom_config = gun_config['momentum']
    eta_config = gun_config['eta']
    phi_config = gun_config['phi']

    # Create EventGenerator with one Generator per species
    # Pattern from acts/CI/physmon/workflows/physmon_simulation.py
    evGen = acts.examples.EventGenerator(
        level=get_log_level(sim_config['log_level']),
        generators=[
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=multiplicity),
                vertex=acts.examples.GaussianVertexGenerator(
                    mean=acts.Vector4(0, 0, 0, 0),
                    stddev=acts.Vector4(0, 0, 0, 0),
                ),
                particles=acts.examples.ParametricParticleGenerator(
                    p=(mom_config['min'] * u.GeV, mom_config['max'] * u.GeV),
                    pTransverse=mom_config['transverse'],
                    eta=(eta_config['min'], eta_config['max']),
                    phi=(phi_config['min'] * u.degree, phi_config['max'] * u.degree),
                    etaUniform=eta_config['uniform'],
                    numParticles=particles_per_vertex,
                    pdg=pdg,
                    randomizeCharge=gun_config['randomize_charge'],
                ),
            )
            for pdg in particle_pdg_list
        ],
        randomNumbers=rnd,
        outputEvent="particle_gun_event",
    )
    s.addReader(evGen)

    # Convert EventGenerator output to particles format expected by FATRAS
    hepmc3Converter = acts.examples.hepmc3.HepMC3InputConverter(
        level=get_log_level(sim_config['log_level']),
        inputEvent=evGen.config.outputEvent,
        outputParticles="particles_generated",
        outputVertices="vertices_generated",
        mergePrimaries=False,
    )
    s.addAlgorithm(hepmc3Converter)

    # Set whiteboard aliases for downstream algorithms
    s.addWhiteboardAlias("particles", hepmc3Converter.config.outputParticles)
    s.addWhiteboardAlias("vertices_truth", hepmc3Converter.config.outputVertices)

    # Alias for FATRAS (expects particles_generated_selected)
    s.addWhiteboardAlias(
        "particles_generated_selected", hepmc3Converter.config.outputParticles
    )

    print(f"  ✓ Configured {len(species_list)} particle species")

    # Write initial particles to CSV (optional)
    if config['output']['write_particles_initial']:
        s.addWriter(
            acts.examples.CsvParticleWriter(
                level=get_log_level(sim_config['log_level']),
                inputParticles="particles_generated",
                outputDir=str(csv_dir),
                outputStem="particles_initial",
            )
        )

    # Add FATRAS simulation
    print("Adding FATRAS simulation...")
    p_min = sim_config.get('p_min', 0.01)
    if p_min is not None:
        p_min = p_min * u.GeV  # Convert to ACTS units

    addFatras(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        outputDirCsv=csv_dir,
        enableInteractions=sim_config['enable_interactions'],
        pMin=p_min,
        loopFraction=sim_config.get('loop_fraction', None),
        maxSteps=sim_config.get('max_steps', None),
    )

    # Add digitization
    print("Adding digitization...")
    digi_config = config['digitization']

    # ACTS is at /data/alice/bkuipers/acts/ (sibling of the pipeline folder)
    srcdir = Path(__file__).resolve().parent.parent.parent / "acts"
    digiConfigFile = srcdir / digi_config['config_file'] if digi_config['config_file'] else None

    if digiConfigFile and not digiConfigFile.exists():
        print(f"Warning: Config file not found at {digiConfigFile}")
        digiConfigFile = None

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        outputDirCsv=csv_dir,
        rnd=rnd,
    )

    # Write detector geometry
    if config['output']['write_geometry']:
        print("Writing detector geometry...")
        s.addWriter(
            acts.examples.CsvTrackingGeometryWriter(
                level=get_log_level(sim_config['log_level']),
                trackingGeometry=trackingGeometry,
                outputDir=str(csv_dir),
                writePerEvent=False,
            )
        )

    # Run the simulation
    print("Running simulation...")
    s.run()

    print(f"\n✓ Data generation complete!")
    print(f"Output files in: {csv_dir}")
    print("\nGenerated files:")
    for f in sorted(csv_dir.glob("*.csv")):
        print(f"  - {f.name}")

    return csv_dir


if __name__ == "__main__":
    # Load configuration from YAML file
    config = load_simulation_config()

    # Generate data using config file settings
    generate_minimal_training_data(config=config)
