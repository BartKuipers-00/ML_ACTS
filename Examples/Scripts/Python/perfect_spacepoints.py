#!/usr/bin/env python3

from pathlib import Path
from typing import Optional, Dict, Any
import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
import uproot
import awkward as ak

import acts
from acts import UnitConstants as u
from acts.examples import GenericDetector, RootParticleReader, RootParticleWriter


from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addFatras,
    addDigitization,
    ParticleSelectorConfig,
    addDigiParticleSelection,
)
from acts.examples.reconstruction import (
    addSeeding,
    TrackSmearingSigmas,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedingAlgorithm,
    TruthEstimatedSeedingAlgorithmConfigArg,
    addTrackWriters,
    TrackSelectorConfig,
    trackSelectorDefaultKWArgs,
)


PDG_MAP = {
    "electron": acts.PdgParticle.eElectron,
    "positron": acts.PdgParticle.ePositron,
    "muon": acts.PdgParticle.eMuon,
    "antimuon": acts.PdgParticle.eAntiMuon,
    "tau": acts.PdgParticle.eTau,
    "antitau": acts.PdgParticle.eAntiTau,
    "pionplus": acts.PdgParticle.ePionPlus,
    "pionminus": acts.PdgParticle.ePionMinus,
    "pion0": acts.PdgParticle.ePionZero,
    "kaonplus": acts.PdgParticle.eKaonPlus,
    "kaonminus": acts.PdgParticle.eKaonMinus,
    "proton": acts.PdgParticle.eProton,
    "antiproton": acts.PdgParticle.eAntiProton,
    "gamma": acts.PdgParticle.eGamma,          #note photon does not work yet (mass problem)
}


def load_config(config_file: Path) -> Dict[str, Any]:
    """Load configuration from JSON file."""
    with open(config_file, 'r') as f:
        return json.load(f)



def runPerfectSpacepoints(
    config: Dict[str, Any],
    trackingGeometry,
    decorators,
    geometrySelection: Path,
    digiConfigFile: Path,
    field,
    outputDir: Path,
    inputParticlePath: Optional[Path] = None,
    s=None,
    loop_protection: Optional[bool] = None,
    loop_fraction: Optional[float] = None,
):


    # Load configuration sections
    sim = config["simulation"]
    pg = sim["particleGun"]
    output = config["output"]
    logging_cfg = config["logging"]

    s = s or acts.examples.Sequencer(
        events=sim["events"],
        numThreads=sim["numThreads"],
        logLevel=getattr(acts.logging, logging_cfg["level"]),
    )
    for d in decorators:
        s.addContextDecorator(d)
    rnd = acts.examples.RandomNumbers(seed=sim["randomSeed"])
    outputDir = Path(outputDir)

    # Normalize particle types from config; accept single string, list, or stringified list
    import ast
    ptype_raw = pg["particles"]["type"]
    particle_types: list[str]
    if isinstance(ptype_raw, list):
        particle_types = [str(sp).lower().replace(" ", "") for sp in ptype_raw]
    elif isinstance(ptype_raw, str):
        s_pt = ptype_raw.strip()
        if (s_pt.startswith("[") and s_pt.endswith("]")) or ("," in s_pt):
            try:
                # Try JSON first
                parsed = json.loads(s_pt)
            except Exception:
                # Fallback to Python literal (handles single quotes)
                parsed = ast.literal_eval(s_pt)
            if not isinstance(parsed, list):
                parsed = [parsed]
            particle_types = [str(sp).lower().replace(" ", "") for sp in parsed]
        else:
            particle_types = [s_pt.lower().replace(" ", "")]
    else:
        particle_types = [str(ptype_raw).lower().replace(" ", "")]

    # Select source mode
    source_mode = "reader" if inputParticlePath is not None else ("multi" if len(particle_types) > 1 else "gun")

    # Build particle source
    customLogLevel = acts.examples.defaultLogging(
        s, getattr(acts.logging, logging_cfg["level"])
    )
    if source_mode == "multi":
        particles_cfg = pg["particles"]
        momentum = pg["momentum"]
        eta_cfg = pg["eta"]
        phi_cfg = pg["phi"]
        vtxGen = acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(0, 0, 0, 0),
        )
        #create one generator per particle type and combine them in one event generator later. Then pass to the sequencer
        generators = []
        for sp in particle_types:
            if sp not in PDG_MAP:
                raise ValueError(
                    f"Unknown particle type '{sp}'. Allowed: {', '.join(PDG_MAP.keys())}"
                )
            pdg_enum = PDG_MAP[sp]
            generators.append(
                acts.examples.EventGenerator.Generator(
                    multiplicity=acts.examples.FixedMultiplicityGenerator(n=pg["multiplicity"]),
                    vertex=vtxGen,
                    particles=acts.examples.ParametricParticleGenerator(
                        **acts.examples.defaultKWArgs(
                            p=(momentum["min"], momentum["max"]),
                            pTransverse=momentum["transverse"],
                            eta=(eta_cfg["min"], eta_cfg["max"]),
                            phi=(phi_cfg["min"], phi_cfg["max"]),
                            etaUniform=eta_cfg["uniform"],
                            numParticles=particles_cfg["count"],
                            pdg=pdg_enum,
                            randomizeCharge=particles_cfg["randomizeCharge"],
                        )
                    ),
                )
            )
        evGen = acts.examples.EventGenerator(
            level=customLogLevel(),
            generators=generators,
            randomNumbers=rnd,
            outputEvent="particle_gun_event",
        )
        s.addReader(evGen)
        hepmc3Converter = acts.examples.hepmc3.HepMC3InputConverter(
            level=customLogLevel(),
            inputEvent=evGen.config.outputEvent,
            outputParticles="particles_generated",
            outputVertices="vertices_generated",
            mergePrimaries=False,
        )
        s.addAlgorithm(hepmc3Converter)
        s.addWhiteboardAlias("particles", hepmc3Converter.config.outputParticles)
        s.addWhiteboardAlias("vertices_truth", hepmc3Converter.config.outputVertices)
        s.addWhiteboardAlias("particles_generated_selected", hepmc3Converter.config.outputParticles)

        # Write generated particles to ROOT file
        s.addWriter(
            RootParticleWriter(
                level=customLogLevel(),
                inputParticles=hepmc3Converter.config.outputParticles,
                filePath=str(outputDir / "particles.root"),
            )
        )

    # Add particle source: either read from ROOT file or generate via particle gun
    if inputParticlePath is not None:
        acts.logging.getLogger("PerfectSpacepointsExample").info(
            "Reading particles from %s", inputParticlePath.resolve()
        )
        assert inputParticlePath.exists(), f"Input file not found: {inputParticlePath}"
        s.addReader(
            RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                outputParticles="particles_generated",
            )
        )
        # Ensure downstream expects the selected alias
        s.addWhiteboardAlias(
            "particles_generated_selected", "particles_generated"
        )
    elif source_mode == "gun":
        # Single-species particle gun
        ptype_key = particle_types[0]
        if ptype_key not in PDG_MAP:
            raise ValueError(
                f"Unknown particle type '{ptype_key}'. Allowed: {', '.join(PDG_MAP.keys())}"
            )
        pdg_enum = PDG_MAP[ptype_key]
        momentum = pg["momentum"]
        eta = pg["eta"]
        phi = pg["phi"]
        particles_cfg = pg["particles"]
        addParticleGun(
            s,
            MomentumConfig(
                momentum["min"] * u.GeV,
                momentum["max"] * u.GeV,
                transverse=momentum["transverse"],
            ),
            EtaConfig(eta["min"], eta["max"], uniform=eta["uniform"]),
            PhiConfig(phi["min"] * u.degree, phi["max"] * u.degree),
            ParticleConfig(
                particles_cfg["count"],
                pdg_enum,
                randomizeCharge=particles_cfg["randomizeCharge"],
            ),
            multiplicity=pg["multiplicity"],
            rnd=rnd,
        )

        # Write generated particles to ROOT file
        s.addWriter(
            RootParticleWriter(
                level=customLogLevel(),
                inputParticles="particles_generated",
                filePath=str(outputDir / "particles.root"),
            )
        )
    # else: reader branch already handled

    # Physics simulation and digitization
    physics = config["simulation_physics"]
    selection = config["particleSelection"]
    pMin = physics["pMin"]
    if pMin is not None:
        pMin = pMin * u.GeV

    addFatras(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        enableInteractions=physics['enableInteractions'],
        pMin=pMin,
        # forward optional simulation hyperparameters from JSON if present
        maxSteps=sim.get("maxSteps", None),
        loopProtection=sim.get("loopProtection", None),
        loopFraction=sim.get("loopFraction", None),
        # pass debug interval from logging config (if present)
        debugStepInterval=logging_cfg.get("n", None),
        outputDirRoot=outputDir,
    )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
    )

    pt_sel = selection["pt"]
    measurements = selection["measurements"]
    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(pt_sel["min"] * u.GeV, None if pt_sel["max"] is None else pt_sel["max"] * u.GeV),
            measurements=(
                measurements["min"],
                None if measurements["max"] is None else measurements["max"],
            ),
            removeNeutral=selection["removeNeutral"],
        ),
    )

    # Seeding configuration
    seed = config["seeding"]
    smearing = seed["smearing"]
    finder = seed["finderConfig"]
    options = seed["finderOptions"]
    truth_est = seed["truthEstimated"]
    initial = seed["initialSigmas"]

    # Derive bFieldInZ from detector field and warn if config differs
    detector_field = config["detector"]["field"]
    auto_bFieldInZ = detector_field["direction"][2] * detector_field["strength"]
    bFieldInZ = auto_bFieldInZ
    if "bFieldInZ" in options and abs(options["bFieldInZ"] - auto_bFieldInZ) > 1e-6:
        print(
            f"\033[93m[WARNING]\033[0m bFieldInZ in seeding.finderOptions ({options['bFieldInZ']}) differs from detector.field.strength*direction ({auto_bFieldInZ})! Using detector value."
        )
        bFieldInZ = auto_bFieldInZ

    matching_fraction_threshold = config.get("matchingFractionThreshold", 0.5)

    # Print hyperparameters summary
    print("\n========== HYPERPARAMETERS ==========")
    print(
        f"Magnetic field (detector.field.strength): {detector_field['strength']} T"
    )
    print(f"Magnetic field direction: {detector_field['direction']}")
    print(f"bFieldInZ (used): {bFieldInZ} T")
    print(f"Matching fraction threshold: {matching_fraction_threshold}")
    print("\n[Seeding: finderConfig]")
    for k, v in finder.items():
        print(f"  {k}: {v}")
    print("\n[Seeding: filterConfig]")
    filter_cfg = seed.get("filterConfig", {})
    for k, v in filter_cfg.items():
        print(f"  {k}: {v}")
    print("\n[Track finding: selectorConfig]")
    tf_cfg = config.get("trackFinding", {})
    selector_cfg = tf_cfg.get("selectorConfig", {})
    for k, v in selector_cfg.items():
        print(f"  {k}: {v}")
    print("\n[Track finding: ckfConfig]")
    ckf_cfg = tf_cfg.get("ckfConfig", {})
    for k, v in ckf_cfg.items():
        print(f"  {k}: {v}")
    print("====================================\n")

    # Create seeding config
    addSeeding(
        s,
        trackingGeometry,
        field,
        TrackSmearingSigmas(
            loc0=smearing["loc0"],
            loc0PtA=smearing["loc0PtA"],
            loc0PtB=smearing["loc0PtB"],
            loc1=smearing["loc1"],
            loc1PtA=smearing["loc1PtA"],
            loc1PtB=smearing["loc1PtB"],
            time=smearing["time"],
            phi=smearing["phi"],
            theta=smearing["theta"],
            ptRel=smearing["ptRel"],
        ),
        SeedFinderConfigArg(
            r=(
                None if finder["r"]["min"] is None else finder["r"]["min"] * u.mm,
                finder["r"]["max"] * u.mm,
            ),
            deltaR=(finder["deltaR"]["min"] * u.mm, finder["deltaR"]["max"] * u.mm),
            collisionRegion=(
                finder["collisionRegion"]["min"] * u.mm,
                finder["collisionRegion"]["max"] * u.mm,
            ),
            z=(finder["z"]["min"] * u.mm, finder["z"]["max"] * u.mm),
            maxSeedsPerSpM=finder["maxSeedsPerSpM"],
            sigmaScattering=finder["sigmaScattering"],
            radLengthPerSeed=finder["radLengthPerSeed"],
            minPt=finder["minPt"] * u.GeV,
            impactMax=finder["impactMax"] * u.mm,
        ),
        SeedFinderOptionsArg(
            bFieldInZ=bFieldInZ * u.T,
            beamPos=(options["beamPos"][0], options["beamPos"][1]),
        ),
        TruthEstimatedSeedingAlgorithmConfigArg(
            deltaR=(
                truth_est["deltaR"]["min"] * u.mm,
                None
                if truth_est["deltaR"]["max"] is None
                else truth_est["deltaR"]["max"] * u.mm,
            )
        ),
        seedingAlgorithm=(
            SeedingAlgorithm.TruthSmeared
            if seed["truthSmearedSeeded"]
            else (
                SeedingAlgorithm.TruthEstimated
                if seed["truthEstimatedSeeded"]
                else SeedingAlgorithm.GridTriplet
            )
        ),
        initialSigmas=[
            initial["loc0"] * u.mm,
            initial["loc1"] * u.mm,
            initial["phi"] * u.degree,
            initial["theta"] * u.degree,
            initial["qop"] * u.e / u.GeV,
            initial["time"] * u.ns,
        ],
        initialSigmaQoverPt=seed["initialSigmaQoverPt"] * u.e / u.GeV,
        initialSigmaPtRel=seed["initialSigmaPtRel"],
        initialVarInflation=seed["initialVarInflation"],
        geoSelectionConfigFile=geometrySelection,
        outputDirRoot=outputDir,
        rnd=rnd,
    )

    # Track finding configuration and sanity checks
    tf = config["trackFinding"]
    selector = tf["selectorConfig"]
    ckf = tf["ckfConfig"]

    pt_gun_min, pt_gun_max = pg["momentum"]["min"], pg["momentum"]["max"]
    pt_sel_min = pt_sel["min"]
    pt_seed_min = finder["minPt"]
    pt_track_min = selector["pt"]["min"]
    measurements_sel_min = measurements["min"]
    measurements_track_min = selector["nMeasurementsMin"]

    print(
        f"pT [GeV]: Gun({pt_gun_min:.1f}-{pt_gun_max:.1f}), Selection({pt_sel_min:.3f}), Seeding({pt_seed_min:.3f}), Tracking({pt_track_min:.3f})"
    )

    warnings = []
    if pt_gun_min < pt_sel_min:
        warnings.append(
            f"Generating particles below selection threshold ({pt_gun_min:.3f} < {pt_sel_min:.3f} GeV)"
        )
    if pt_gun_min < pt_seed_min:
        warnings.append(
            f"Generating particles below seeding threshold ({pt_gun_min:.3f} < {pt_seed_min:.3f} GeV)"
        )
    if measurements_track_min > measurements_sel_min:
        warnings.append(
            f"Track nMeasurements({measurements_track_min}) > Particle selection({measurements_sel_min})"
        )
    if warnings:
        print(f"\033[93m[WARNING]\033[0m {'; '.join(warnings)}")

    customLogLevel = acts.examples.defaultLogging(
        s, getattr(acts.logging, logging_cfg["level"])
    )

    # Setup track selector configuration
    trackSelectorConfig = TrackSelectorConfig(
        pt=(
            selector["pt"]["min"] * u.GeV,
            None if selector["pt"]["max"] is None else selector["pt"]["max"] * u.GeV,
        ),
        absEta=(
            None if selector["absEta"]["min"] is None else selector["absEta"]["min"],
            None if selector["absEta"]["max"] is None else selector["absEta"]["max"],
        ),
        loc0=(selector["loc0"]["min"] * u.mm, selector["loc0"]["max"] * u.mm),
        nMeasurementsMin=selector["nMeasurementsMin"],
        maxHoles=selector["maxHoles"],
        maxOutliers=selector["maxOutliers"],
    )

    cutSets = [
        acts.TrackSelector.Config(
            **(trackSelectorDefaultKWArgs(trackSelectorConfig))
        )
    ]
    trkSelCfg = cutSets[0]

    # Setup the track finding algorithm
    trackFinder = acts.examples.TrackFindingAlgorithm(
        level=customLogLevel(),
        measurementSelectorCfg=acts.MeasurementSelector.Config(
            [
                (
                    acts.GeometryIdentifier(),
                    (
                        [],
                        [ckf["chi2CutOffMeasurement"]],
                        [ckf["chi2CutOffOutlier"]],
                        [ckf["numMeasurementsCutOff"]],
                    ),
                )
            ]
        ),
        inputMeasurements="measurements",
        inputInitialTrackParameters="estimatedparameters",
        inputSeeds=(
            "estimatedseeds" if ckf["seedDeduplication"] or ckf["stayOnSeed"] else ""
        ),
        outputTracks="ckf_tracks",
        findTracks=acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction(
            trackingGeometry, field, customLogLevel()
        ),
        **acts.examples.defaultKWArgs(
            trackingGeometry=trackingGeometry,
            magneticField=field,
            trackSelectorCfg=trkSelCfg,
            seedDeduplication=ckf["seedDeduplication"],
            stayOnSeed=ckf["stayOnSeed"],
        ),
        # optional loop-protection overrides (only pass if explicitly set)
        **({} if loop_protection is None and loop_fraction is None else {
            **({} if loop_protection is None else {"loopProtection": loop_protection}),
            **({} if loop_fraction is None else {"loopFraction": loop_fraction}),
        }),
    )
    s.addAlgorithm(trackFinder)
    s.addWhiteboardAlias("tracks", trackFinder.config.outputTracks)

    # Add TrackTruthMatcher with configurable matching ratio
    matchAlg = acts.examples.TrackTruthMatcher(
        level=customLogLevel(),
        inputTracks=trackFinder.config.outputTracks,
        inputParticles="particles_selected",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputTrackParticleMatching="ckf_track_particle_matching",
        outputParticleTrackMatching="ckf_particle_track_matching",
        doubleMatching=True,
        matchingRatio=matching_fraction_threshold,
    )
    s.addAlgorithm(matchAlg)
    s.addWhiteboardAlias(
        "track_particle_matching", matchAlg.config.outputTrackParticleMatching
    )
    s.addWhiteboardAlias(
        "particle_track_matching", matchAlg.config.outputParticleTrackMatching
    )

    # Add track writers
    addTrackWriters(
        s,
        name="ckf",
        tracks=trackFinder.config.outputTracks,
        outputDirCsv=outputDir / output["csvSubdir"] if output["outputCsv"] else None,
        outputDirRoot=outputDir,
        writeSummary=True,
        writeStates=output["writeTrackStates"],
        writeFitterPerformance=True,
        writeFinderPerformance=True,
        writeCovMat=False,
        logLevel=getattr(acts.logging, logging_cfg["level"]),
    )


    # writeMatchingDetails=True so we get the `matchingdetails` TTree.
    s.addWriter(
        acts.examples.TrackFinderPerformanceWriter(
            level=getattr(acts.logging, logging_cfg["level"]),
            inputTracks=trackFinder.config.outputTracks,
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            inputParticleTrackMatching="particle_track_matching",
            inputParticleMeasurementsMap="particle_measurements_map",
            filePath=str(outputDir / f"performance_finding_ckf_matchingdetails.root"),
            writeMatchingDetails=True,
        )
    )


    return s


if "__main__" == __name__:
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Perfect Spacepoints Track Finding")
    parser.add_argument(
        "--config", 
        type=str, 
        default="Examples/Configs/perfect-spacepoints-config.json",
        help="Path to configuration file"
    )
    parser.add_argument(
        "--output-dir", 
        type=str, 
        default=".",
        help="Output directory for results"
    )
    parser.add_argument(
        "--input-particles", 
        type=str, 
        default=None,
        help="Path to input particles file "
    )
    # Loop protection is configured in the JSON config under trackFinding.ckfConfig
    parser.add_argument(
        "run_name", 
        nargs='?',
        type=str, 
        default=None,
        help="Name for the output directory (will create a directory with this name)"
    )
    
    args = parser.parse_args()
    
    # Load configuration
    config_file = Path(args.config)
    if not config_file.is_absolute():
        srcdir = Path(__file__).resolve().parent.parent.parent.parent
        config_file = srcdir / config_file
    
    config = load_config(config_file)
    
    # Setup detector
    detector = config['detector']
    field_config = detector['field']
    detector_obj = GenericDetector()
    trackingGeometry = detector_obj.trackingGeometry()
    decorators = detector_obj.contextDecorators()
    field = acts.ConstantBField(acts.Vector3(
        field_config['direction'][0] * field_config['strength'] * u.T,
        field_config['direction'][1] * field_config['strength'] * u.T,
        field_config['direction'][2] * field_config['strength'] * u.T
    ))
    
    # Setup paths
    srcdir = Path(__file__).resolve().parent.parent.parent.parent
    geometrySelection = srcdir / "Examples/Configs" / detector['geometry']['seedingConfig']
    
    #digiConfigFile = srcdir / "Examples/Configs/generic-digi-perfect-config.json"
    digiConfigFile = srcdir / "Examples/Configs/generic-digi-smearing-config.json"
    
    inputParticlePath = None
    if args.input_particles:
        inputParticlePath = Path(args.input_particles)
        if not inputParticlePath.exists():
            print(f"Warning: Input particles file {inputParticlePath} not found")
            inputParticlePath = None
    
    # Setup output directory - prompt for run_name if not provided
    if args.run_name is None:
        run_name = input("Enter a name for this run (output directory): ").strip()
        if not run_name:
            run_name = "output"
        output_path = Path(run_name)
    else:
        output_path = Path(args.run_name)
    
    # Create output directory if it doesn't exist
    output_path.mkdir(parents=True, exist_ok=True)
    print(f"Output will be saved to: {output_path.resolve()}")
    
    # Read loop-protection settings from the JSON config (if present)
    tf_cfg = config.get("trackFinding", {})
    ckf_cfg = tf_cfg.get("ckfConfig", {})
    loop_protection_cfg = ckf_cfg.get("loopProtection", None)
    loop_fraction_cfg = ckf_cfg.get("loopFraction", None)

    # Run simulation
    sequencer = runPerfectSpacepoints(
        config=config,
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        geometrySelection=geometrySelection,
        digiConfigFile=digiConfigFile,
        field=field,
        outputDir=output_path,
        inputParticlePath=inputParticlePath,
        loop_protection=loop_protection_cfg,
        loop_fraction=loop_fraction_cfg,
    )
    sequencer.run()

