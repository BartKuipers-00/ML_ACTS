#!/usr/bin/env python3
"""
Configuration is loaded from Examples/Configs/perfect-spacepoints-multigen-config.json
"""

from pathlib import Path
from typing import Optional, Dict, Any
import argparse
import json

import acts
import acts.examples
import acts.examples.hepmc3
from acts import UnitConstants as u
from acts.examples import GenericDetector, RootParticleWriter

from acts.examples.simulation import (
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
    "gamma": acts.PdgParticle.eGamma,
}


def load_config(config_file: Path) -> Dict[str, Any]:
    """Load configuration from JSON file."""
    with open(config_file, "r") as f:
        return json.load(f)


def runPerfectSpacepointsMultiGen(
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
    sim = config["simulation"]
    eg_cfg = sim["eventGenerator"]
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

    customLogLevel = acts.examples.defaultLogging(
        s, getattr(acts.logging, logging_cfg["level"])
    )

    # -------------------------------------------------------------------------
    # Particle source


# Build one EventGenerator.Generator per requested particle species
    particle_types = eg_cfg["particleTypes"]
    if isinstance(particle_types, str):
        particle_types = [particle_types]

    pdg_list = []
    for pt in particle_types:
        key = pt.lower().replace(" ", "")
        if key not in PDG_MAP:
            raise ValueError(
                f"Unknown particle type '{key}'. Allowed: {', '.join(PDG_MAP.keys())}"
            )
        pdg_list.append(PDG_MAP[key])

    momentum = eg_cfg["momentum"]
    eta_cfg = eg_cfg["eta"]
    phi_cfg = eg_cfg["phi"]

    vtxGen = acts.examples.GaussianVertexGenerator(
        mean=acts.Vector4(0, 0, 0, 0),
        stddev=acts.Vector4(0, 0, 0, 0),
    )

    evGen = acts.examples.EventGenerator(
        level=customLogLevel(),
        generators=[
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(
                    n=eg_cfg["multiplicity"]
                ),
                vertex=vtxGen,
                particles=acts.examples.ParametricParticleGenerator(
                    **acts.examples.defaultKWArgs(
                        p=(momentum["min"], momentum["max"]),
                        pTransverse=momentum["transverse"],
                        eta=(eta_cfg["min"], eta_cfg["max"]),
                        phi=(phi_cfg["min"] * u.degree, phi_cfg["max"] * u.degree),
                        etaUniform=eta_cfg["uniform"],
                        numParticles=eg_cfg["particlesPerVertex"],
                        pdg=pdg,
                        randomizeCharge=eg_cfg["randomizeCharge"],
                    )
                ),
            )
            for pdg in pdg_list
        ],
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
    s.addWhiteboardAlias(
        "particles_generated_selected", hepmc3Converter.config.outputParticles
    )

    s.addWriter(
        RootParticleWriter(
            level=customLogLevel(),
            inputParticles=hepmc3Converter.config.outputParticles,
            filePath=str(outputDir / "particles.root"),
        )
    )

    # -------------------------------------------------------------------------
    # Physics simulation and digitization
    # -------------------------------------------------------------------------
    physics = config["simulation_physics"]
    selection = config["particleSelection"]
    pMin = physics["pMin"]
    if pMin is not None:
        pMin = pMin * u.GeV

    nav_cfg = config.get("navigation", {})
    max_surface_retargets = int(nav_cfg.get("max_surface_retargets", 3))

    addFatras(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        enableInteractions=physics["enableInteractions"],
        pMin=pMin,
        maxSteps=sim.get("maxSteps", None),
        loopFraction=sim.get("loopFraction", None),
        debugStepInterval=logging_cfg.get("n", None),
        maxSurfaceRetargets=max_surface_retargets,
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
            pt=(
                pt_sel["min"] * u.GeV,
                None if pt_sel["max"] is None else pt_sel["max"] * u.GeV,
            ),
            measurements=(
                measurements["min"],
                None if measurements["max"] is None else measurements["max"],
            ),
            removeNeutral=selection["removeNeutral"],
        ),
    )

    # -------------------------------------------------------------------------
    # Seeding
    # -------------------------------------------------------------------------
    seed = config["seeding"]
    smearing = seed["smearing"]
    finder = seed["finderConfig"]
    options = seed["finderOptions"]
    truth_est = seed["truthEstimated"]
    initial = seed["initialSigmas"]

    detector_field = config["detector"]["field"]
    auto_bFieldInZ = detector_field["direction"][2] * detector_field["strength"]
    bFieldInZ = auto_bFieldInZ
    if "bFieldInZ" in options and abs(options["bFieldInZ"] - auto_bFieldInZ) > 1e-6:
           print(
            f"[WARNING] bFieldInZ in seeding.finderOptions "
            f"({options['bFieldInZ']}) differs from detector.field "
            f"({auto_bFieldInZ})! Using detector value."
        )

    matching_fraction_threshold = config.get("matchingFractionThreshold", 0.5)

    # Hyperparameter summary
    print("\n========== HYPERPARAMETERS ==========")
    print(f"Magnetic field strength: {detector_field['strength']} T")
    print(f"Magnetic field direction: {detector_field['direction']}")
    print(f"bFieldInZ (used): {bFieldInZ} T")
    print(f"Matching fraction threshold: {matching_fraction_threshold}")
    if inputParticlePath is None:
        print(f"Particle species: {eg_cfg['particleTypes']}")
        print(
            f"Particles per vertex per species: {eg_cfg['particlesPerVertex']}  "
            f"×  multiplicity {eg_cfg['multiplicity']}"
        )
        print(
            f"pT: {eg_cfg['momentum']['min']}-{eg_cfg['momentum']['max']} GeV"
        )
    print("\n[Seeding: finderConfig]")
    for k, v in finder.items():
        print(f"  {k}: {v}")
    filter_cfg = seed.get("filterConfig", {})
    print("\n[Seeding: filterConfig]")
    for k, v in filter_cfg.items():
        print(f"  {k}: {v}")
    tf_cfg = config.get("trackFinding", {})
    print("\n[Track finding: selectorConfig]")
    for k, v in tf_cfg.get("selectorConfig", {}).items():
        print(f"  {k}: {v}")
    print("\n[Track finding: ckfConfig]")
    for k, v in tf_cfg.get("ckfConfig", {}).items():
        print(f"  {k}: {v}")
    print("====================================\n")

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
        # Skip the framework auto-writer; we add our own writers below with
        # fine pT binning suitable for low-pT seeding-efficiency analysis.
        outputDirRoot=None,
        rnd=rnd,
    )

    # -------------------------------------------------------------------------
    # Seeding performance output (fine pT binning, low-pT analysis)
    # -------------------------------------------------------------------------
    if not seed["truthSmearedSeeded"]:
        _seedingEffBinning = {
            "Pt":     acts.examples.Binning("pT [GeV/c]", 40, 0.0, 0.5),
            "Eta":    acts.examples.Binning("#eta", 40, -4, 4),
            "Phi":    acts.examples.Binning("#phi", 100, -3.15, 3.15),
            "Z0":     acts.examples.Binning("z_0 [mm]", 50, -200, 200),
            "DeltaR": acts.examples.Binning("#Delta R", 100, 0, 0.3),
            "prodR":  acts.examples.Binning("prod_R [mm]", 100, 0, 200),
        }
        s.addWriter(
            acts.examples.TrackFinderPerformanceWriter(
                level=getattr(acts.logging, logging_cfg["level"]),
                inputTracks="seed-tracks",
                inputParticles="particles_selected",
                inputTrackParticleMatching="seed_particle_matching",
                inputParticleTrackMatching="particle_seed_matching",
                inputParticleMeasurementsMap="particle_measurements_map",
                filePath=str(outputDir / "performance_seeding.root"),
                effPlotToolConfig=acts.examples.EffPlotToolConfig(_seedingEffBinning),
                writeMatchingDetails=True,
            )
        )
        s.addWriter(
            acts.examples.RootTrackParameterWriter(
                level=getattr(acts.logging, logging_cfg["level"]),
                inputTrackParameters="estimatedparameters",
                inputProtoTracks="seed-prototracks",
                inputParticles="particles",
                inputSimHits="simhits",
                inputMeasurementParticlesMap="measurement_particles_map",
                inputMeasurementSimHitsMap="measurement_simhits_map",
                filePath=str(outputDir / "estimatedparams.root"),
                treeName="estimatedparams",
            )
        )

    # -------------------------------------------------------------------------
    # Track finding (CKF)
    # -------------------------------------------------------------------------
    tf = config["trackFinding"]
    selector = tf["selectorConfig"]
    ckf = tf["ckfConfig"]

    pt_seed_min = finder["minPt"]
    pt_track_min = selector["pt"]["min"]
    pt_sel_min = pt_sel["min"]
    measurements_sel_min = measurements["min"]
    measurements_track_min = selector["nMeasurementsMin"]

    warnings = []
    if not inputParticlePath:
        pt_gun_min = eg_cfg["momentum"]["min"]
        if pt_gun_min < pt_sel_min:
            warnings.append(
                f"Generating particles below selection threshold "
                f"({pt_gun_min:.3f} < {pt_sel_min:.3f} GeV)"
            )
        if pt_gun_min < pt_seed_min:
            warnings.append(
                f"Generating particles below seeding threshold "
                f"({pt_gun_min:.3f} < {pt_seed_min:.3f} GeV)"
            )
    if measurements_track_min > measurements_sel_min:
        warnings.append(
            f"Track nMeasurements({measurements_track_min}) > "
            f"Particle selection({measurements_sel_min})"
        )
    if warnings:
        print("[WARNING]", "; ".join(warnings))

    customLogLevel = acts.examples.defaultLogging(
        s, getattr(acts.logging, logging_cfg["level"])
    )

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
        acts.TrackSelector.Config(**(trackSelectorDefaultKWArgs(trackSelectorConfig)))
    ]
    trkSelCfg = cutSets[0]

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
            trackingGeometry,
            field,
            customLogLevel(),
            maxSurfaceRetargets=max_surface_retargets,
        ),
        **acts.examples.defaultKWArgs(
            trackingGeometry=trackingGeometry,
            magneticField=field,
            trackSelectorCfg=trkSelCfg,
            seedDeduplication=ckf["seedDeduplication"],
            stayOnSeed=ckf["stayOnSeed"],
            loopProtection=loop_protection,
            loopFraction=loop_fraction,
        ),
    )
    s.addAlgorithm(trackFinder)
    s.addWhiteboardAlias("tracks", trackFinder.config.outputTracks)

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

    s.addWriter(
        acts.examples.TrackFinderPerformanceWriter(
            level=getattr(acts.logging, logging_cfg["level"]),
            inputTracks=trackFinder.config.outputTracks,
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            inputParticleTrackMatching="particle_track_matching",
            inputParticleMeasurementsMap="particle_measurements_map",
            filePath=str(
                outputDir / "performance_finding_ckf_matchingdetails.root"
            ),
            writeMatchingDetails=True,
        )
    )

    return s


if "__main__" == __name__:
    parser = argparse.ArgumentParser(
        description="Perfect Spacepoints Track Finding — multi-species EventGenerator"
    )
    parser.add_argument(
        "output_dir",
        nargs="?",
        type=str,
        default=None,
        help="Output directory for results",
    )
    args = parser.parse_args()

    output_dir = args.output_dir if args.output_dir is not None else "."

    srcdir = Path(__file__).resolve().parent.parent.parent.parent
    config_file = srcdir / "Examples/Configs/perfect-spacepoints-multigen-config.json"
    config = load_config(config_file)

    detector = config["detector"]
    field_config = detector["field"]
    detector_obj = GenericDetector()
    trackingGeometry = detector_obj.trackingGeometry()
    decorators = detector_obj.contextDecorators()
    field = acts.ConstantBField(
        acts.Vector3(
            field_config["direction"][0] * field_config["strength"] * u.T,
            field_config["direction"][1] * field_config["strength"] * u.T,
            field_config["direction"][2] * field_config["strength"] * u.T,
        )
    )

    geometrySelection = (
        srcdir / "Examples/Configs" / detector["geometry"]["seedingConfig"]
    )
    digiConfigFile = srcdir / "Examples/Configs/generic-digi-smearing-config.json"

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    print(f"Output will be saved to: {output_path.resolve()}")

    sim_config = config.get("simulation", {})
    sequencer = runPerfectSpacepointsMultiGen(
        config=config,
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        geometrySelection=geometrySelection,
        digiConfigFile=digiConfigFile,
        field=field,
        outputDir=output_path,
        inputParticlePath=None,
        loop_fraction=sim_config.get("loopFraction", None),
    )
    # Snapshot the per-stage step-limit drop counters: how many particles
    # (Fatras) / seed tracks (CKF) were terminated because the propagator
    # hit PropagatorError::StepCountLimitReached. try/finally so the
    # diagnostic still prints if the Sequencer aborts (e.g. unmasked FPEs).
    acts.resetStepLimitCounts()
    acts.resetApexInsideShellCounts()
    n_events_run = sim_config.get("events", "?")
    try:
        sequencer.run()
    finally:
        sim_failed   = acts.fatrasFailedCount()
        sim_steplim  = acts.fatrasStepLimitCount()
        reco_failed  = acts.trackFindingFailedCount()
        reco_steplim = acts.trackFindingStepLimitCount()
        print(f"[failures diagnostic] events={n_events_run}  "
              f"fatras_failed={sim_failed} (of which step_limit={sim_steplim})  "
              f"ckf_failed={reco_failed} (of which step_limit={reco_steplim})  "
              f"total={sim_failed + reco_failed}")
        # Navigator-internal sensitive-surface retarget counters.
        # invoked   = how many times Layer::compatibleSurfaces was re-run after
        #             a sensitive module returned unreachable mid-step.
        # succeeded = subset where the propagator subsequently reached at least
        #             one sensitive module on the same layer entry (counted
        #             at most once per layer per track).
        fatras_rt_inv = acts.fatrasRetargetInvokedCount()
        fatras_rt_ok  = acts.fatrasRetargetSucceededCount()
        ckf_rt_inv    = acts.ckfRetargetInvokedCount()
        ckf_rt_ok     = acts.ckfRetargetSucceededCount()
        def _pct(n, d):
            return 0.0 if d == 0 else 100.0 * n / d
        print(f"[retarget diagnostic] "
              f"fatras_retarget_invoked={fatras_rt_inv} succeeded={fatras_rt_ok} "
              f"({_pct(fatras_rt_ok, fatras_rt_inv):.1f}%)  "
              f"ckf_retarget_invoked={ckf_rt_inv} succeeded={ckf_rt_ok} "
              f"({_pct(ckf_rt_ok, ckf_rt_inv):.1f}%)")
        # Apex-inside-shell diagnostic (temporary instrumentation). Per
        # barrel sensor layer: how many times a turning point fired
        # geometrically inside that layer's apr=1/apr=2 shell, and how
        # many of those fire-events were followed by a sensor hit on the
        # same layer (recovery). not_recovered = inside - recovered.
        snap = acts.apexInsideShellSnapshot()
        if snap:
            print(f"[apex-inside-shell diagnostic]")
            print(f"  {'layer':<10} {'inside':>8} {'recovered':>10} "
                  f"{'not_recovered':>14} {'recovery%':>10}")
            t_in = t_rec = 0
            for vol, lay, n_inside, n_rec in sorted(snap):
                t_in += n_inside
                t_rec += n_rec
                print(f"  V{vol}L{lay:<7} {n_inside:>8} {n_rec:>10} "
                      f"{n_inside - n_rec:>14} "
                      f"{_pct(n_rec, n_inside):>9.1f}%")
            print(f"  {'TOTAL':<10} {t_in:>8} {t_rec:>10} "
                  f"{t_in - t_rec:>14} {_pct(t_rec, t_in):>9.1f}%")
