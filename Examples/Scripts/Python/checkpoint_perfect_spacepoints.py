#!/usr/bin/env python3

from pathlib import Path
from typing import Optional, Dict, Any
import argparse
import json
import numpy as np
import matplotlib.pyplot as plt

import acts
from acts import UnitConstants as u
from acts.examples import GenericDetector, RootParticleReader

# Mapping from config string to PDG enum (lowercase keys)
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
):
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
        addCKFTracks,
        TrackSelectorConfig,
        CkfConfig,
    )

    # Get configuration directly from JSON
    sim = config['simulation']
    pg = sim['particleGun']
    output = config['output']
    logging = config['logging']
    
    s = s or acts.examples.Sequencer(
        events=sim['events'], 
        numThreads=sim['numThreads'], 
        logLevel=getattr(acts.logging, logging['level'])
    )
    for d in decorators:
        s.addContextDecorator(d)
    rnd = acts.examples.RandomNumbers(seed=sim['randomSeed'])
    outputDir = Path(outputDir)

    if inputParticlePath is None:
        # Create particle gun config from JSON
        momentum = pg['momentum']
        eta = pg['eta']
        phi = pg['phi']
        particles = pg['particles']
        ptype_key = particles['type'].lower().replace(" ", "")
        if ptype_key not in PDG_MAP:
            raise ValueError(f"Unknown particle type '{particles['type']}'. Allowed: {', '.join(PDG_MAP.keys())}")
        pdg = PDG_MAP[ptype_key]
        addParticleGun(
            s,
            MomentumConfig(
                momentum['min'] * u.GeV,
                momentum['max'] * u.GeV,
                transverse=momentum['transverse']
            ),
            EtaConfig(eta['min'], eta['max'], uniform=eta['uniform']),
            PhiConfig(phi['min'] * u.degree, phi['max'] * u.degree),
            ParticleConfig(
                particles['count'],
                pdg,
                randomizeCharge=particles['randomizeCharge']
            ),
            multiplicity=pg['multiplicity'],
            rnd=rnd,
        )
    else:
        acts.logging.getLogger("PerfectSpacepointsExample").info(
            "Reading particles from %s", inputParticlePath.resolve()
        )
        assert inputParticlePath.exists()
        s.addReader(
            RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                outputParticles="particles_generated",
            )
        )

    # Get physics and selection config directly from JSON
    physics = config['simulation_physics']
    selection = config['particleSelection']
    
    # Perfect simulation - no material interactions
    addFatras(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        enableInteractions=physics['enableInteractions'],
    )
 

    # Perfect digitization - no smearing, no noise
    # Note: doSmearing, doNoise, doInefficiency are controlled by digiConfigFile
    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
    )

    # Create particle selection config from JSON
    pt = selection['pt']
    measurements = selection['measurements']
    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(pt['min'] * u.GeV, None if pt['max'] is None else pt['max'] * u.GeV),
            measurements=(measurements['min'], None if measurements['max'] is None else measurements['max']),
            removeNeutral=selection['removeNeutral']
        )
    )

    # Get seeding configuration directly from JSON
    seed = config['seeding']
    smearing = seed['smearing']
    finder = seed['finderConfig']
    options = seed['finderOptions']
    truth_est = seed['truthEstimated']
    initial = seed['initialSigmas']
    
    # Create seeding config from JSON
    addSeeding(
        s,
        trackingGeometry,
        field,
        TrackSmearingSigmas(
            loc0=smearing['loc0'],
            loc0PtA=smearing['loc0PtA'],
            loc0PtB=smearing['loc0PtB'],
            loc1=smearing['loc1'],
            loc1PtA=smearing['loc1PtA'],
            loc1PtB=smearing['loc1PtB'],
            time=smearing['time'],
            phi=smearing['phi'],
            theta=smearing['theta'],
            ptRel=smearing['ptRel']
        ),
        SeedFinderConfigArg(
            r=(None if finder['r']['min'] is None else finder['r']['min'] * u.mm,
               finder['r']['max'] * u.mm),
            deltaR=(finder['deltaR']['min'] * u.mm, finder['deltaR']['max'] * u.mm),
            collisionRegion=(finder['collisionRegion']['min'] * u.mm, 
                           finder['collisionRegion']['max'] * u.mm),
            z=(finder['z']['min'] * u.mm, finder['z']['max'] * u.mm),
            maxSeedsPerSpM=finder['maxSeedsPerSpM'],
            sigmaScattering=finder['sigmaScattering'],
            radLengthPerSeed=finder['radLengthPerSeed'],
            minPt=finder['minPt'] * u.MeV,
            impactMax=finder['impactMax'] * u.mm
        ),
        SeedFinderOptionsArg(
            bFieldInZ=options['bFieldInZ'] * u.T,
            beamPos=(options['beamPos'][0], options['beamPos'][1])
        ),
        TruthEstimatedSeedingAlgorithmConfigArg(
            deltaR=(truth_est['deltaR']['min'] * u.mm, 
                   None if truth_est['deltaR']['max'] is None else truth_est['deltaR']['max'] * u.mm)
        ),
        seedingAlgorithm=(
            SeedingAlgorithm.TruthSmeared
            if seed['truthSmearedSeeded']
            else (
                SeedingAlgorithm.TruthEstimated
                if seed['truthEstimatedSeeded']
                else SeedingAlgorithm.GridTriplet
            )
        ),
        initialSigmas=[
            initial['loc0'] * u.mm,
            initial['loc1'] * u.mm,
            initial['phi'] * u.degree,
            initial['theta'] * u.degree,
            initial['qop'] * u.e / u.GeV,
            initial['time'] * u.ns
        ],
        initialSigmaQoverPt=seed['initialSigmaQoverPt'] * u.e / u.GeV,
        initialSigmaPtRel=seed['initialSigmaPtRel'],
        initialVarInflation=seed['initialVarInflation'],
        geoSelectionConfigFile=geometrySelection,
        outputDirRoot=outputDir,
        rnd=rnd,
    )

    # Get track finding configuration from JSON
    tf = config['trackFinding']
    selector = tf['selectorConfig']
    ckf = tf['ckfConfig']
    
    # Track finding with perfect conditions
    addCKFTracks(
        s,
        trackingGeometry,
        field,
        TrackSelectorConfig(
            pt=(selector['pt']['min'] * u.MeV, 
               None if selector['pt']['max'] is None else selector['pt']['max'] * u.MeV),
            absEta=(None if selector['absEta']['min'] is None else selector['absEta']['min'],
                   None if selector['absEta']['max'] is None else selector['absEta']['max']),
            loc0=(selector['loc0']['min'] * u.mm, selector['loc0']['max'] * u.mm),
            nMeasurementsMin=selector['nMeasurementsMin'],
            maxHoles=selector['maxHoles'],
            maxOutliers=selector['maxOutliers']
        ),
        CkfConfig(
            chi2CutOffMeasurement=ckf['chi2CutOffMeasurement'],
            chi2CutOffOutlier=ckf['chi2CutOffOutlier'],
            numMeasurementsCutOff=ckf['numMeasurementsCutOff'],
            seedDeduplication=ckf['seedDeduplication'],
            stayOnSeed=ckf['stayOnSeed']
        ),
        outputDirRoot=outputDir,
        outputDirCsv=outputDir / output['csvSubdir'] if output['outputCsv'] else None,
        writeTrackStates=output['writeTrackStates'],
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
        help="Path to input particles file (optional)"
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
    
    # Run simulation
    sequencer = runPerfectSpacepoints(
        config=config,
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        geometrySelection=geometrySelection,
        digiConfigFile=digiConfigFile,
        field=field,
        outputDir=Path(args.output_dir),
        inputParticlePath=inputParticlePath,
    )
    sequencer.run()


