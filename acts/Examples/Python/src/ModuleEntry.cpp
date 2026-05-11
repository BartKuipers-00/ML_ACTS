// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/ActsVersion.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Propagator/StepLimitDiagnostics.hpp"

#include <tuple>
#include <unordered_map>

#include <pybind11/detail/common.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pyerrors.h>

namespace py = pybind11;
using namespace Acts::Python;

namespace Acts::Python {
void addContext(Context& ctx);
void addAny(Context& ctx);
void addUnits(Context& ctx);
void addFramework(Context& ctx);
void addLogging(Context& ctx);
void addPdgParticle(Context& ctx);
void addAlgebra(Context& ctx);
void addBinning(Context& ctx);
void addEventData(Context& ctx);

void addPropagation(Context& ctx);
void addNavigation(Context& ctx);

void addAlignment(Context& ctx);
void addGeometry(Context& ctx);
void addGeometryBuildingGen1(Context& ctx);
void addExperimentalGeometry(Context& ctx);

void addMagneticField(Context& ctx);

void addMaterial(Context& ctx);
void addOutput(Context& ctx);
void addDetector(Context& ctx);
void addExampleAlgorithms(Context& ctx);
void addInput(Context& ctx);
void addGenerators(Context& ctx);
void addTruthTracking(Context& ctx);
void addTrackFitting(Context& ctx);
void addTrackFinding(Context& ctx);
void addTruthJet(Context& ctx);
void addVertexing(Context& ctx);
void addAmbiguityResolution(Context& ctx);
void addUtilities(Context& ctx);

void addRootInput(Context& ctx);
void addRootOutput(Context& ctx);

// Plugins
void addDigitization(Context& ctx);
void addPythia8(Context& ctx);
void addGeoModel(Context& ctx);
void addTGeo(Context& ctx);
void addJson(Context& ctx);
void addDetray(Context& ctx);
void addHepMC3(Context& ctx);
void addGnnTrackFinding(Context& ctx);
void addSvg(Context& ctx);
void addObj(Context& ctx);
void addOnnx(Context& ctx);
void addOnnxNeuralCalibrator(Context& ctx);
void addCovfie(Context& ctx);
void addTraccc(Context& ctx);
void addHashing(Context& ctx);

}  // namespace Acts::Python

PYBIND11_MODULE(ActsPythonBindings, m) {
  Acts::Python::Context ctx;
  ctx.modules["main"] = m;
  auto mex = m.def_submodule("_examples");
  ctx.modules["examples"] = mex;
  auto prop = m.def_submodule("_propagator");
  ctx.modules["propagation"] = prop;
  m.doc() = "Acts";

  m.attr("__version__") =
      std::tuple{Acts::VersionMajor, Acts::VersionMinor, Acts::VersionPatch};

  {
    auto mv = m.def_submodule("version");

    mv.attr("major") = Acts::VersionMajor;
    mv.attr("minor") = Acts::VersionMinor;
    mv.attr("patch") = Acts::VersionPatch;

    mv.attr("commit_hash") = Acts::CommitHash;
    mv.attr("commit_hash_short") = Acts::CommitHashShort;
  }

  // Diagnostic counters split by stage (Fatras vs CKF) and by failure type:
  //   *Failed     — increments on ANY non-ok propagator result
  //   *StepLimit  — subset: only PropagatorError::StepCountLimitReached
  m.def("fatrasFailedCount", &Acts::detail::fatrasFailedCount);
  m.def("fatrasStepLimitCount", &Acts::detail::fatrasStepLimitCount);
  m.def("trackFindingFailedCount", &Acts::detail::trackFindingFailedCount);
  m.def("trackFindingStepLimitCount",
        &Acts::detail::trackFindingStepLimitCount);
  // Sensitive-surface retarget counters (Navigator-internal). "Invoked" =
  // every time the navigator re-runs Layer::compatibleSurfaces because a
  // sensitive surface returned unreachable mid-step. "Succeeded" = subset
  // where the propagator subsequently received onSurface for a sensitive
  // module on the same layer entry (counted at most once per layer per track).
  m.def("fatrasRetargetInvokedCount",
        &Acts::detail::fatrasRetargetInvokedCount);
  m.def("fatrasRetargetSucceededCount",
        &Acts::detail::fatrasRetargetSucceededCount);
  m.def("ckfRetargetInvokedCount", &Acts::detail::ckfRetargetInvokedCount);
  m.def("ckfRetargetSucceededCount",
        &Acts::detail::ckfRetargetSucceededCount);
  // Per-layer apex-inside-shell diagnostic. Returns a list of (vol, lay,
  // n_inside, n_recovered) tuples — one per barrel sensor layer that fired
  // a turning-point inside its apr=1/apr=2 shell. n_inside - n_recovered
  // is "TP fired inside this layer's shell, but no incoming-arc sensor hit
  // followed on the same layer".
  m.def("apexInsideShellSnapshot", &Acts::detail::apexInsideShellSnapshot);
  m.def("resetApexInsideShellCounts",
        &Acts::detail::resetApexInsideShellCounts);
  m.def("resetStepLimitCounts", &Acts::detail::resetStepLimitCounts);

  addContext(ctx);
  addAny(ctx);
  addUnits(ctx);
  addFramework(ctx);
  addLogging(ctx);
  addPdgParticle(ctx);
  addAlgebra(ctx);
  addBinning(ctx);
  addEventData(ctx);
  addOutput(ctx);

  addPropagation(ctx);
  addNavigation(ctx);
  addAlignment(ctx);
  addGeometryBuildingGen1(ctx);
  addGeometry(ctx);
  addExperimentalGeometry(ctx);

  addMagneticField(ctx);
  addMaterial(ctx);
  addDetector(ctx);
  addExampleAlgorithms(ctx);
  addInput(ctx);
  addGenerators(ctx);
  addTruthTracking(ctx);
  addTrackFitting(ctx);
  addTrackFinding(ctx);
  addTruthJet(ctx);
  addVertexing(ctx);
  addAmbiguityResolution(ctx);
  addUtilities(ctx);

  addDigitization(ctx);
  addPythia8(ctx);
  addJson(ctx);
  addGeoModel(ctx);
  addTGeo(ctx);
  addDetray(ctx);
  addHepMC3(ctx);
  addGnnTrackFinding(ctx);
  addObj(ctx);
  addSvg(ctx);
  addOnnx(ctx);
  addOnnxNeuralCalibrator(ctx);
  addCovfie(ctx);
  addTraccc(ctx);
  addHashing(ctx);

  addRootInput(ctx);
  addRootOutput(ctx);
}
