// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/HelixIntersection.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>

namespace Acts::detail {

/// Update surface status - Single component
///
/// This method intersect the provided surface and update the navigation
/// step estimation accordingly (hence it changes the state). It also
/// returns the status of the intersection to trigger onSurface in case
/// the surface is reached.
///
/// @tparam stepper_t The type of stepper used for the propagation
///
/// @param stepper [in] The stepper in use
/// @param state [in,out] The stepping state (thread-local cache)
/// @param surface [in] The surface provided
/// @param index [in] The surface intersection index
/// @param direction [in] The propagation direction
/// @param boundaryTolerance [in] The boundary check for this status update
/// @param surfaceTolerance [in] Surface tolerance used for intersection
/// @param stype [in] The step size type to be set
/// @param isInBarrelVolume [in] Flag indicating if in barrel volume (for radial direction logic)
/// @param logger [in] A @c Logger instance
template <typename stepper_t>
IntersectionStatus updateSingleSurfaceStatus(
    const stepper_t& stepper, typename stepper_t::State& state,
    const Surface& surface, std::uint8_t index, Direction direction,
    const BoundaryTolerance& boundaryTolerance, double surfaceTolerance,
    ConstrainedStep::Type stype, bool isInBarrelVolume, bool shellMode,
    const Logger& logger) {
  ACTS_VERBOSE("Update single surface status for surface: "
               << surface.geometryId() << " index " << static_cast<int>(index));

  // Check for radial momentum flip (turning point detection for spiraling particles)
  Vector3 position = stepper.position(state);
  Vector3 dir = stepper.direction(state);

  // projection on radial direction p_r = p(cartesian) dot r_hat = (p_x * x + p_y * y) / sqrt(x^2 + y^2)
  double r_xy = std::sqrt(position[0] * position[0] + position[1] * position[1]);
  bool radiallyInward_current = false;
  if (r_xy > 1e-6) {  // Avoid division by zero
    double pr = (position[0] * dir[0] + position[1] * dir[1]) / r_xy;
    radiallyInward_current = (pr <= 0.0);  // Inward or tangential (use radial direction)
  }

  // Detect change in radial direction (turning point)
  if (radiallyInward_current != state.radiallyInward_previous) {
    ACTS_VERBOSE("Turning point detected: radiallyInward changed from "
                 << state.radiallyInward_previous << " to "
                 << radiallyInward_current);
    state.radiallyInward_previous = radiallyInward_current;
    state.turningPointDetected = true;  // Set flag for propagator & navigator
    return IntersectionStatus::unreachable;
  }

  // Update for next step
  state.radiallyInward_previous = radiallyInward_current;

  // Predicate / direction selection per regime:
  //   (1) shell mode (gated on by the navigator when the trajectory's
  //       apex sits inside the current layer's shell): the block below
  //       compute the helix arc length to r_ideal and either skip module
  //       probing (FAR branch) or use the BISECTOR direction
  //       normalize(tangent + (-r̂)) for the line-plane intersect
  //       (CLOSE branch).
  //   (2) inward-barrel non-sensitive (approach / boundary / layer-rep)
  //       outside shell mode → pure −r̂ direction with the standard
  //       line predicate. These cylinders are radially well-conditioned
  //       in the barrel (≈⊥ r̂), giving the cleanest line-plane intersect.
  //   (3) inward-barrel sensitive outside shell mode → closed-form
  //       helix-plane intersect using the raw tangent, computed in the
  //       lambda below. Matches the actual curved landing on the module
  //       rectangle.
  //   (4) anything else (outward, non-barrel) → raw tangent + line.
  Vector3 intersectionDirection = direction * stepper.direction(state);
  // NOTE: do NOT use the `isInBarrelVolume` parameter here. It comes from
  // state.navigation.isInBarrelVolume which is only set inside
  // computeEffectiveDirection (called from compatibleLayers). After a TP-fix
  // surfaceTarget restart, compatibleLayers is skipped and the flag stays
  // stale — so a sensor would be mis-classified as non-sensitive, and the
  // CLOSE-branch shell predicate switches to pure-radial-inward,
  // landing at slightly negative pathLength → trajectory steps backward →
  // pr flips → TP cascade. The sensitive flag + r_xy guard alone are
  // sufficient: when shellMode is active, the navigator has already
  // verified barrel + sensor-layer + in-shell; outside shell, the
  // helix-plane predicate has its own internal precondition fallback for
  // non-barrel cases.
  const bool barrelSensitive =
      r_xy > 1e-6 && surface.geometryId().sensitive() != 0;

  // Shell mode: per-step gated module probing.
  //   1. Compute d = helix arc length from (pos, dir) to r_ideal (the
  //      layer's representing-surface cylinder).
  //   2. FAR branch (d > probeTrigger AND r_xy > r_ideal + 0.5): set
  //      navigator step slot to (d - safetyOffset). The trajectory
  //      descends along its curve but stops `safetyOffset` mm of arc
  //      short of r_ideal, so it never lands EXACTLY on the sensor
  //      cylinder where the bisector predicate's path collapses to ~0.
  //   3. CLOSE branch (d <= probeTrigger or trajectory past r_ideal):
  //      swap intersect direction to the BISECTOR normalize(tangent +
  //      (-r̂)) and fall through to the standard surface.intersect path
  //      below. From this point we drop the helix-arc approximation
  //      entirely — bisector predicate determines path length AND the
  //      onSurface check fires naturally from surface.intersect.
  //   4. probeTrigger (3.2) is slightly LARGER than safetyOffset (3.0)
  //      so the CLOSE switch happens reliably without needing FP
  //      precision on `d == safetyOffset` exactly.
  //   5. Safety bound on r_xy: helixBarrelCylinderLanding returns the
  //      smaller FORWARD arc. Once r_xy < r_ideal, the next forward
  //      crossing is the OUTGOING one after the apex turn (~3000 mm).
  //      Force CLOSE in that regime by gating FAR on
  //      `r_xy > r_ideal + 0.5`.
  constexpr double s_shellSafetyOffset = 5.0;
  constexpr double s_shellProbeTrigger = 5.2;
  constexpr double s_shellRxyBuffer = 2.5;
  bool shellClose = false;
  if (shellMode && r_xy > 1e-6) {
    double r_ideal = 0.0;
    double d_to_r_ideal = 0.0;
    bool haveHelixArc = false;
    const Layer* assocLayer = surface.associatedLayer();
    if (assocLayer != nullptr &&
        assocLayer->surfaceRepresentation().type() ==
            Surface::SurfaceType::Cylinder) {
      if constexpr (requires {
                      stepper.qOverP(state);
                      stepper.getField(state, position);
                    }) {
        auto fieldRes = stepper.getField(state, position);
        if (fieldRes.ok()) {
          const auto& cylBounds = static_cast<const CylinderBounds&>(
              assocLayer->surfaceRepresentation().bounds());
          r_ideal = cylBounds.get(CylinderBounds::eR);
          Vector3 landing;
          if (Acts::detail::helixBarrelCylinderLanding(
                  position, dir, stepper.qOverP(state), *fieldRes, r_ideal,
                  landing, &d_to_r_ideal)) {
            haveHelixArc = true;
          }
        }
      }
    }

    std::ostringstream layerIdOss;
    if (assocLayer != nullptr) {
      layerIdOss << assocLayer->surfaceRepresentation().geometryId();
    } else {
      layerIdOss << "<none>";
    }
    const std::string layerIdStr = layerIdOss.str();

    // Direction-aware "headed toward r_ideal" check. Mirrors Layer.cpp's
    // FAR-phase condition exactly so the two stay in lockstep:
    //   outgoing (pr > 0): trajectory below r_ideal by more than buffer
    //   incoming (pr < 0): trajectory above r_ideal by more than buffer
    // Absolute-value-only check would be unsafe: a trajectory past r_ideal
    // moving further away would have |r_xy − r_ideal| > buffer but
    // helixBarrelCylinderLanding returns the next-loop landing (large d),
    // which would trigger a catastrophic FAR step.
    const double prDir =
        (r_xy > 1.0e-6)
            ? (position[0] * dir[0] + position[1] * dir[1]) / r_xy
            : 0.0;
    const bool headedTowardRIdeal =
        (prDir > 0.0 && r_xy < r_ideal - s_shellRxyBuffer) ||
        (prDir < 0.0 && r_xy > r_ideal + s_shellRxyBuffer);
    if (haveHelixArc && d_to_r_ideal > s_shellProbeTrigger &&
        headedTowardRIdeal) {
      // Compute d* = helix arc to the buffer cylinder (r_ideal ± buffer on
      // the trajectory's side). Stepping by d* lands the trajectory exactly
      // on the close-zone boundary — no radial overshoot, regardless of
      // trajectory tangentiality. The path-limit cap on the FAR step
      // (~10 mm arc) can otherwise cause overshoot of the buffer when
      // d − safety_offset > path-limit; using d* lets the step land
      // precisely on the close zone whenever d* is the binding constraint.
      //
      // Fall back to (d − safety_offset) if the buffer-cylinder landing
      // fails (FP edge case — trajectory just past buffer cylinder).
      const double bufferTargetR =
          (prDir > 0.0) ? r_ideal - s_shellRxyBuffer
                        : r_ideal + s_shellRxyBuffer;
      double d_to_buffer = 0.0;
      bool haveBufferArc = false;
      if constexpr (requires {
                      stepper.qOverP(state);
                      stepper.getField(state, position);
                    }) {
        auto fieldRes2 = stepper.getField(state, position);
        if (fieldRes2.ok()) {
          Vector3 bufferLanding;
          if (Acts::detail::helixBarrelCylinderLanding(
                  position, dir, stepper.qOverP(state), *fieldRes2,
                  bufferTargetR, bufferLanding, &d_to_buffer) &&
              d_to_buffer > 0.0) {
            haveBufferArc = true;
          }
        }
      }
      const double stepCandSafety = d_to_r_ideal - s_shellSafetyOffset;
      const double stepFar = haveBufferArc
                                 ? std::min(stepCandSafety, d_to_buffer)
                                 : stepCandSafety;
      ACTS_VERBOSE("[shell] layer=" << layerIdStr
                   << " r_ideal=" << r_ideal
                   << " pos=(" << position[0] << "," << position[1] << ","
                   << position[2] << ") r_xy=" << r_xy
                   << " d_to_r_ideal=" << d_to_r_ideal
                   << " d_to_buffer=" << (haveBufferArc ? d_to_buffer : -1.0)
                   << " trigger=" << s_shellProbeTrigger
                   << " branch=FAR step_set=" << stepFar
                   << " (= min(d-" << s_shellSafetyOffset << ", d*))");
      stepper.releaseStepSize(state, stype);
      stepper.updateStepSize(state, stepFar, stype);
      return IntersectionStatus::reachable;
    }

    // CLOSE branch (or helix-arc unavailable — treat as close, let the
    // standard intersect decide).
    //
    // FAR-to-CLOSE handoff for the layer-rep carrier:
    //   When the target is SPECIFICALLY the layer's representing-surface
    //   cylinder (the sole carrier added to navSurfaces during FAR phase
    //   by Layer::compatibleSurfaces) AND we've just dropped below the
    //   FAR trigger (d <= probeTrigger), force unreachable. The navigator
    //   advances → navSurfaces exhausted → retarget → compatibleSurfaces
    //   re-runs at the now close-zone position. That re-entry sees d <= 3,
    //   enters CLOSE phase, and runs the sensor neighbors() lookup with
    //   the trajectory's actual close-zone phi/z — picking up the sensors
    //   the descending helix is heading into rather than the ones caught
    //   at the entry/apex position (event 12 failure mode in 0.24-0.26 GeV
    //   verbose trace).
    //
    // Identity check is strict: target surface address must equal the
    // layer's representing-surface address. Approach surfaces (apr=1,
    // apr=2) and boundary surfaces are ALSO non-sensitive but are
    // legitimate targets the trajectory must traverse — we MUST NOT mark
    // those unreachable here.
    if (haveHelixArc && !barrelSensitive && assocLayer != nullptr &&
        &surface == &assocLayer->surfaceRepresentation()) {
      ACTS_VERBOSE("[shell] layer=" << layerIdStr
                   << " r_ideal=" << r_ideal
                   << " pos=(" << position[0] << "," << position[1] << ","
                   << position[2] << ") r_xy=" << r_xy
                   << " d_to_r_ideal=" << d_to_r_ideal
                   << " branch=CLOSE-HANDOFF [layer-rep carrier done, "
                   "forcing retarget for close-zone neighbors() lookup]");
      return IntersectionStatus::unreachable;
    }

    // Predicate selection by surface kind:
    //   sensor surface → STRAIGHT LINE (raw tangent). intersectionDirection
    //     already holds dir * stepper.direction(state) from setup, so we
    //     leave it untouched for sensors.
    //   non-sensitive (approach / boundary) when helix-arc unavailable →
    //     pure radial inward (-r̂). Approach cylinders are radially
    //     well-conditioned in the barrel, so pure radial is the cleanest
    //     line-plane predicate, same as outside shell mode.
    shellClose = true;
    if (!barrelSensitive) {
      const double r_hat_x = position[0] / r_xy;
      const double r_hat_y = position[1] / r_xy;
      intersectionDirection[0] = -r_hat_x;
      intersectionDirection[1] = -r_hat_y;
      intersectionDirection[2] = 0.0;
    }
    ACTS_VERBOSE("[shell] layer=" << layerIdStr
                 << " r_ideal=" << r_ideal << " pos=(" << position[0]
                 << "," << position[1] << "," << position[2]
                 << ") r_xy=" << r_xy << " d_to_r_ideal=" << d_to_r_ideal
                 << " trigger=" << s_shellProbeTrigger
                 << " branch=CLOSE "
                 << (barrelSensitive ? "[tangent predicate, sensor]"
                                     : "[radial-inward predicate, non-sensitive]")
                 << (haveHelixArc ? "" : " (helix unavailable)")
                 << (!headedTowardRIdeal ? " (not headed toward r_ideal)" : ""));
  } else if (state.radiallyInward_previous && isInBarrelVolume &&
             r_xy > 1e-6 && !barrelSensitive) {
    // (b) non-sensitive: substitute pure radial inward.
    const double r_hat_x = position[0] / r_xy;
    const double r_hat_y = position[1] / r_xy;
    intersectionDirection[0] = -r_hat_x;
    intersectionDirection[1] = -r_hat_y;
    intersectionDirection[2] = 0.0;
    ACTS_VERBOSE(
        "Approach/boundary surface: using pure radial inward for intersect "
        "in barrel");
  }

  // Runtime A/B switch: set ACTS_DISABLE_HELIX_INTERSECT=1 in the
  // environment to force the line-plane fallback even on sensitive
  // inward-barrel surfaces. Evaluated once on first call.
  static const bool disableHelix =
      std::getenv("ACTS_DISABLE_HELIX_INTERSECT") != nullptr;

  // Compute the intersection.
  //   - barrelSensitive (both outside shell AND in shell
  //     CLOSE-branch) uses the closed-form helix-plane intersect with
  //     strict bounds. Forcing strict bounds here gives an authoritative
  //     bounded check on the curved landing — the navigator's candidate
  //     list deliberately uses lenient bounds upstream when building.
  //     The shell d-3 FAR branch above returns early (reachable
  //     with step = d−3) before reaching this lambda, so this lambda
  //     handles only the CLOSE-branch case for shell sensors.
  //     Using helixPlaneIntersection here (instead of a straight-line
  //     tangent intersect) gives the EXACT arc length to where the
  //     curve crosses the sensor plane, which eliminates the negative-
  //     pathLength ping-pong observed with the straight-line variant.
  //   - For non-sensitive surfaces in shell CLOSE-branch, the
  //     `intersectionDirection` was already substituted to pure radial
  //     inward above; we use that with surface.intersect, forced to
  //     strict bounds via the shellClose ternary.
  //   - All other regimes use surface.intersect with the navigator's
  //     boundaryTolerance.
  auto sIntersection = [&]() {
    if constexpr (requires {
                    stepper.qOverP(state);
                    stepper.getField(state, position);
                  }) {
      if (barrelSensitive && !disableHelix) {
        auto fieldRes = stepper.getField(state, position);
        if (fieldRes.ok()) {
          ACTS_VERBOSE(
              "Sensitive surface: using closed-form helix-plane intersect "
              "for barrel intersect"
              << (shellClose ? " (shell CLOSE)" : ""));
          return Acts::detail::helixPlaneIntersection(
              state.options.geoContext, surface, position, dir,
              stepper.qOverP(state), *fieldRes, BoundaryTolerance::None(),
              surfaceTolerance);
        }
      }
    }
    return surface
        .intersect(state.options.geoContext, stepper.position(state),
                   intersectionDirection,
                   shellClose ? BoundaryTolerance::None()
                                   : boundaryTolerance,
                   surfaceTolerance)[index];
  }();

  // The intersection is on surface already
  if (sIntersection.status() == IntersectionStatus::onSurface) {
    ACTS_VERBOSE("Intersection: state is ON SURFACE");
    if (shellClose) {
      ACTS_VERBOSE("[shell] close-branch surface=" << surface.geometryId()
                   << " pathLen=" << sIntersection.pathLength()
                   << " status=onSurface");
    }
    state.stepSize.release(stype);
    stepper.updateStepSize(state, sIntersection.pathLength(), stype);
    return IntersectionStatus::onSurface;
  }

  const double nearLimit = std::numeric_limits<double>::lowest();
  const double farLimit = std::numeric_limits<double>::max();

  if (sIntersection.isValid() &&
      detail::checkPathLength(sIntersection.pathLength(), nearLimit, farLimit,
                              logger)) {
    ACTS_VERBOSE("Surface is reachable");

    if (shellClose) {
      ACTS_VERBOSE("[shell] close-branch surface=" << surface.geometryId()
                   << " pathLen=" << sIntersection.pathLength()
                   << " status=reachable");
    }

    stepper.releaseStepSize(state, stype);
    stepper.updateStepSize(state, sIntersection.pathLength(), stype);
    return IntersectionStatus::reachable;
  }

  ACTS_VERBOSE("Surface is NOT reachable");
  if (shellClose) {
    ACTS_VERBOSE("[shell] close-branch surface=" << surface.geometryId()
                 << " status=unreachable");
  }
  return IntersectionStatus::unreachable;
}

}  // namespace Acts::detail
