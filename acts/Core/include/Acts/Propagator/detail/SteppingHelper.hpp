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

#include <cstdlib>
#include <limits>

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
    ConstrainedStep::Type stype, bool isInBarrelVolume, bool radialDownMode,
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

  // Determine which direction / predicate to use for intersection calculation.
  // Two regimes for inward-barrel motion:
  //   (a) sensitive sensor surfaces  → closed-form helix-plane intersect,
  //       using the *raw* tangent direction. The helix predicate matches
  //       the actual curved trajectory's bounded landing on the module's
  //       rectangle, which a line predicate (any line, including the
  //       bisector substitution) cannot do near the apex where the curve
  //       and its tangent disagree at sub-mm scales. Falls back to the
  //       line predicate internally if the helix preconditions don't hold
  //       (non-axial B, tilted surface, degenerate kinematics).
  //   (b) approach / boundary / layer-representing surfaces → pure −r̂
  //       direction with the standard line predicate. These planes are
  //       radially well-conditioned in the barrel (≈⊥ r̂), so pure radial
  //       gives the cleanest line-plane intersect.
  // Outside the inward-barrel regime, we use the raw stepper tangent and
  // the standard line predicate everywhere.
  Vector3 intersectionDirection = direction * stepper.direction(state);
  // The helix-plane predicate is direction-agnostic (closed-form arc-length
  // to a plane works for outward and inward motion). Use it for ALL
  // sensitive surfaces in the barrel — both arcs of a looping trajectory
  // benefit, mitigating the apex-inside-shell failure mode that drops
  // outgoing AND incoming hits symmetrically.
  const bool barrelSensitive =
      isInBarrelVolume && r_xy > 1e-6 &&
      surface.geometryId().sensitive() != 0;

  // Radial-down mode forces the line predicate with pure radial-inward
  // direction for ALL surfaces (sensitive included). No helix in this
  // regime — the navigator just looks down at the next module each step.
  // Outside radial-down mode, the old logic remains: pure-radial-inward
  // for non-sensitive in inward barrel, helix for sensitive in barrel.
  if (radialDownMode && r_xy > 1e-6) {
    const double r_hat_x = position[0] / r_xy;
    const double r_hat_y = position[1] / r_xy;
    intersectionDirection[0] = -r_hat_x;
    intersectionDirection[1] = -r_hat_y;
    intersectionDirection[2] = 0.0;
    ACTS_VERBOSE(
        "Radial-down mode: using pure radial inward for intersect");
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

  // Compute the intersection. Helix predicate when (a) holds and the stepper
  // exposes q/p and field access; otherwise the standard surface.intersect.
  auto sIntersection = [&]() {
    if constexpr (requires {
                    stepper.qOverP(state);
                    stepper.getField(state, position);
                  }) {
      if (barrelSensitive && !disableHelix && !radialDownMode) {
        auto fieldRes = stepper.getField(state, position);
        if (fieldRes.ok()) {
          ACTS_VERBOSE(
              "Sensitive surface: using closed-form helix-plane intersect "
              "for barrel intersect");
          // Override the candidate's boundaryTolerance to None for the
          // helix's bounded check. Layer::compatibleSurfaces deliberately
          // tags sensitive candidates with Infinite so the list is built
          // leniently; here we want the strict bounds test against the
          // actual *curved* landing point.
          return Acts::detail::helixPlaneIntersection(
              state.options.geoContext, surface, position, dir,
              stepper.qOverP(state), *fieldRes, BoundaryTolerance::None(),
              surfaceTolerance);
        }
      }
    }
    return surface
        .intersect(state.options.geoContext, stepper.position(state),
                   intersectionDirection, boundaryTolerance,
                   surfaceTolerance)[index];
  }();

  // The intersection is on surface already
  if (sIntersection.status() == IntersectionStatus::onSurface) {
    ACTS_VERBOSE("Intersection: state is ON SURFACE");
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

    // Radial-down mode: override the step size with the helix arc length to
    // where the trajectory crosses the layer's r_m cylinder. The line +
    // radial-inward intersect gave us the bounded-check verdict (reachable
    // / onSurface / unreachable) and its pathLength is the radial distance
    // to the sensor's plane — a few mm at apex. Taking that as the step
    // size causes tangentially-moving trajectories at apex to make near-
    // zero r-progress per step, triggering hundreds of TP oscillations
    // before the trajectory descends. Using the helix arc length to the
    // r_m cylinder ensures each step advances meaningfully along the arc.
    double stepPath = sIntersection.pathLength();
    if (radialDownMode) {
      if constexpr (requires {
                      stepper.qOverP(state);
                      stepper.getField(state, position);
                    }) {
        auto fieldRes = stepper.getField(state, position);
        if (fieldRes.ok()) {
          const Layer* layer = surface.associatedLayer();
          if (layer != nullptr &&
              layer->surfaceRepresentation().type() ==
                  Surface::SurfaceType::Cylinder) {
            const auto& cylBounds = static_cast<const CylinderBounds&>(
                layer->surfaceRepresentation().bounds());
            const double rLayer = cylBounds.get(CylinderBounds::eR);
            Vector3 helixLanding;
            double arcLen = 0.0;
            if (detail::helixBarrelCylinderLanding(
                    position, dir, stepper.qOverP(state), *fieldRes, rLayer,
                    helixLanding, &arcLen) &&
                arcLen > stepPath) {
              ACTS_VERBOSE(
                  "Radial-down mode: stepping helix arc to r_m cylinder ("
                  << arcLen << " mm) instead of radial chord (" << stepPath
                  << " mm)");
              stepPath = arcLen;
            }
          }
        }
      }
    }

    stepper.releaseStepSize(state, stype);
    stepper.updateStepSize(state, stepPath, stype);
    return IntersectionStatus::reachable;
  }

  ACTS_VERBOSE("Surface is NOT reachable");
  return IntersectionStatus::unreachable;
}

}  // namespace Acts::detail
