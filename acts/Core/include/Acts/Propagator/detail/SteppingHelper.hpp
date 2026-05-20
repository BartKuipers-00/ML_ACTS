// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/HelixIntersection.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <cmath>
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
    ConstrainedStep::Type stype, bool isInBarrelVolume,
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
  //   (1) inward-barrel non-sensitive (approach / boundary / layer-rep)
  //       → pure −r̂ direction with the standard line predicate. These
  //       cylinders are radially well-conditioned in the barrel (≈⊥ r̂),
  //       giving the cleanest line-plane intersect.
  //   (2) inward-barrel sensitive → closed-form helix-plane intersect
  //       using the raw tangent, computed in the lambda below. Matches
  //       the actual curved landing on the module rectangle.
  //   (3) anything else (outward, non-barrel) → raw tangent + line.
  Vector3 intersectionDirection = direction * stepper.direction(state);
  const bool barrelSensitive =
      r_xy > 1e-6 && surface.geometryId().sensitive() != 0;

  if (state.radiallyInward_previous && isInBarrelVolume &&
      r_xy > 1e-6 && !barrelSensitive) {
    // Non-sensitive while moving inward: substitute pure radial inward.
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
  //   - barrelSensitive uses the closed-form helix-plane intersect with
  //     strict bounds. Forcing strict bounds here gives an authoritative
  //     bounded check on the curved landing — the navigator's candidate
  //     list deliberately uses lenient bounds upstream when building.
  //     Using helixPlaneIntersection here (instead of a straight-line
  //     tangent intersect) gives the EXACT arc length to where the curve
  //     crosses the sensor plane.
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
              "for barrel intersect");
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
    stepper.releaseStepSize(state, stype);
    stepper.updateStepSize(state, sIntersection.pathLength(), stype);
    return IntersectionStatus::reachable;
  }

  ACTS_VERBOSE("Surface is NOT reachable");
  return IntersectionStatus::unreachable;
}

}  // namespace Acts::detail
