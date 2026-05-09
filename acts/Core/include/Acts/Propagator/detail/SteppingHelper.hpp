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
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"

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
    ConstrainedStep::Type stype, bool isInBarrelVolume, const Logger& logger) {
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

  // Determine which direction to use for intersection calculation
  Vector3 intersectionDirection = direction * stepper.direction(state);

  // When going radially inward in barrel regions, substitute the
  // intersection direction. The substitution depends on the target surface:
  //   - sensitive sensor surfaces: bisector(direction, -r̂). Keeps a
  //     well-conditioned radial component for the line-plane denominator
  //     while preserving forward motion, so the per-step intersect lands
  //     on the actual module the trajectory is heading toward (raw tangent
  //     mis-targets across phi-bins on the return arc).
  //   - approach / boundary / layer-rep surfaces: pure -r̂. These planes
  //     are radially well-conditioned (≈⊥ r̂ in the barrel) so pure radial
  //     gives the cleanest intersect; the bisector here can mis-pick
  //     across layers near the apex where tangent ⊥ radial-inward.
  if (state.radiallyInward_previous && isInBarrelVolume && r_xy > 1e-6) {
    const double r_hat_x = position[0] / r_xy;
    const double r_hat_y = position[1] / r_xy;
    const bool isSensitive = surface.geometryId().sensitive() != 0;

    if (isSensitive) {
      const double sum_x = intersectionDirection[0] - r_hat_x;
      const double sum_y = intersectionDirection[1] - r_hat_y;
      const double sum_z = intersectionDirection[2];
      const double sum_mag =
          std::sqrt(sum_x * sum_x + sum_y * sum_y + sum_z * sum_z);
      if (sum_mag > 1e-9) {
        intersectionDirection[0] = sum_x / sum_mag;
        intersectionDirection[1] = sum_y / sum_mag;
        intersectionDirection[2] = sum_z / sum_mag;
        ACTS_VERBOSE(
            "Sensitive surface: using bisector(tangent, radial-inward) for "
            "intersect in barrel");
      } else {
        intersectionDirection[0] = -r_hat_x;
        intersectionDirection[1] = -r_hat_y;
        intersectionDirection[2] = 0.0;
        ACTS_VERBOSE(
            "Sensitive surface: bisector degenerate — falling back to radial "
            "inward");
      }
    } else {
      intersectionDirection[0] = -r_hat_x;
      intersectionDirection[1] = -r_hat_y;
      intersectionDirection[2] = 0.0;
      ACTS_VERBOSE(
          "Approach/boundary surface: using pure radial inward for intersect "
          "in barrel");
    }
  }

  auto sIntersection =
      surface.intersect(state.options.geoContext, stepper.position(state),
                        intersectionDirection, boundaryTolerance,
                        surfaceTolerance)[index];

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
