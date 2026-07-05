// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {
class Surface;
}  // namespace Acts

#include <atomic>
#include <cstdint>

namespace Acts::detail {

/// [ckf-timing] diagnostic counters: total helixPlaneIntersection calls and
/// how many of those fell back to the straight-line intersect.
std::atomic<std::uint64_t>& helixIntersectCallCounter();
std::atomic<std::uint64_t>& helixIntersectFallbackCounter();

/// Closed-form helix-plane intersection for a charged particle in a uniform
/// magnetic field with a planar surface (e.g. a barrel sensitive module).
///
/// The trajectory is a perfect helix specified by `position`, `direction`,
/// signed inverse momentum `qOverP`, and the local magnetic field `bField`.
/// Currently supports the common axial-B / barrel-module case (B parallel to
/// the global z axis, surface normal in the transverse plane). For any other
/// configuration, or when the helix's circle does not reach the surface's
/// plane, the function returns `SurfaceIntersection::invalid()` so the caller
/// can fall back to the standard line-plane intersection.
///
/// The smaller of the two forward arc-length solutions is returned.
///
/// @param gctx              geometry context
/// @param surface           planar surface to intersect
/// @param position          current 3D position on the helix
/// @param direction         current unit tangent to the helix
/// @param qOverP            signed inverse momentum (1/p with charge sign)
/// @param bField            magnetic field vector at the current position
/// @param boundaryTolerance bounds-check policy for the planar landing point
/// @param surfaceTolerance  on-surface tolerance for the path-length test
SurfaceIntersection helixPlaneIntersection(
    const GeometryContext& gctx, const Surface& surface,
    const Vector3& position, const Vector3& direction, double qOverP,
    const Vector3& bField,
    const BoundaryTolerance& boundaryTolerance = BoundaryTolerance::None(),
    double surfaceTolerance = s_onSurfaceTolerance);

/// Compute the 3D landing point of the helix (parametrised by `position`,
/// `direction`, `qOverP`, axial `bField`) on a barrel cylinder of radius
/// `cylinderRadius` centred on the global z axis. Returns the smallest
/// positive forward arc-length solution. If the helix's transverse circle
/// does not reach the cylinder, returns false (and `landing` is unchanged).
///
/// Used by `Layer::compatibleSurfaces` to compute the lookupPosition for
/// `m_surfaceArray->neighbors(...)` so that the candidate bin is the one
/// the *curved* trajectory will enter, not the one a near-tangent line
/// extrapolation predicts.
/// @param arcLength  if non-null, set to the helix arc length from `position`
///                   to `landing` on success. (Used by the stepper in
///                   shell mode to size the step along the helix
///                   instead of along the tiny radial chord to a sensor's
///                   plane.)
bool helixBarrelCylinderLanding(const Vector3& position,
                                const Vector3& direction, double qOverP,
                                const Vector3& bField, double cylinderRadius,
                                Vector3& landing,
                                double* arcLength = nullptr);

}  // namespace Acts::detail
