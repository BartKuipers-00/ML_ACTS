// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/HelixIntersection.hpp"

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/Result.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <utility>

namespace Acts::detail {

namespace {

// Wrap a delta-α into a forward branch consistent with the rotation sign.
// rotSign = +1: forward motion ⇒ α increases ⇒ desired Δα ∈ [0, 2π).
// rotSign = -1: forward motion ⇒ α decreases ⇒ desired Δα ∈ (-2π, 0].
constexpr double wrapForward(double dalpha, double rotSign) {
  constexpr double twoPi = 2.0 * std::numbers::pi;
  if (rotSign > 0.0) {
    while (dalpha < 0.0) {
      dalpha += twoPi;
    }
    while (dalpha >= twoPi) {
      dalpha -= twoPi;
    }
  } else {
    while (dalpha > 0.0) {
      dalpha -= twoPi;
    }
    while (dalpha <= -twoPi) {
      dalpha += twoPi;
    }
  }
  return dalpha;
}

}  // namespace

SurfaceIntersection helixPlaneIntersection(
    const GeometryContext& gctx, const Surface& surface,
    const Vector3& position, const Vector3& direction, double qOverP,
    const Vector3& bField, const BoundaryTolerance& boundaryTolerance,
    double surfaceTolerance) {
  // Helper: when preconditions for the closed-form helix model are not
  // satisfied (non-axial B, degenerate kinematics, tilted surface), we
  // fall back to the standard line-plane intersect so the caller always
  // gets one authoritative answer.
  auto lineIntersect = [&]() {
    return surface
        .intersect(gctx, position, direction, boundaryTolerance,
                   surfaceTolerance)
        .closest();
  };

  // Guard against degenerate inputs that would make the helix model
  // meaningless.
  const double bMag = bField.norm();
  if (bMag < 1.0e-12 || std::abs(qOverP) < 1.0e-12) {
    return lineIntersect();
  }

  // Restrict to axial B (B parallel to global z). For non-axial fields,
  // bail out — the closed-form below is specific to barrel + axial.
  const double bAxialFraction = std::abs(bField[2]) / bMag;
  if (bAxialFraction < 0.999) {
    return lineIntersect();
  }
  const double bz = bField[2];

  // Decompose direction into transverse and z components. cos(λ) is the
  // transverse magnitude; sin(λ) the z component. direction is unit-norm.
  const double dT_sq =
      direction[0] * direction[0] + direction[1] * direction[1];
  if (dT_sq < 1.0e-12) {
    // Direction parallel to z axis ⇒ no transverse circle; fall back.
    return lineIntersect();
  }
  const double cosLambda = std::sqrt(dT_sq);
  const double sinLambda = direction[2];

  // Helix radius and rotation sign in the transverse plane.
  //   p     = 1 / |qOverP|       (native units)
  //   p_T   = p · cos(λ)
  //   R     = p_T / |q · bz|     (in ACTS native units)
  // Using the identity |q|/|qOverP| = p, the unsigned radius simplifies to
  //   R = cos(λ) / (|qOverP| · |bz|).
  const double R = cosLambda / (std::abs(qOverP) * std::abs(bz));

  // Sign convention. F = q v × B with B = bz ẑ gives F_xy =
  // q·bz·(v_y, -v_x, 0). The trajectory curves toward +(v_y, -v_x) when
  // q·bz > 0, which puts the helix center to the right of v_T. Let
  //   signC = sign(qOverP · bz)
  // be the sign that puts the center at position + signC · R · (d̂_y, -d̂_x).
  // The angular coordinate α (around the center) then evolves with
  //   sign(dα/ds) = -signC ≡ rotSign.
  const double signC = (qOverP * bz > 0.0) ? 1.0 : -1.0;
  const double rotSign = -signC;

  // Helix center in xy.
  const double Cx = position[0] + signC * R * direction[1] / cosLambda;
  const double Cy = position[1] - signC * R * direction[0] / cosLambda;

  // Initial angular coordinate.
  const double alpha0 = std::atan2(position[1] - Cy, position[0] - Cx);

  // Surface normal at the current point. Restrict to barrel-like (n̂_z ≈ 0).
  const Vector3 n3 = surface.normal(gctx, position, direction);
  const double nT_sq = n3[0] * n3[0] + n3[1] * n3[1];
  if (nT_sq < 1.0e-9 || std::abs(n3[2]) > 0.05) {
    // Either degenerate transverse normal or significantly tilted out of
    // the (x,y) plane. Closed form below assumes barrel orientation.
    return lineIntersect();
  }
  const double nT_inv = 1.0 / std::sqrt(nT_sq);
  const double nx = n3[0] * nT_inv;
  const double ny = n3[1] * nT_inv;

  // The plane equation in the transverse plane is n̂ · X_xy = d_plane,
  // where d_plane = n̂ · P_S_xy (P_S = surface center). Because n̂_z ≈ 0,
  // this captures the planar surface in 3D as long as the surface extends
  // along z (true for barrel modules).
  const Vector3 sCenter = surface.center(gctx);
  const double d_plane = nx * sCenter[0] + ny * sCenter[1];

  // Helix-plane intersection condition reduces to
  //   n_x · cos(α) + n_y · sin(α) = (d_plane - n̂ · C_xy) / R
  // with |left side| ≤ 1 (since n̂_xy is unit). Real solutions iff |rhs| ≤ 1.
  const double rhs = (d_plane - (nx * Cx + ny * Cy)) / R;
  if (std::abs(rhs) > 1.0) {
    // Helix's circle never reaches the module's plane (apex below the
    // module's r). Genuinely unreachable — not a fallback.
    return SurfaceIntersection::invalid(surface);
  }

  // Solve.
  const double phi_n = std::atan2(ny, nx);
  const double dphi = std::acos(std::clamp(rhs, -1.0, 1.0));
  const double alpha1 = phi_n + dphi;
  const double alpha2 = phi_n - dphi;

  // Convert each candidate α to a forward arc length. We parametrize by
  // total arc length s; angular advance per unit s is rotSign · cos(λ) / R.
  const double dalpha_to_s = R / (rotSign * cosLambda);

  const double dalpha1 = wrapForward(alpha1 - alpha0, rotSign);
  const double dalpha2 = wrapForward(alpha2 - alpha0, rotSign);

  const double s1 = dalpha1 * dalpha_to_s;  // ≥ 0
  const double s2 = dalpha2 * dalpha_to_s;  // ≥ 0

  // Pick the smaller forward arc length. (s1, s2 are both ≥ 0 by wrapForward.)
  // We do NOT try the far solution if the near one fails the bounded check —
  // the far solution is on the *other side* of the helix loop (after a long
  // arc), where energy losses and material make the local-helix model
  // unreliable, and the trajectory will have crossed many other surfaces
  // first. If the near landing isn't on the module's silicon, the curve
  // physically doesn't hit the module on its forward path.
  const bool pickFirst = (s1 <= s2);
  const double sChosen = pickFirst ? s1 : s2;
  const double alphaChosen = pickFirst ? alpha1 : alpha2;

  const Vector3 landing{Cx + R * std::cos(alphaChosen),
                        Cy + R * std::sin(alphaChosen),
                        position[2] + sChosen * sinLambda};

  IntersectionStatus status = IntersectionStatus::reachable;
  if (sChosen < surfaceTolerance) {
    status = IntersectionStatus::onSurface;
  }
  if (!boundaryTolerance.isInfinite()) {
    auto localResult = surface.globalToLocal(gctx, landing, direction);
    if (!localResult.ok() ||
        !surface.bounds().inside(localResult.value(), boundaryTolerance)) {
      return SurfaceIntersection::invalid(surface);
    }
  }

  return SurfaceIntersection{Intersection3D{landing, sChosen, status}, surface,
                             0, boundaryTolerance};
}

bool helixBarrelCylinderLanding(const Vector3& position,
                                const Vector3& direction, double qOverP,
                                const Vector3& bField, double cylinderRadius,
                                Vector3& landing, double* arcLength) {
  if (cylinderRadius <= 0.0) {
    return false;
  }
  const double bMag = bField.norm();
  if (bMag < 1.0e-12 || std::abs(qOverP) < 1.0e-12) {
    return false;
  }
  // Axial-B only (matches the helixPlaneIntersection precondition).
  if (std::abs(bField[2]) / bMag < 0.999) {
    return false;
  }
  const double bz = bField[2];

  const double dT_sq =
      direction[0] * direction[0] + direction[1] * direction[1];
  if (dT_sq < 1.0e-12) {
    return false;
  }
  const double cosLambda = std::sqrt(dT_sq);
  const double sinLambda = direction[2];

  const double R = cosLambda / (std::abs(qOverP) * std::abs(bz));
  const double signC = (qOverP * bz > 0.0) ? 1.0 : -1.0;
  const double rotSign = -signC;

  const double Cx = position[0] + signC * R * direction[1] / cosLambda;
  const double Cy = position[1] - signC * R * direction[0] / cosLambda;

  // The helix's transverse circle (radius R, center (Cx, Cy)) intersects
  // the cylinder x² + y² = R_layer² when the equation
  //   (Cx + R cos α)² + (Cy + R sin α)² = R_layer²
  // has a real solution. Expanding gives
  //   Cx · cos α + Cy · sin α = (R_layer² - Cx² - Cy² - R²) / (2 R).
  const double Cmag2 = Cx * Cx + Cy * Cy;
  const double rhs = (cylinderRadius * cylinderRadius - Cmag2 - R * R) /
                     (2.0 * R);
  const double rhsScale = std::sqrt(Cmag2);
  if (rhsScale < 1.0e-12) {
    // Helix center on the z-axis: every α is at the same r. Either the
    // cylinder coincides with the helix circle (ill-defined) or never
    // intersects. Bail out.
    return false;
  }
  const double rhsNorm = rhs / rhsScale;
  if (std::abs(rhsNorm) > 1.0) {
    // Helix circle and cylinder do not intersect.
    return false;
  }

  // Solve A · cos α + B · sin α = rhs, with A = Cx, B = Cy.
  const double phi_n = std::atan2(Cy, Cx);
  const double dphi = std::acos(std::clamp(rhsNorm, -1.0, 1.0));
  const double alpha1 = phi_n + dphi;
  const double alpha2 = phi_n - dphi;

  const double alpha0 = std::atan2(position[1] - Cy, position[0] - Cx);
  const double dalpha_to_s = R / (rotSign * cosLambda);

  const double dalpha1 = wrapForward(alpha1 - alpha0, rotSign);
  const double dalpha2 = wrapForward(alpha2 - alpha0, rotSign);

  const double s1 = dalpha1 * dalpha_to_s;
  const double s2 = dalpha2 * dalpha_to_s;

  const bool pickFirst = (s1 <= s2);
  const double sChosen = pickFirst ? s1 : s2;
  const double alphaChosen = pickFirst ? alpha1 : alpha2;

  landing = Vector3{Cx + R * std::cos(alphaChosen),
                    Cy + R * std::sin(alphaChosen),
                    position[2] + sChosen * sinLambda};
  if (arcLength != nullptr) {
    *arcLength = sChosen;
  }
  return true;
}

}  // namespace Acts::detail
