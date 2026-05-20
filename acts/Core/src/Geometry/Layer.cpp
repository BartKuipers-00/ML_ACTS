// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/Layer.hpp"

#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/HelixIntersection.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <iostream>
#include "Acts/Utilities/StringHelpers.hpp"

#include <algorithm>
#include <vector>

namespace Acts {

Layer::Layer(std::unique_ptr<SurfaceArray> surfaceArray, double thickness,
             std::unique_ptr<ApproachDescriptor> ades, LayerType laytyp)
    : m_nextLayers(NextLayers(nullptr, nullptr)),
      m_surfaceArray(surfaceArray.release()),
      m_layerThickness(thickness),
      m_approachDescriptor(nullptr),
      m_representingVolume(nullptr),
      m_layerType(laytyp),
      m_ssRepresentingSurface(1) {
  if (ades) {
    ades->registerLayer(*this);
    m_approachDescriptor = std::move(ades);
    m_ssApproachSurfaces = 1;  // indicates existence
  }
  // indicates existence of sensitive surfaces
  if (m_surfaceArray) {
    m_ssSensitiveSurfaces = 1;
  }
}

Layer::~Layer() noexcept = default;

const ApproachDescriptor* Layer::approachDescriptor() const {
  return m_approachDescriptor.get();
}

ApproachDescriptor* Layer::approachDescriptor() {
  return const_cast<ApproachDescriptor*>(m_approachDescriptor.get());
}

void Layer::closeGeometry(const IMaterialDecorator* materialDecorator,
                          const GeometryIdentifier& layerID,
                          const GeometryIdentifierHook& hook,
                          const Logger& logger) {
  // set the volumeID of this
  assignGeometryId(layerID);
  // assign to the representing surface
  Surface* rSurface = const_cast<Surface*>(&surfaceRepresentation());
  if (materialDecorator != nullptr) {
    materialDecorator->decorate(*rSurface);
  }
  ACTS_DEBUG("layerID: " << layerID);

  rSurface->assignGeometryId(layerID);

  // also find out how the sub structure is defined
  if (surfaceRepresentation().surfaceMaterial() != nullptr) {
    m_ssRepresentingSurface = 2;
  }
  // loop over the approach surfaces
  if (m_approachDescriptor) {
    // indicates the existence of approach surfaces
    m_ssApproachSurfaces = 1;
    // loop through the approachSurfaces and assign unique GeomeryID
    GeometryIdentifier::Value iasurface = 0;
    for (auto& aSurface : m_approachDescriptor->containedSurfaces()) {
      auto asurfaceID = GeometryIdentifier(layerID).withApproach(++iasurface);
      auto mutableASurface = const_cast<Surface*>(aSurface);
      mutableASurface->assignGeometryId(asurfaceID);
      if (materialDecorator != nullptr) {
        materialDecorator->decorate(*mutableASurface);
      }
      // if any of the approach surfaces has material
      if (aSurface->surfaceMaterial() != nullptr) {
        m_ssApproachSurfaces = 2;
      }
    }
  }
  // check if you have sensitive surfaces
  if (m_surfaceArray) {
    // indicates the existence of sensitive surfaces
    m_ssSensitiveSurfaces = 1;
    // loop sensitive surfaces and assign unique GeometryIdentifier
    GeometryIdentifier::Value issurface = 0;
    for (auto& sSurface : m_surfaceArray->surfaces()) {
      auto ssurfaceID = GeometryIdentifier(layerID).withSensitive(++issurface);
      ssurfaceID = hook.decorateIdentifier(ssurfaceID, *sSurface);
      auto mutableSSurface = const_cast<Surface*>(sSurface);
      mutableSSurface->assignGeometryId(ssurfaceID);
      if (materialDecorator != nullptr) {
        materialDecorator->decorate(*mutableSSurface);
      }
      // if any of the sensitive surfaces has material
      if (sSurface->surfaceMaterial() != nullptr) {
        m_ssSensitiveSurfaces = 2;
      }
    }
  }
}

boost::container::small_vector<SurfaceIntersection, 10>
Layer::compatibleSurfaces(const GeometryContext& gctx, const Vector3& position,
                          const Vector3& direction,
                          const NavigationOptions<Surface>& options) const {
  // the list of valid intersection
  boost::container::small_vector<SurfaceIntersection, 10> sIntersections;

  // fast exit - there is nothing to
  if (!m_surfaceArray || !m_approachDescriptor) {
    return sIntersections;
  }

  double nearLimit = options.nearLimit;
  double farLimit = options.farLimit;

  // True iff the navigator has populated helix info (i.e. we're in the
  // inward-barrel regime where the trajectory's curvature significantly
  // differs from a line over the candidate-list construction's relevant
  // path lengths). Used both for the bin-lookup raycast and the per-surface
  // intersect predicate selection.
  const bool helixRegime = options.helixBField.norm() > 0.0 &&
                           options.helixQOverP != 0.0 &&
                           options.helixDirection.squaredNorm() > 0.0;

  // Shell phase: FAR vs CLOSE based on helix arc d to r_ideal.
  //
  // FAR (d > 3 mm of helix arc to the layer's representing-surface
  // cylinder): defer the sensor neighbors() lookup; the navSurfaces list
  // is populated with only the layer's representing cylinder as a
  // "carrier" target. The SteppingHelper FAR branch steps (d - 3) along
  // the curve and then forces a retarget so compatibleSurfaces is
  // re-invoked at the close-zone position.
  //
  // CLOSE (d <= 3 OR helix can't reach r_ideal): run the standard sensor
  // neighbors() + helix-plane lookup. This is the actual sensor probe,
  // running with the trajectory at most 3 mm above r_ideal so the
  // bin-projected lookupPosition matches where the helix will land.
  //
  // Effect: neighbors() is called exactly once per layer-entry in the
  // close zone, regardless of how far above r_ideal the trajectory
  // entered the shell. Sensors that the FAR descent drifted-in toward
  // are now caught by the close-zone bin lookup, instead of relying on
  // the entry-position bin which can phi-drift out of range during the
  // descent (event 12 failure mode in 0.24-0.26 GeV verbose trace).
  // PROBE: when true, skip the shell-mode FAR descent entirely. Section
  // (B) (sensor neighbors lookup with helix-cylinder-landing projection)
  // runs immediately on shell-mode entry, and the layer-rep carrier is
  // never added (so the SteppingHelper FAR-step block stays in place but
  // never fires — no target ever points at the layer-rep cylinder).
  //
  // Motivation: for strip layers the apr envelope and the silicon shell
  // coincide (apr half-thickness ≈ shell half-thickness ≈ 5–6 mm), so
  // there is no clean "before silicon" position for the FAR descent to
  // park at. This probe answers whether the close-zone retarget alone
  // (with the helix-plane intersect for each sensor, still using the
  // helix-landing projection for neighbors()) is sufficient to recover
  // the missing hits, without any FAR-phase intermediate.
  constexpr bool kShellSkipFarPhase = false;
  bool shellFar = false;
  if (!kShellSkipFarPhase && options.shellMode && helixRegime &&
      surfaceRepresentation().type() == Surface::SurfaceType::Cylinder) {
    const auto& cylBounds = static_cast<const CylinderBounds&>(
        surfaceRepresentation().bounds());
    const double r_ideal = cylBounds.get(CylinderBounds::eR);
    Vector3 landing;
    double d_to_r_ideal = 0.0;
    if (detail::helixBarrelCylinderLanding(
            position, options.helixDirection, options.helixQOverP,
            options.helixBField, r_ideal, landing, &d_to_r_ideal)) {
      // FAR phase conditions match SteppingHelper's FAR-step gate exactly.
      // Symmetric for incoming and outgoing motion:
      //   incoming (pr < 0):  d > probeTrigger AND rxy > r_ideal + buffer
      //   outgoing (pr > 0):  d > probeTrigger AND rxy < r_ideal - buffer
      // The radial-momentum sign (pr) selects which side of r_ideal the
      // trajectory is approaching from; absolute-value-only check is
      // unsafe — a trajectory past r_ideal going further away would
      // re-trigger FAR with a huge d from the next-loop landing.
      // Strict alignment ("if and only if") with SteppingHelper guarantees
      // consistency: Layer.cpp says FAR ⇒ SteppingHelper's FAR step
      // fires, so the trajectory makes progress.
      constexpr double kShellProbeTrigger = 5.2;   // mm of helix arc
      constexpr double kShellRxyBuffer    = 2.5;   // mm
      const double rxyAtPos = std::sqrt(position[0] * position[0] +
                                         position[1] * position[1]);
      const double prDir =
          (rxyAtPos > 1.0e-6)
              ? (position[0] * options.helixDirection[0] +
                 position[1] * options.helixDirection[1]) / rxyAtPos
              : 0.0;
      const bool headedTowardRIdeal =
          (prDir > 0.0 && rxyAtPos < r_ideal - kShellRxyBuffer) ||
          (prDir < 0.0 && rxyAtPos > r_ideal + kShellRxyBuffer);
      if (d_to_r_ideal > kShellProbeTrigger && headedTowardRIdeal) {
        shellFar = true;
      }
    }
  }

  auto isUnique = [&](const SurfaceIntersection& b) {
    return std::ranges::none_of(sIntersections, [&b](const auto& a) {
      return &a.surface() == &b.surface() && a.index() == b.index();
    });
  };

  // lemma 0 : accept the surface
  auto acceptSurface = [&options](const Surface& sf,
                                  bool sensitive = false) -> bool {
    // surface is sensitive and you're asked to resolve
    if (sensitive && options.resolveSensitive) {
      return true;
    }
    // next option: it's a material surface, and you want to have it
    if (options.resolveMaterial && sf.surfaceMaterial() != nullptr) {
      return true;
    }
    // last option: resolve all
    return options.resolvePassive;
  };

  // lemma 1 : check and fill the surface
  // [&sIntersections, &options, &parameters
  auto processSurface = [&](const Surface& sf, bool sensitive = false) {
    // veto if it's start surface
    if (options.startObject == &sf) {
      return;
    }
    // veto if it doesn't fit the prescription
    if (!acceptSurface(sf, sensitive)) {
      return;
    }
    BoundaryTolerance boundaryTolerance = options.boundaryTolerance;
    if (rangeContainsValue(options.externalSurfaces, sf.geometryId())) {
      boundaryTolerance = BoundaryTolerance::Infinite();
    }
    // Predicate selection in compatibleSurfaces (four regimes):
    //   (1) shellMode + sensitive → line predicate with BISECTOR
    //       direction (normalize(tangent + (-r̂))) and strict bounds.
    //       This branch is reached only when the SteppingHelper per-step
    //       gate determined the trajectory is within ~1.5 mm of r_ideal.
    //   (2) helixRegime + sensitive (no shell) → closed-form
    //       helix-plane intersect using the raw tangent, bounded check
    //       on the actual curved landing.
    //   (3) helixRegime + non-sensitive (inward motion) → line predicate
    //       with PURE radial inward direction. Approach / boundary /
    //       layer-rep cylinders are radially well-conditioned in the
    //       barrel, so pure radial gives the cleanest intersect.
    //   (4) anything else (outward, non-barrel) → line with the supplied
    //       `direction` (raw tangent).
    // The path-length mismatch between helix-arc (sensors) and line-radial
    // (approaches) is resolved by the sensor-priority sort in
    // Navigator::resolveSurfaces: sensors come before approach surfaces in
    // navSurfaces regardless of path length.
    SurfaceIntersection sfi = SurfaceIntersection::invalid();
    // Compute radial motion sign at this position+direction. The pure-
    // radial-inward substitution for non-sensitive surfaces makes sense
    // only when the trajectory is heading inward; for outward motion the
    // raw tangent is correct.
    const double rxyAtPos = std::sqrt(position[0] * position[0] +
                                       position[1] * position[1]);
    const double prDir =
        (rxyAtPos > 1.0e-6)
            ? (position[0] * direction[0] + position[1] * direction[1]) /
                  rxyAtPos
            : 0.0;
    const bool inwardMotionLocal = prDir < 0.0;
    if (options.shellMode && sensitive && rxyAtPos > 1.0e-6) {
      // Shell mode for sensors: use the closed-form HELIX-PLANE
      // intersect (same predicate as the standard outward sensor path).
      // The straight-line tangent and bisector variants both produced
      // ping-pong loops on near-apex states: the line extrapolation
      // diverges from the curve over ~5–10 mm, predicting negative
      // pathLength to sensors the curve will actually reach by curving
      // forward. helixPlaneIntersection returns the exact arc length
      // to the curve's landing on the sensor's plane plus a strict
      // bounded check against the silicon. Strict bounds forced here
      // (the navigator candidate list uses lenient bounds upstream).
      //
      // STRAIGHT-LINE / BISECTOR variants kept commented for revert:
      //   sfi = sf.intersect(gctx, position, direction,
      //                      BoundaryTolerance::None()).closest();
      //   // bisector:
      //   const double r_hat_x = position[0] / rxyAtPos;
      //   const double r_hat_y = position[1] / rxyAtPos;
      //   const double sum_x = direction[0] - r_hat_x;
      //   const double sum_y = direction[1] - r_hat_y;
      //   const double sum_z = direction[2];
      //   const double sum_mag = std::sqrt(sum_x*sum_x + sum_y*sum_y +
      //                                    sum_z*sum_z);
      //   Vector3 dirToUse =
      //       (sum_mag > 1.0e-9)
      //           ? Vector3{sum_x/sum_mag, sum_y/sum_mag, sum_z/sum_mag}
      //           : Vector3{-r_hat_x, -r_hat_y, 0.0};
      //   sfi = sf.intersect(gctx, position, dirToUse,
      //                      BoundaryTolerance::None()).closest();
      sfi = detail::helixPlaneIntersection(
          gctx, sf, position, options.helixDirection, options.helixQOverP,
          options.helixBField, BoundaryTolerance::None());
    } else if (helixRegime && sensitive) {
      sfi = detail::helixPlaneIntersection(
          gctx, sf, position, options.helixDirection, options.helixQOverP,
          options.helixBField, boundaryTolerance);
    } else if (helixRegime && inwardMotionLocal) {
      // Pure radial inward (only for non-sensitive in inward motion).
      const Vector3 inwardDir{-position[0] / rxyAtPos,
                              -position[1] / rxyAtPos, 0.0};
      sfi = sf.intersect(gctx, position, inwardDir, boundaryTolerance).closest();
    } else {
      sfi = sf.intersect(gctx, position, direction, boundaryTolerance).closest();
    }

    if (!sfi.isValid()) {
      return;
    }

    bool pathLengthOk = detail::checkPathLength(sfi.pathLength(), nearLimit, farLimit);
    if (!pathLengthOk) {
      return;
    }

    if (!isUnique(sfi)) {
      return;
    }

    sIntersections.push_back(sfi);
  };

  // (A) approach descriptor section
  //
  // the approach surfaces are in principle always testSurfaces
  // - the surface on approach is excluded via the veto
  // - the surfaces are only collected if needed
  //
  // SKIP this section when shell mode is active. In that mode the
  // SteppingHelper FAR-branch already steps the helix arc toward
  // r_ideal; including the approach surfaces here would let the
  // navigator switch its target away from the sensor list and abort
  // the per-step bisector probing. apr=1 re-enters the candidate pool
  // via the standard layerTarget path once shell mode exits.
  if (m_approachDescriptor && !options.shellMode &&
      (options.resolveMaterial || options.resolvePassive)) {
    // the approach surfaces
    const std::vector<const Surface*>& approachSurfaces =
        m_approachDescriptor->containedSurfaces();
    // we loop through and veto
    // - if the approach surface is the parameter surface
    // - if the surface is not compatible with the collect
    for (auto& aSurface : approachSurfaces) {
      processSurface(*aSurface);
    }
  }

  // (B) sensitive surface section
  //
  // check the sensitive surfaces if you have some.
  //
  // SKIP in shell FAR phase: neighbors() is deferred until the
  // trajectory descends into the close zone (handled by the retarget
  // path after the FAR-step + layer-rep unreachable handoff).
  if (m_surfaceArray && !shellFar &&
      (options.resolveMaterial || options.resolvePassive ||
                         options.resolveSensitive)) {
    // ── PATCH: at()-based lookup (currently DISABLED — upstream-clean state)
    // Project the trajectory onto the layer's representing surface and use
    // that landing point for the surface-array bin lookup. at() returns
    // only the surfaces in the single bin, so it MUST be paired with
    // BoundaryTolerance::Infinite() above (the bin's primary sensor's
    // plane is treated as unbounded so out-of-bin crossings still register).
    // To re-enable: uncomment this block AND the "if (sensitive) Infinite()"
    // block in processSurface, AND comment out the neighbors() line below.
    // Vector3 lookupPosition = position;
    // if (SurfaceIntersection intersection =
    //         surfaceRepresentation()
    //             .intersect(gctx, position, direction)
    //             .closest();
    //     intersection.isValid()) {
    //   lookupPosition = intersection.position();
    // }
    // // get the candidates: at() returns only surfaces in the exact bin
    // // containing lookupPosition (~1 candidate per layer crossing).
    // const std::vector<const Surface*>& sensitiveSurfaces =
    //     m_surfaceArray->at(lookupPosition);
    // ── /PATCH ───────────────────────────────────────────────────────────
    //
    // Upstream-clean: neighbors() returns the bin and its surrounding
    // neighbors (3x3 grid), so trajectories crossing module seams or stave
    // overlaps are still found via finite-bounds intersect() rejection.
    //
    // ── PATCH: lookupPosition projection ────────────────────────────────
    // The default neighbors(position) uses bare current position for the
    // bin-lookup. That is mathematically equivalent to a RADIAL projection
    // of the current point onto the representing surface — fine for the
    // immediate next layer (propagator is right next to it), but wrong for
    // inner-layer return-arc crossings where the trajectory is far from
    // the layer in r and has azimuthal drift before arriving. We instead
    // ray-cast the trajectory's line onto the layer's representing surface
    // and use that crossing point for the bin lookup. To revert, comment
    // this PATCH block out (lookupPosition stays unused) and pass `position`
    // back to neighbors() below.
    Vector3 lookupPosition = position;
    // Bin lookup direction: in the helix regime, ray-cast the helix's
    // transverse circle onto the layer's representing cylinder so the
    // 3×3 phi-z bin neighbourhood centres on the *curved* trajectory's
    // actual layer crossing, not where a near-tangent line tangent
    // would land. Outside the helix regime, fall back to the line
    // ray-cast onto the representing surface.
    //
    // Same logic for shell mode: at the close-zone retarget the
    // trajectory is typically 2-3 mm of helix arc below silicon at
    // r_ideal, and over that short arc the phi drifts by ~arc/R (a few
    // mm at the cylinder). The bare current position would center the
    // bin lookup on the trajectory's CURRENT phi, missing the bin that
    // contains the sensors the helix will actually land on. Projecting
    // via helixBarrelCylinderLanding shifts the lookup to where the
    // trajectory lands on the layer rep cylinder, so neighbors() returns
    // the sensors phi-aligned with the actual landing.
    bool helixLookupOk = false;
    if (helixRegime &&
        surfaceRepresentation().type() == Surface::SurfaceType::Cylinder) {
      const auto& cylBounds = static_cast<const CylinderBounds&>(
          surfaceRepresentation().bounds());
      const double rLayer = cylBounds.get(CylinderBounds::eR);
      Vector3 helixLanding;
      if (detail::helixBarrelCylinderLanding(
              position, options.helixDirection, options.helixQOverP,
              options.helixBField, rLayer, helixLanding)) {
        lookupPosition = helixLanding;
        helixLookupOk = true;
      }
    }
    if (!helixLookupOk) {
      if (SurfaceIntersection sIntersection =
              surfaceRepresentation()
                  .intersect(gctx, position, direction)
                  .closest();
          sIntersection.isValid()) {
        lookupPosition = sIntersection.position();
      }
    }
    // ── /PATCH ──────────────────────────────────────────────────────────
    const std::vector<const Surface*>& sensitiveSurfaces =
        m_surfaceArray->neighbors(lookupPosition);
    // loop through and veto
    // - if the approach surface is the parameter surface
    // - if the surface is not compatible with the type(s) that are collected
    for (auto& sSurface : sensitiveSurfaces) {
      processSurface(*sSurface, true);
    }
  }

  // (C) representing surface section
  //
  // the layer surface itself is a testSurface.
  //
  // Behavior in shell mode:
  //   FAR phase (d > 3 from r_ideal): force-add layer rep as the synthetic
  //     "carrier" target. Bypasses acceptSurface: the layer rep is a
  //     passive geometric reference (no readout, no material) so the
  //     general filter correctly rejects it for normal navigation; in
  //     FAR phase we deliberately reuse it as a stop-target so the
  //     SteppingHelper FAR branch has something to step toward (d - 3
  //     along the curve). After the FAR descent finishes (d <= 3),
  //     SteppingHelper forces unreachable on this carrier so the navigator
  //     retargets and compatibleSurfaces re-runs in CLOSE phase.
  //   CLOSE phase: SKIP layer rep. Only sensor candidates are in
  //     navSurfaces (from section B), so the trajectory targets a real
  //     sensor next.
  if (!options.shellMode) {
    const Surface* layerSurface = &surfaceRepresentation();
    processSurface(*layerSurface);
  } else if (shellFar) {
    // Force-inject the layer rep, bypassing acceptSurface. Use a pure-
    // radial line predicate against the cylinder, with the sign chosen
    // by the trajectory's radial direction:
    //   outgoing (pr > 0): +r̂ — the cylinder is outward from us
    //   incoming (pr < 0): -r̂ — the cylinder is inward from us
    // Using the wrong sign would yield no forward landing on the layer-rep
    // (the closest intersection would be in the opposite direction, with
    // negative path, failing checkPathLength).
    // The path length stored here is informational only; the SteppingHelper
    // FAR branch recomputes d_to_r_ideal each step and overrides the step
    // size to (d - 3).
    const Surface* layerSurface = &surfaceRepresentation();
    const double rxyAtPos = std::sqrt(position[0] * position[0] +
                                       position[1] * position[1]);
    SurfaceIntersection sfi = SurfaceIntersection::invalid();
    if (rxyAtPos > 1.0e-6) {
      const double prDir =
          (position[0] * options.helixDirection[0] +
           position[1] * options.helixDirection[1]) / rxyAtPos;
      const double sign = (prDir > 0.0) ? +1.0 : -1.0;
      const Vector3 radialDir{sign * position[0] / rxyAtPos,
                              sign * position[1] / rxyAtPos, 0.0};
      sfi = layerSurface
                ->intersect(gctx, position, radialDir,
                            options.boundaryTolerance)
                .closest();
    } else {
      sfi = layerSurface
                ->intersect(gctx, position, direction,
                            options.boundaryTolerance)
                .closest();
    }
    if (sfi.isValid() &&
        detail::checkPathLength(sfi.pathLength(), nearLimit, farLimit)) {
      sIntersections.push_back(sfi);
    }
  }

  return sIntersections;
}

SurfaceIntersection Layer::surfaceOnApproach(
    const GeometryContext& gctx, const Vector3& position,
    const Vector3& direction, const NavigationOptions<Layer>& options) const {
  // resolve directive based by options
  // - options.resolvePassive is on -> always
  // - options.resolveSensitive is on -> always
  // - options.resolveMaterial is on
  //   && either sensitive or approach surfaces have material
  bool resolvePS = options.resolveSensitive || options.resolvePassive;
  bool resolveMS = options.resolveMaterial &&
                   (m_ssSensitiveSurfaces > 1 || m_ssApproachSurfaces > 1 ||
                    (surfaceRepresentation().surfaceMaterial() != nullptr));

  // The Limits
  double nearLimit = options.nearLimit;
  double farLimit = options.farLimit;

  // Helper function to find valid intersection
  auto findValidIntersection =
      [&](const SurfaceMultiIntersection& sfmi) -> SurfaceIntersection {
    for (const auto& sfi : sfmi.split()) {
      if (sfi.isValid() &&
          detail::checkPathLength(sfi.pathLength(), nearLimit, farLimit)) {
        return sfi;
      }
    }

    // Return an invalid one
    return SurfaceIntersection::invalid();
  };

  // Approach descriptor present and resolving is necessary
  if (m_approachDescriptor && (resolvePS || resolveMS)) {
    SurfaceIntersection aSurface = m_approachDescriptor->approachSurface(
        gctx, position, direction, options.boundaryTolerance, nearLimit,
        farLimit);
    return aSurface;
  }

  // Intersect and check the representing surface
  const Surface& rSurface = surfaceRepresentation();
  auto sIntersection =
      rSurface.intersect(gctx, position, direction, options.boundaryTolerance);
  return findValidIntersection(sIntersection);
}

}  // namespace Acts
