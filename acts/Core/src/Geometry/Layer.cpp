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
    // Predicate selection in compatibleSurfaces:
    //   sensitive modules in the helix regime → helix-plane (with raw
    //                                            tangent), bounded check on
    //                                            the actual curved landing.
    //   approach / boundary / layer-rep in the helix regime → line predicate
    //                                            with PURE radial inward as
    //                                            direction. These cylinders
    //                                            are radially well-conditioned
    //                                            in the barrel; pure radial
    //                                            gives the cleanest intersect.
    //   anything else (outward motion, non-barrel) → line with the supplied
    //                                            `direction` (raw tangent).
    // The path-length mismatch between helix-arc (sensors) and line-radial
    // (approaches) is resolved by the sensor-priority sort in
    // Navigator::resolveSurfaces: sensors come before approach surfaces in
    // navSurfaces regardless of path length, so the navigator targets every
    // sensor before falling through to the layer-exit approach.
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
    if (options.radialDownMode && sensitive && rxyAtPos > 1.0e-6) {
      // Radial-down mode for sensors: use line + pure radial-inward
      // direction (no helix). The bin lookup already used the radial-down
      // projection; the per-surface intersect now also uses the radial-
      // inward line for the bounded check. Sensors directly below the
      // trajectory's current (phi, z) pass; sensors offset in phi by more
      // than the silicon width fail. Strict BoundaryTolerance::None.
      const Vector3 radialDownDir{-position[0] / rxyAtPos,
                                  -position[1] / rxyAtPos, 0.0};
      sfi = sf.intersect(gctx, position, radialDownDir, boundaryTolerance)
                .closest();
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
  // SKIP this section when radial-down mode is active. apr=1 is always
  // reachable from anywhere in the shell with a small path length; once
  // it's in navSurfaces, the sensor-priority sort keeps sensors first
  // but the navigator falls through to apr=1 the moment all helix-bin
  // sensors get filtered out — short-circuiting the per-step
  // re-resolve that the radial-down mode is designed to drive. apr=1
  // re-enters the candidate pool via the standard layerTarget path once
  // radial-down mode exits at r < apr1_r + 3 mm.
  if (m_approachDescriptor && !options.radialDownMode &&
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
  // check the sensitive surfaces if you have some
  if (m_surfaceArray && (options.resolveMaterial || options.resolvePassive ||
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
    // RADIAL-DOWN MODE: project the current position radially onto the
    // layer's representing cylinder (r=r_m). The 3×3 bin neighbourhood
    // centres on the trajectory's CURRENT (phi, z) at r_m — which is
    // the bin the trajectory will pass through next as it descends. Re-
    // run every retarget as the trajectory steps, so each query samples
    // the local bin.
    if (options.radialDownMode &&
        surfaceRepresentation().type() == Surface::SurfaceType::Cylinder) {
      const auto& cylBounds = static_cast<const CylinderBounds&>(
          surfaceRepresentation().bounds());
      const double rLayer = cylBounds.get(CylinderBounds::eR);
      const double rxy = std::sqrt(position[0] * position[0] +
                                   position[1] * position[1]);
      if (rxy > 1.0e-6) {
        lookupPosition = Vector3{position[0] * (rLayer / rxy),
                                 position[1] * (rLayer / rxy),
                                 position[2]};
      }
    } else {
      // Bin lookup direction: in the helix regime, ray-cast the helix's
      // transverse circle onto the layer's representing cylinder so the
      // 3×3 phi-z bin neighbourhood centres on the *curved* trajectory's
      // actual layer crossing, not where a near-tangent line tangent
      // would land. Outside the helix regime, fall back to the line
      // ray-cast onto the representing surface as before.
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
  // the layer surface itself is a testSurface. SKIP in radial-down mode
  // for the same reason as the approach surfaces (section A): the layer
  // rep is always reachable from inside the shell and would short-circuit
  // the per-step sensor re-resolve.
  if (!options.radialDownMode) {
    const Surface* layerSurface = &surfaceRepresentation();
    processSurface(*layerSurface);
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
