// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/NavigatorError.hpp"
#include "Acts/Propagator/NavigatorOptions.hpp"
#include "Acts/Propagator/NavigatorStatistics.hpp"
#include "Acts/Propagator/StepLimitDiagnostics.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <string>

#include <boost/container/small_vector.hpp>

namespace Acts {

/// @brief The navigation options for the tracking geometry
///
/// @tparam object_t Type of the object for navigation to check against
template <typename object_t>
struct NavigationOptions {
  /// The boundary check directive
  BoundaryTolerance boundaryTolerance = BoundaryTolerance::None();

  // How to resolve the geometry
  /// Always look for sensitive
  bool resolveSensitive = true;
  /// Always look for material
  bool resolveMaterial = true;
  /// always look for passive
  bool resolvePassive = false;

  /// Hint for start object
  const object_t* startObject = nullptr;
  /// Hint for end object
  const object_t* endObject = nullptr;

  /// External surface identifier for which the boundary check is ignored
  std::vector<GeometryIdentifier> externalSurfaces = {};

  /// The minimum distance for a surface to be considered
  double nearLimit = 0;
  /// The maximum distance for a surface to be considered
  double farLimit = std::numeric_limits<double>::max();

  /// Helix model for the per-surface intersect predicate inside
  /// Layer::compatibleSurfaces. When `helixBField.norm() > 0` and
  /// `helixQOverP != 0`, the closed-form helix-plane intersect is used for
  /// sensitive surfaces (with bounded check on the actual curved landing
  /// point) instead of the line predicate. Outside this regime — and for
  /// non-sensitive surfaces — the line predicate is used unchanged.
  /// Set by the navigator from stepper state when it has reason to believe
  /// the trajectory's curvature significantly differs from a line over the
  /// candidate-list construction's relevant path lengths.
  double helixQOverP = 0.0;
  Vector3 helixBField = Vector3::Zero();
  /// The actual stepper tangent (raw direction). Distinct from the
  /// `direction` argument to `compatibleSurfaces`, which the navigator
  /// substitutes (e.g. with the radial bisector) to make the line-tangent
  /// predicate well-behaved at the bin lookup. The helix predicate needs
  /// the raw tangent to compute the helix center correctly.
  Vector3 helixDirection = Vector3::Zero();
};

/// @brief Steers the propagation through the geometry by providing the next
/// surface to be targeted.
///
/// The Navigator is part of the propagation and responsible for steering
/// the surface sequence to encounter all the relevant surfaces which are
/// intersected by the trajectory.
///
/// The current navigation stage is cached in the state struct and updated
/// when necessary. If any surface in the extrapolation flow is hit, it is
/// set to the navigation state, such that other actors can deal with it.
///
/// The current target surface is referenced by an index which points into
/// the navigation candidates. The navigation candidates are ordered by the
/// path length to the surface. If a surface is hit, the
/// `state.currentSurface` pointer is set. This actors to observe
/// that we are on a surface.
///
class Navigator {
 public:
  using NavigationSurfaces =
      boost::container::small_vector<SurfaceIntersection, 10>;

  using NavigationLayers =
      boost::container::small_vector<LayerIntersection, 10>;

  using NavigationBoundaries =
      boost::container::small_vector<BoundaryIntersection, 4>;

  using ExternalSurfaces = std::multimap<std::uint64_t, GeometryIdentifier>;

  using GeometryVersion = TrackingGeometry::GeometryVersion;

  /// The navigation stage
  enum struct Stage : int {
    initial = 0,
    surfaceTarget = 1,
    layerTarget = 2,
    boundaryTarget = 3,
  };

  /// The navigator configuration
  struct Config {
    /// Tracking Geometry for this Navigator
    std::shared_ptr<const TrackingGeometry> trackingGeometry{nullptr};

    /// stop at every sensitive surface (whether it has material or not)
    bool resolveSensitive = true;
    /// stop at every material surface (whether it is passive or not)
    bool resolveMaterial = true;
    /// stop at every surface regardless what it is
    bool resolvePassive = false;
  };

  /// The navigator options
  struct Options : public NavigatorPlainOptions {
    explicit Options(const GeometryContext& gctx)
        : NavigatorPlainOptions(gctx) {}

    /// The surface tolerance
    double surfaceTolerance = s_onSurfaceTolerance;

    /// The near limit to resolve surfaces
    double nearLimit = s_onSurfaceTolerance;

    /// The far limit to resolve surfaces
    double farLimit = std::numeric_limits<double>::max();

    /// Externally provided surfaces - these are tried to be hit
    ExternalSurfaces externalSurfaces = {};

    void insertExternalSurface(GeometryIdentifier geoid) {
      externalSurfaces.insert(
          std::pair<std::uint64_t, GeometryIdentifier>(geoid.layer(), geoid));
    }

    void setPlainOptions(const NavigatorPlainOptions& options) {
      static_cast<NavigatorPlainOptions&>(*this) = options;
    }
  };

  /// @brief Nested State struct
  ///
  /// It acts as an internal state which is created for every propagation and
  /// meant to keep thread-local navigation information.
  struct State {
    explicit State(const Options& options_) : options(options_) {}

    Options options;

    // Navigation on surface level
    /// the vector of navigation surfaces to work through
    NavigationSurfaces navSurfaces = {};
    /// the current surface index of the navigation state
    std::optional<std::size_t> navSurfaceIndex;

    // Navigation on layer level
    /// the vector of navigation layers to work through
    NavigationLayers navLayers = {};
    /// the current layer index of the navigation state
    std::optional<std::size_t> navLayerIndex;

    // Navigation on volume level
    /// the vector of boundary surfaces to work through
    NavigationBoundaries navBoundaries = {};
    /// the current boundary index of the navigation state
    std::optional<std::size_t> navBoundaryIndex;

    SurfaceIntersection& navSurface() {
      return navSurfaces.at(navSurfaceIndex.value());
    }
    LayerIntersection& navLayer() {
      return navLayers.at(navLayerIndex.value());
    }
    BoundaryIntersection& navBoundary() {
      return navBoundaries.at(navBoundaryIndex.value());
    }

    const TrackingVolume* startVolume = nullptr;
    const Layer* startLayer = nullptr;
    const Surface* startSurface = nullptr;
    const TrackingVolume* currentVolume = nullptr;
    const Layer* currentLayer = nullptr;
    const Surface* currentSurface = nullptr;
    const Surface* targetSurface = nullptr;

    bool navigationBreak = false;
    Stage navigationStage = Stage::initial;

    /// Flag indicating radial inward motion (toward beam axis)
    /// When true, use pure radial direction for layer/boundary resolution at turning points
    bool radiallyInward = false;

    /// Flag indicating if current volume is barrel (cylindrical layers)
    /// true = barrel region, false = endcap or unknown
    /// Only checked/updated when radiallyInward is true
    bool isInBarrelVolume = false;

    /// Stepper-side helix model snapshot, refreshed by the propagator each
    /// step. The navigator forwards these into NavigationOptions::helix*
    /// when calling Layer::compatibleSurfaces inside the inward-barrel
    /// regime, so the candidate-list construction can use the closed-form
    /// helix-plane intersect for sensitive surfaces. helixBField == 0 means
    /// "no helix info available, fall back to line predicate".
    double currentQOverP = 0.0;
    Vector3 currentBField = Vector3::Zero();

    NavigatorStatistics statistics;

    NavigationStream stream;

    void resetAfterLayerSwitch() {
      navSurfaces.clear();
      navSurfaceIndex.reset();
    }

    void resetAfterVolumeSwitch() {
      resetAfterLayerSwitch();

      navLayers.clear();
      navLayerIndex.reset();
      navBoundaries.clear();
      navBoundaryIndex.reset();

      currentLayer = nullptr;
    }

    void reset() {
      resetAfterVolumeSwitch();

      currentVolume = nullptr;
      currentSurface = nullptr;

      navigationBreak = false;
      navigationStage = Stage::initial;
    }
  };

  /// Constructor with configuration object
  ///
  /// @param cfg The navigator configuration
  /// @param _logger a logger instance
  explicit Navigator(Config cfg,
                     std::shared_ptr<const Logger> _logger =
                         getDefaultLogger("Navigator", Logging::Level::INFO))
      : m_cfg{std::move(cfg)}, m_logger{std::move(_logger)} {
    if (m_cfg.trackingGeometry == nullptr) {
      throw std::invalid_argument("Navigator: No tracking geometry provided.");
    }
    m_geometryVersion = m_cfg.trackingGeometry->geometryVersion();
  }

  State makeState(const Options& options) const {
    State state(options);
    return state;
  }

  const Surface* currentSurface(const State& state) const {
    return state.currentSurface;
  }

  const TrackingVolume* currentVolume(const State& state) const {
    return state.currentVolume;
  }

  const IVolumeMaterial* currentVolumeMaterial(const State& state) const {
    if (state.currentVolume == nullptr) {
      return nullptr;
    }
    return state.currentVolume->volumeMaterial();
  }

  const Surface* startSurface(const State& state) const {
    return state.startSurface;
  }

  const Surface* targetSurface(const State& state) const {
    return state.targetSurface;
  }

  bool endOfWorldReached(const State& state) const {
    return state.currentVolume == nullptr;
  }

  bool navigationBreak(const State& state) const {
    return state.navigationBreak;
  }

  /// @brief Initialize the navigator state
  ///
  /// This function initializes the navigator state for a new propagation.
  ///
  /// @param state The navigation state
  /// @param position The start position
  /// @param direction The start direction
  /// @param propagationDirection The propagation direction
  ///
  /// @return Indication if the initialization was successful
  [[nodiscard]] Result<void> initialize(State& state, const Vector3& position,
                                        const Vector3& direction,
                                        Direction propagationDirection) const {
    (void)propagationDirection;

    ACTS_VERBOSE(volInfo(state) << "Initialization.");

    auto printGeometryVersion = [](auto ver) {
      using enum TrackingGeometry::GeometryVersion;
      switch (ver) {
        case Gen1:
          return "Gen1";
        case Gen3:
          return "Gen3";
        default:
          throw std::runtime_error("Unknown geometry version.");
      }
    };
    ACTS_VERBOSE(volInfo(state) << "Geometry version is: "
                                << printGeometryVersion(m_geometryVersion));

    state.reset();

    // Empirical pre-allocation of candidates for the next navigation iteration.
    // @TODO: Make this user configurable through the configuration
    state.stream.candidates().reserve(50);

    state.startSurface = state.options.startSurface;
    state.targetSurface = state.options.targetSurface;

    // @TODO: Implement fast initialization with Gen3. This requires the volume lookup to work properly

    // Fast Navigation initialization for start condition:
    // - short-cut through object association, saves navigation in the
    // - geometry and volume tree search for the lowest volume
    if (state.startSurface != nullptr &&
        state.startSurface->associatedLayer() != nullptr) {
      ACTS_VERBOSE(
          volInfo(state)
          << "Fast start initialization through association from Surface.");

      state.startLayer = state.startSurface->associatedLayer();
      state.startVolume = state.startLayer->trackingVolume();
    } else if (state.startVolume != nullptr) {
      ACTS_VERBOSE(
          volInfo(state)
          << "Fast start initialization through association from Volume.");

      state.startLayer = state.startVolume->associatedLayer(
          state.options.geoContext, position);
    } else {
      ACTS_VERBOSE(volInfo(state)
                   << "Slow start initialization through search.");
      ACTS_VERBOSE(volInfo(state)
                   << "Starting from position " << toString(position)
                   << " and direction " << toString(direction));

      // current volume and layer search through global search
      state.startVolume = m_cfg.trackingGeometry->lowestTrackingVolume(
          state.options.geoContext, position);

      if (state.startVolume != nullptr) {
        state.startLayer = state.startVolume->associatedLayer(
            state.options.geoContext, position);
      } else {
        ACTS_ERROR(volInfo(state)
                   << "No start volume resolved. Nothing left to do.");
        state.navigationBreak = true;
      }
    }

    state.currentVolume = state.startVolume;
    state.currentLayer = state.startLayer;
    state.currentSurface = state.startSurface;

    if (state.currentVolume != nullptr) {
      ACTS_VERBOSE(volInfo(state) << "Start volume resolved "
                                  << state.currentVolume->geometryId());

      if (!state.currentVolume->inside(position,
                                       state.options.surfaceTolerance)) {
        ACTS_DEBUG(
            volInfo(state)
            << "We did not end up inside the expected volume. position = "
            << position.transpose());

        return Result<void>::failure(NavigatorError::NotInsideExpectedVolume);
      }
    }
    if (state.currentLayer != nullptr) {
      ACTS_VERBOSE(volInfo(state) << "Start layer resolved "
                                  << state.currentLayer->geometryId());
    }
    if (state.currentSurface != nullptr) {
      ACTS_VERBOSE(volInfo(state) << "Start surface resolved "
                                  << state.currentSurface->geometryId());

      if (!state.currentSurface->isOnSurface(
              state.options.geoContext, position, direction,
              BoundaryTolerance::Infinite(), state.options.surfaceTolerance)) {
        ACTS_DEBUG(volInfo(state)
                   << "We did not end up on the expected surface. surface = "
                   << state.currentSurface->geometryId()
                   << " position = " << position.transpose()
                   << " direction = " << direction.transpose());

        return Result<void>::failure(NavigatorError::NotOnExpectedSurface);
      }
    }

    return Result<void>::success();
  }

  /// @brief Get the next target surface
  ///
  /// This function gets the next target surface for the propagation.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return The next target surface
  NavigationTarget nextTarget(State& state, const Vector3& position,
                              const Vector3& direction) const {
    // Reset the current surface
    state.currentSurface = nullptr;

    if (inactive(state)) {
      return NavigationTarget::None();
    }

    ACTS_VERBOSE(volInfo(state) << "Entering Navigator::nextTarget.");

    auto tryGetNextTarget = [&]() -> NavigationTarget {
      // Try targeting the surfaces - then layers - then boundaries

      if (state.navigationStage == Stage::initial) {
        ACTS_VERBOSE(volInfo(state) << "Target surfaces.");
        state.navigationStage = Stage::surfaceTarget;
      }

      if (state.navigationStage == Stage::surfaceTarget) {
        if (!state.navSurfaceIndex.has_value()) {
          // First time (or after a retarget below), resolve the surfaces
          resolveSurfaces(state, position, direction);
          state.navSurfaceIndex = 0;
          // Log surface list contents
          ACTS_VERBOSE(volInfo(state) << "navSurfaces list has " << state.navSurfaces.size() << " items:");
          for (std::size_t i = 0; i < state.navSurfaces.size(); ++i) {
            ACTS_VERBOSE("  [" << i << "] geoId=" << state.navSurfaces[i].surface().geometryId()
                         << " pathLength=" << state.navSurfaces[i].pathLength());
          }
        } else {
          ++state.navSurfaceIndex.value();
        }
        if (state.navSurfaceIndex.value() < state.navSurfaces.size()) {
          ACTS_VERBOSE(volInfo(state) << "Target set to next surface.");
          return NavigationTarget(state.navSurface().surface(),
                                  state.navSurface().index(),
                                  state.navSurface().boundaryTolerance());
        }

        // No (more) candidates — switch to layers
        ACTS_VERBOSE(volInfo(state) << "Target layers.");
        if (m_geometryVersion == GeometryVersion::Gen1) {
          state.navigationStage = Stage::layerTarget;
        } else {
          state.navigationStage = Stage::boundaryTarget;
        }
      }

      if (state.navigationStage == Stage::layerTarget) {
        if (!state.navLayerIndex.has_value()) {
          // First time, resolve the layers
          resolveLayers(state, position, direction);
          state.navLayerIndex = 0;
          // Log layer list contents
          ACTS_VERBOSE(volInfo(state) << "navLayers list has " << state.navLayers.size() << " items:");
          for (std::size_t i = 0; i < state.navLayers.size(); ++i) {
            ACTS_VERBOSE("  [" << i << "] geoId=" << state.navLayers[i].first.surface().geometryId()
                         << " pathLength=" << state.navLayers[i].first.pathLength());
          }
        } else {
          ++state.navLayerIndex.value();
        }
        if (state.navLayerIndex.value() < state.navLayers.size()) {
          ACTS_VERBOSE(volInfo(state) << "Target set to next layer.");
          return NavigationTarget(state.navLayer().first.surface(),
                                  state.navLayer().first.index(),
                                  state.navLayer().first.boundaryTolerance());
        } else {
          // This was the last layer, switch to boundaries
          ACTS_VERBOSE(volInfo(state) << "Target boundaries.");
          state.navigationStage = Stage::boundaryTarget;
        }
      }

      if (state.navigationStage == Stage::boundaryTarget) {
        if (!state.navBoundaryIndex.has_value()) {
          // First time, resolve the boundaries
          resolveBoundaries(state, position, direction);
          state.navBoundaryIndex = 0;
          // Log boundary list contents
          ACTS_VERBOSE(volInfo(state) << "navBoundaries list has " << state.navBoundaries.size() << " items:");
          for (std::size_t i = 0; i < state.navBoundaries.size(); ++i) {
            ACTS_VERBOSE("  [" << i << "] geoId=" << state.navBoundaries[i].intersection.surface().geometryId()
                         << " pathLength=" << state.navBoundaries[i].intersection.pathLength());
          }
        } else {
          ++state.navBoundaryIndex.value();
        }
        if (state.navBoundaryIndex.value() < state.navBoundaries.size()) {
          ACTS_VERBOSE(volInfo(state) << "Target set to next boundary.");
          return NavigationTarget(
              state.navBoundary().intersection.surface(),
              state.navBoundary().intersection.index(),
              state.navBoundary().intersection.boundaryTolerance());
        } else {
          // This was the last boundary, we have to leave the volume somehow,
          // renavigate
          ACTS_VERBOSE(volInfo(state)
                       << "Boundary targets exhausted. Renavigate.");
        }
      }

      ACTS_VERBOSE(volInfo(state)
                   << "Unknown state. No target found. Renavigate.");
      return NavigationTarget::None();
    };

    NavigationTarget nextTarget = tryGetNextTarget();
    if (!nextTarget.isNone()) {
      return nextTarget;
    }

    state.reset();
    ++state.statistics.nRenavigations;

    // We might have punched through a boundary and entered another volume
    // so we have to reinitialize
    state.currentVolume = m_cfg.trackingGeometry->lowestTrackingVolume(
        state.options.geoContext, position);

    if (state.currentVolume == nullptr) {
      ACTS_VERBOSE(volInfo(state) << "No volume found, stop navigation.");
      state.navigationBreak = true;
      return NavigationTarget::None();
    }

    state.currentLayer = state.currentVolume->associatedLayer(
        state.options.geoContext, position);

    ACTS_VERBOSE(volInfo(state) << "Resolved volume and layer.");

    // Rerun the targeting
    nextTarget = tryGetNextTarget();
    if (!nextTarget.isNone()) {
      return nextTarget;
    }

    ACTS_VERBOSE(
        volInfo(state)
        << "No targets found again, we got really lost! Stop navigation.");
    state.navigationBreak = true;
    return NavigationTarget::None();
  }

  /// @brief Check if the current target is still valid
  ///
  /// This function checks if the target is valid.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  ///
  /// @return True if the target is valid
  bool checkTargetValid(const State& state, const Vector3& position,
                        const Vector3& direction) const {
    (void)position;
    (void)direction;

    return state.navigationStage != Stage::initial;
  }

  /// @brief Handle the surface reached
  ///
  /// This function handles the surface reached.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  /// @param surface The surface reached
  void handleSurfaceReached(State& state, const Vector3& position,
                            const Vector3& direction,
                            const Surface& surface) const {
    if (inactive(state)) {
      return;
    }

    ACTS_VERBOSE(volInfo(state) << "Entering Navigator::handleSurfaceReached.");

    state.currentSurface = &surface;

    ACTS_VERBOSE(volInfo(state)
                 << "Current surface: " << state.currentSurface->geometryId());

    if (state.navigationStage == Stage::surfaceTarget &&
        &state.navSurface().surface() == &surface) {
      ACTS_VERBOSE(volInfo(state) << "Handling surface status.");
      return;
    }

    if (state.navigationStage == Stage::layerTarget &&
        &state.navLayer().first.surface() == &surface) {
      ACTS_VERBOSE(volInfo(state) << "Handling layer status.");

      // Switch to the next layer
      state.currentLayer = state.navLayer().second;
      state.navigationStage = Stage::surfaceTarget;

      // partial reset
      state.resetAfterLayerSwitch();

      return;
    }

    if (state.navigationStage == Stage::boundaryTarget &&
        &state.navBoundary().intersection.surface() == &surface) {
      ACTS_VERBOSE(volInfo(state) << "Handling boundary status.");

      if (m_geometryVersion == GeometryVersion::Gen1) {
        // Switch to the next volume using the boundary
        const BoundarySurface* boundary = state.navBoundary().boundarySurface;
        assert(boundary != nullptr && "Retrieved boundary surface is nullptr");
        state.currentVolume = boundary->attachedVolume(state.options.geoContext,
                                                       position, direction);
      } else {
        const Portal* portal = state.navBoundary().portal;
        assert(portal != nullptr && "Retrieved portal is nullptr");
        auto res = portal->resolveVolume(state.options.geoContext, position,
                                         direction);
        if (!res.ok()) {
          ACTS_ERROR(volInfo(state)
                     << "Failed to resolve volume through portal: "
                     << res.error().message());
          return;
        }

        state.currentVolume = res.value();
      }

      // partial reset
      state.resetAfterVolumeSwitch();

      if (state.currentVolume != nullptr) {
        ACTS_VERBOSE(volInfo(state) << "Volume updated.");
        if (m_geometryVersion == GeometryVersion::Gen1) {
          state.navigationStage = Stage::layerTarget;
        } else {
          state.navigationStage = Stage::surfaceTarget;
        }
      } else {
        ACTS_VERBOSE(volInfo(state)
                     << "No more volume to progress to, stopping navigation.");
        state.navigationBreak = true;
      }

      return;
    }

    ACTS_ERROR(volInfo(state) << "Surface reached but unknown state.");
  }

 private:

  /// @brief Resolve compatible surfaces
  ///
  /// This function resolves the compatible surfaces for the navigation.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void resolveSurfaces(State& state, const Vector3& position,
                       const Vector3& direction) const {
    ACTS_VERBOSE(volInfo(state) << "Searching for compatible surfaces.");

    if (m_geometryVersion == GeometryVersion::Gen1) {
      const Layer* currentLayer = state.currentLayer;

      if (currentLayer == nullptr) {
        ACTS_VERBOSE(volInfo(state) << "No layer to resolve surfaces.");
        return;
      }

      const Surface* layerSurface = &currentLayer->surfaceRepresentation();

      NavigationOptions<Surface> navOpts;
      navOpts.resolveSensitive = m_cfg.resolveSensitive;
      navOpts.resolveMaterial = m_cfg.resolveMaterial;
      navOpts.resolvePassive = m_cfg.resolvePassive;
      navOpts.startObject = state.currentSurface;
      navOpts.endObject = state.targetSurface;
      navOpts.nearLimit = state.options.nearLimit;
      navOpts.farLimit = state.options.farLimit;
      // Forward helix model into compatibleSurfaces for the inward-barrel
      // regime so the candidate-list construction's per-surface intersect
      // uses the closed-form helix-plane predicate (with bounded check on
      // the actual curved landing) instead of the line tangent. The raw
      // tangent (`direction`) is passed via helixDirection because the
      // `direction` argument to compatibleSurfaces is the bisector below.
      // Set ACTS_DISABLE_HELIX_INTERSECT=1 in the environment for a clean
      // line-everywhere baseline (matches the SteppingHelper switch).
      static const bool helixDisabled =
          std::getenv("ACTS_DISABLE_HELIX_INTERSECT") != nullptr;
      // Forward helix info for compatibleSurfaces unconditionally (unless
      // the env override disables helix). The closed-form helix-plane
      // intersect is direction-agnostic and self-validating: it checks
      // axial-B and barrel-orientation preconditions internally and
      // falls back if they fail.
      if (!helixDisabled) {
        navOpts.helixQOverP = state.currentQOverP;
        navOpts.helixBField = state.currentBField;
        navOpts.helixDirection = direction;
      }

      if (!state.options.externalSurfaces.empty()) {
        auto layerId = layerSurface->geometryId().layer();
        auto externalSurfaceRange =
            state.options.externalSurfaces.equal_range(layerId);
        navOpts.externalSurfaces.reserve(
            state.options.externalSurfaces.count(layerId));
        for (auto itSurface = externalSurfaceRange.first;
             itSurface != externalSurfaceRange.second; itSurface++) {
          navOpts.externalSurfaces.push_back(itSurface->second);
        }
      }
      // Pass the RAW stepper direction (no bisector substitution). The
      // bisector heuristic was masking inconsistent path-length predictions
      // between cylindrical approach surfaces (under-estimated) and planar
      // sensors (helix-correct), causing the navigator to target approach
      // surfaces ahead of the sensors. Layer::compatibleSurfaces handles its
      // own substitutions internally:
      //   sensors in helix regime → helix-plane (with raw tangent),
      //   approach/boundary in helix regime → line + pure-radial-inward,
      //   anything else → line + raw tangent.
      state.navSurfaces = currentLayer->compatibleSurfaces(
          state.options.geoContext, position, direction, navOpts);
      // Sort the surfaces. Order of preference (high to low):
      //   1. external surfaces (infinite boundary tolerance) — always first.
      //   2. sensitive sensor modules — before approach / boundary / layer
      //      surfaces, regardless of path length. The path-length predicates
      //      for sensors (helix-plane arc) and approaches (line+pure-radial)
      //      use different conventions, so direct path-length ordering is
      //      unreliable; sensors-first ensures the navigator visits every
      //      sensor candidate before falling through to a layer-exit
      //      approach surface and skipping a sensor by mistake.
      //   3. by path length within the same priority class.
      std::ranges::sort(
          state.navSurfaces,
          [&state](const SurfaceIntersection& a, const SurfaceIntersection& b) {
            const bool aIsExternal = a.boundaryTolerance().isInfinite();
            const bool bIsExternal = b.boundaryTolerance().isInfinite();
            if (aIsExternal != bIsExternal) {
              return aIsExternal;  // external first
            }
            const bool aIsSensitive = a.surface().geometryId().sensitive() != 0;
            const bool bIsSensitive = b.surface().geometryId().sensitive() != 0;
            if (aIsSensitive != bIsSensitive) {
              return aIsSensitive;  // sensors before non-sensors
            }
            // Same priority class — sort by path length, breaking ties by
            // geometry id.
            if (std::abs(a.pathLength() - b.pathLength()) >
                state.options.surfaceTolerance) {
              return SurfaceIntersection::pathLengthOrder(a, b);
            }
            return a.surface().geometryId() < b.surface().geometryId();
          });
      // For now we implicitly remove overlapping surfaces.
      // For track finding it might be useful to discover overlapping surfaces
      // and check for compatible measurements. This is under investigation
      // and might be implemented in the future.
      auto toBeRemoved = std::ranges::unique(
          state.navSurfaces, [&](const auto& a, const auto& b) {
            return std::abs(a.pathLength() - b.pathLength()) <
                   state.options.surfaceTolerance;
          });
      if (toBeRemoved.begin() != toBeRemoved.end()) {
        ACTS_VERBOSE(volInfo(state)
                     << "Removing "
                     << std::distance(toBeRemoved.begin(), toBeRemoved.end())
                     << " overlapping surfaces.");
      }
      state.navSurfaces.erase(toBeRemoved.begin(), toBeRemoved.end());

      
    } else {
      // @TODO: What to do with external surfaces?
      // Gen 3 !
      state.stream.reset();
      AppendOnlyNavigationStream appendOnly{state.stream};
      NavigationArguments args;
      args.position = position;
      args.direction = direction;
      args.wantsPortals = false;
      args.wantsSurfaces = true;
      state.currentVolume->initializeNavigationCandidates(args, appendOnly,
                                                          logger());

      // Filter out portals before intersection

      ACTS_VERBOSE(volInfo(state)
                   << "Found " << state.stream.candidates().size()
                   << " navigation candidates.");

      state.stream.initialize(state.options.geoContext, {position, direction},
                              BoundaryTolerance::None(),
                              state.options.surfaceTolerance);
      ACTS_VERBOSE(volInfo(state)
                   << "Now " << state.stream.candidates().size()
                   << " navigation candidates after initialization");

      state.navSurfaces.clear();

      auto it = std::ranges::find_if(
          state.stream.candidates(), [&](const auto& candidate) {
            return detail::checkPathLength(candidate.intersection.pathLength(),
                                           state.options.nearLimit,
                                           state.options.farLimit, logger());
          });

      for (; it != state.stream.candidates().end(); ++it) {
        state.navSurfaces.emplace_back(it->intersection);
      }
    }

    // Print surface information
    if (logger().doPrint(Logging::VERBOSE)) {
      std::ostringstream os;
      os << state.navSurfaces.size();
      os << " surface candidates found at path(s): ";
      for (auto& sfc : state.navSurfaces) {
        os << sfc.pathLength() << "  ";
      }
      logger().log(Logging::VERBOSE, os.str());

      // Detailed per-candidate dump (geometry id, layer/approach/sensitive,
      // path length, hasMaterial, index). This helps debugging missing
      // return-leg layers.
      std::ostringstream ods;
      ods << "Detailed surface candidates: ";
      for (auto& sfc : state.navSurfaces) {
        const auto& surf = sfc.surface();
        ods << "geo=" << surf.geometryId() << " ";
        // try to detect if this is an approach/sensitive/representing
        ods << "path=" << sfc.pathLength() << " ";
        ods << "idx=" << sfc.index() << " ";
        ods << "mat=" << (surf.surfaceMaterial() != nullptr) << " ";
        ods << "| ";
      }
      logger().log(Logging::VERBOSE, ods.str());
    }

    if (state.navSurfaces.empty()) {
      ACTS_VERBOSE(volInfo(state) << "No surface candidates found.");
    }
  }

  /// @brief Resolve compatible layers
  ///
  /// This function resolves the compatible layers for the navigation.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void resolveLayers(State& state, const Vector3& position,
                     const Vector3& direction) const {
    ACTS_VERBOSE(volInfo(state) << "Searching for compatible layers.");

   
    NavigationOptions<Layer> navOpts;
    navOpts.resolveSensitive = m_cfg.resolveSensitive;
    navOpts.resolveMaterial = m_cfg.resolveMaterial;
    navOpts.resolvePassive = m_cfg.resolvePassive;
    navOpts.startObject = state.currentLayer;
    navOpts.nearLimit = state.options.nearLimit;
    navOpts.farLimit = state.options.farLimit;
    
    // Use pure radial direction for layer resolution when going inward in barrel regions
    Vector3 effectiveDirection = computeEffectiveDirection(state, position, direction);
    

    // Debug: Check what associatedLayer returns for this position
    const Layer* associatedLayerResult = state.currentVolume->associatedLayer(state.options.geoContext, position);
    ACTS_VERBOSE(volInfo(state) << "associatedLayer() returned: " << (associatedLayerResult ? associatedLayerResult->geometryId() : GeometryIdentifier{}));

    // Request the compatible layers
    state.navLayers = state.currentVolume->compatibleLayers(
        state.options.geoContext, position, effectiveDirection, navOpts);
    std::ranges::sort(state.navLayers, [](const auto& a, const auto& b) {
      return SurfaceIntersection::pathLengthOrder(a.first, b.first);
    });

    // Print layer information
    if (logger().doPrint(Logging::VERBOSE)) {
      std::ostringstream os;
      os << state.navLayers.size();
      os << " layer candidates found at path(s): ";
      for (auto& lc : state.navLayers) {
        os << lc.first.pathLength() << "  ";
      }
      logger().log(Logging::VERBOSE, os.str());
    }

    if (state.navLayers.empty()) {
      ACTS_VERBOSE(volInfo(state) << "No layer candidates found.");
    }
  }

  /// @brief Resolve compatible boundaries
  ///
  /// This function resolves the compatible boundaries for the navigation.
  ///
  /// @param state The navigation state
  /// @param position The current position
  /// @param direction The current direction
  void resolveBoundaries(State& state, const Vector3& position,
                         const Vector3& direction) const {
    ACTS_VERBOSE(volInfo(state) << "Searching for compatible boundaries.");

    NavigationOptions<Surface> navOpts;
    navOpts.startObject = state.currentSurface;
    navOpts.nearLimit = state.options.nearLimit;
    navOpts.farLimit = state.options.farLimit;

    ACTS_VERBOSE(volInfo(state)
                 << "Try to find boundaries, we are at: " << toString(position)
                 << ", dir: " << toString(direction));

    // Use pure radial direction for boundary resolution when going inward in barrel regions
    Vector3 effectiveDirection = computeEffectiveDirection(state, position, direction);

    if (m_geometryVersion == GeometryVersion::Gen1) {
      // Request the compatible boundaries
      state.navBoundaries = state.currentVolume->compatibleBoundaries(
          state.options.geoContext, position, effectiveDirection, navOpts, logger());
      std::ranges::sort(state.navBoundaries, [](const auto& a, const auto& b) {
        return SurfaceIntersection::pathLengthOrder(a.intersection,
                                                    b.intersection);
      });
        // ...existing code...
    } else {
      // Gen 3 !
      state.stream.reset();
      AppendOnlyNavigationStream appendOnly{state.stream};
      NavigationArguments args;
      args.position = position;
      args.direction = effectiveDirection;
      args.wantsPortals = true;
      args.wantsSurfaces = false;
      state.currentVolume->initializeNavigationCandidates(args, appendOnly,
                                                          logger());

      ACTS_VERBOSE(volInfo(state)
                   << "Found " << state.stream.candidates().size()
                   << " navigation candidates.");

      state.stream.initialize(state.options.geoContext, {position, effectiveDirection},
                              BoundaryTolerance::None(),
                              state.options.surfaceTolerance);

      state.navBoundaries.clear();
      for (auto& candidate : state.stream.candidates()) {
        if (!detail::checkPathLength(candidate.intersection.pathLength(),
                                     state.options.nearLimit,
                                     state.options.farLimit, logger())) {
          continue;
        }

        state.navBoundaries.emplace_back(candidate.intersection, nullptr,
                                         candidate.portal);
      }
    }

    // Print boundary information
    // if (logger().doPrint(Logging::VERBOSE)) {
    //   std::ostringstream os;
    //   os << state.navBoundaries.size();
    //   os << " boundary candidates found at path(s): ";
    //   for (auto& bc : state.navBoundaries) {
    //     os << bc.intersection.pathLength() << "  ";
    //   }
    //   logger().log(Logging::VERBOSE, os.str());

    //   // Detailed per-candidate dump (geometry id, type, path length, index)
    //   std::ostringstream ods;
    //   ods << "Detailed boundary candidates: ";
    //   for (auto& bc : state.navBoundaries) {
    //     const auto& surf = bc.intersection.surface();
    //     std::string typeStr;
    //     switch (surf.type()) {
    //       case Surface::SurfaceType::Cylinder: typeStr = "Cylinder"; break;
    //       case Surface::SurfaceType::Disc: typeStr = "Disc"; break;
    //       case Surface::SurfaceType::Plane: typeStr = "Plane"; break;
    //       case Surface::SurfaceType::Cone: typeStr = "Cone"; break;
    //       default: typeStr = "Unknown"; break;
    //     }
    //     // intersection position
    //     auto ipos = bc.intersection.position();
    //     // surface center
    //     auto center = surf.center(state.options.geoContext);
    //     // surface normal at intersection (use effectiveDirection for sign)
    //     Vector3 normalVec = Vector3::Zero();
    //     try {
    //       normalVec = surf.normal(state.options.geoContext, ipos, effectiveDirection);
    //     } catch (...) {
    //       normalVec = Vector3::Zero();
    //     }

    //     ods << "geo=" << surf.geometryId() << " ";
    //     ods << "type=" << typeStr << " ";
    //     ods << "path=" << bc.intersection.pathLength() << " ";
    //     ods << "idx=" << bc.intersection.index() << " ";
    //     ods << "ipos=" << toString(ipos) << " ";
    //     ods << "center=" << toString(center) << " ";
    //     // Print shape-specific sizing info: cylinder radius/zBounds or disc rmin/rmax
    //     auto bvals = surf.bounds().values();
    //     if (surf.type() == Surface::SurfaceType::Cylinder && bvals.size() >= 2) {
    //       ods << "radius=" << bvals[0] << " ";
    //       ods << "halfLengthZ=" << bvals[1] << " ";
    //     } else if ((surf.type() == Surface::SurfaceType::Disc) && bvals.size() >= 2) {
    //       ods << "rMin=" << bvals[0] << " rMax=" << bvals[1] << " ";
    //     }
    //     ods << "normal=" << toString(normalVec) << " | ";
    //   }
    //   logger().log(Logging::VERBOSE, ods.str());
    // }

    // if (state.navBoundaries.empty()) {
    //   ACTS_VERBOSE(volInfo(state) << "No boundary candidates found.");
    // }
  }

  /// @brief Check if the navigator is inactive
  ///
  /// This function checks if the navigator is inactive.
  ///
  /// @param state The navigation state
  ///
  /// @return True if the navigator is inactive
  bool inactive(const State& state) const {
    // Turn the navigator into void when you are instructed to do nothing
    if (!m_cfg.resolveSensitive && !m_cfg.resolveMaterial &&
        !m_cfg.resolvePassive) {
      return true;
    }

    if (state.navigationBreak) {
      return true;
    }

    return false;
  }

 private:
  template <typename propagator_state_t>
  std::string volInfo(const propagator_state_t& state) const {
    return (state.currentVolume != nullptr ? state.currentVolume->volumeName()
                                           : "No Volume") +
           " | ";
  }

  /// @brief Check if volume contains barrel-like (cylindrical) layers
  /// 
  /// @param volume The tracking volume to check
  /// @return true if volume has cylindrical layers (barrel), false otherwise
  bool isBarrelVolume(const TrackingVolume* volume) const {
    if (volume == nullptr || volume->confinedLayers() == nullptr) {
      return false;
    }
    
    const auto& layers = volume->confinedLayers()->arrayObjects();
    if (layers.empty()) {
      return false;
    }
    
    // Check first layer's surface type - all layers in a volume are the same type
    auto surfType = layers.front()->surfaceRepresentation().type();
    return (surfType == Surface::SurfaceType::Cylinder);
  }

  /// @brief Compute effective direction for navigation at turning points
  ///
  /// When particles turn inward in barrel regions, the raw momentum direction
  /// is near-tangential at the apex and the line-plane intersect mis-targets
  /// inner layers and modules. Two substitutions are needed:
  ///
  ///  - Pure radial (`useBisector = false`): for compatibleLayers /
  ///    compatibleBoundaries volume-wide candidate searches. A bisector
  ///    here is mostly tangential at the turning point and drops inner
  ///    layers from the candidate list.
  ///  - Bisector (`useBisector = true`): for compatibleSurfaces on the
  ///    current layer (bin lookup + per-surface intersect). The bisector
  ///    keeps a well-conditioned radial component for the line-plane
  ///    denominator while preserving forward motion, so the bin lookup
  ///    lands on the right phi-bin.
  ///
  /// @param state The navigation state
  /// @param position Current position
  /// @param direction Original propagation direction
  /// @param useBisector If true, return bisector(direction, -r̂); else -r̂.
  /// @return Effective direction for intersection calculations
  Vector3 computeEffectiveDirection(State& state, const Vector3& position,
                                    const Vector3& direction,
                                    bool useBisector = false) const {
    Vector3 effectiveDirection = direction;

    if (!state.radiallyInward) {
      return effectiveDirection;
    }

    state.isInBarrelVolume = isBarrelVolume(state.currentVolume);
    if (!state.isInBarrelVolume) {
      ACTS_VERBOSE(volInfo(state) << "Endcap volume - keeping original direction with z-component");
      return effectiveDirection;
    }

    double r_xy = std::sqrt(position[0] * position[0] + position[1] * position[1]);
    if (r_xy <= 1e-6) {
      return effectiveDirection;
    }

    double r_hat_x = position[0] / r_xy;
    double r_hat_y = position[1] / r_xy;

    if (useBisector) {
      double sum_x = direction[0] - r_hat_x;
      double sum_y = direction[1] - r_hat_y;
      double sum_z = direction[2];
      double sum_mag = std::sqrt(sum_x * sum_x + sum_y * sum_y + sum_z * sum_z);
      if (sum_mag > 1e-9) {
        effectiveDirection[0] = sum_x / sum_mag;
        effectiveDirection[1] = sum_y / sum_mag;
        effectiveDirection[2] = sum_z / sum_mag;
        ACTS_VERBOSE("Using bisector(tangent, radial-inward) effective direction in barrel");
        return effectiveDirection;
      }
      ACTS_VERBOSE("Bisector degenerate — falling back to pure radial inward");
    }

    effectiveDirection[0] = -r_hat_x;
    effectiveDirection[1] = -r_hat_y;
    effectiveDirection[2] = 0.0;
    ACTS_VERBOSE("Using pure radial inward direction in barrel");
    return effectiveDirection;
  }

  const Logger& logger() const { return *m_logger; }

  Config m_cfg;

  // Cached so we don't have to query the TrackingGeometry constantly.
  TrackingGeometry::GeometryVersion m_geometryVersion;

  std::shared_ptr<const Logger> m_logger;
};

}  // namespace Acts
