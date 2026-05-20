// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/Propagator.hpp"

#include "Acts/EventData/TrackParametersConcept.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"

#include <cstdlib>
#include <sstream>
#include "Acts/Propagator/ActorList.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/PropagatorError.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/StepLimitDiagnostics.hpp"
#include "Acts/Propagator/detail/LoopProtection.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <concepts>

template <typename S, typename N>
template <typename propagator_state_t>
Acts::Result<void> Acts::Propagator<S, N>::propagate(
    propagator_state_t& state) const {
  ACTS_VERBOSE("Entering propagation.");

  state.stage = PropagatorStage::prePropagation;

  // Pre-Propagation: call to the actor list, abort condition check
  state.options.actorList.act(state, m_stepper, m_navigator, logger());

  if (state.options.actorList.checkAbort(state, m_stepper, m_navigator,
                                         logger())) {
    ACTS_VERBOSE("Propagation terminated without going into stepping loop.");

    state.stage = PropagatorStage::postPropagation;

    state.options.actorList.act(state, m_stepper, m_navigator, logger());

    return Result<void>::success();
  }

  auto getNextTarget = [&]() -> Result<NavigationTarget> {
    for (unsigned int i = 0; i < state.options.maxTargetSkipping; ++i) {
      NavigationTarget nextTarget = m_navigator.nextTarget(
          state.navigation, state.position, state.direction);
      if (nextTarget.isNone()) {
        return NavigationTarget::None();
      }
      IntersectionStatus preStepSurfaceStatus = m_stepper.updateSurfaceStatus(
          state.stepping, *nextTarget.surface,
          nextTarget.surfaceIntersectionIndex, state.options.direction,
          nextTarget.boundaryTolerance, state.options.surfaceTolerance,
          ConstrainedStep::Type::Navigator, state.navigation.isInBarrelVolume,
          logger());
      if (preStepSurfaceStatus == IntersectionStatus::onSurface) {
        // This indicates a geometry overlap which is not handled by the
        // navigator, so we skip this target.
        // This can also happen in a well-behaved geometry with external
        // surfaces.
        ACTS_VERBOSE("Pre-step surface status is onSurface, skipping target "
                     << nextTarget.surface->geometryId());
        continue;
      }
      if (preStepSurfaceStatus == IntersectionStatus::reachable) {
        return nextTarget;
      }
    }

    ACTS_ERROR("getNextTarget failed to find a valid target surface after "
               << state.options.maxTargetSkipping << " attempts.");
    return Result<NavigationTarget>::failure(
        PropagatorError::NextTargetLimitReached);
  };

  // priming error condition
  bool terminatedNormally = false;

  // Pre-Stepping: target setting
  state.stage = PropagatorStage::preStep;

  Result<NavigationTarget> nextTargetResult = getNextTarget();
  if (!nextTargetResult.ok()) {
    return nextTargetResult.error();
  }
  NavigationTarget nextTarget = *nextTargetResult;

  ACTS_VERBOSE("Starting stepping loop.");

  // Stepping loop
  for (; state.steps < state.options.maxSteps; ++state.steps) {
    // Perform a step
    Result<double> res =
        m_stepper.step(state.stepping, state.options.direction,
                       m_navigator.currentVolumeMaterial(state.navigation));
    if (!res.ok()) {
      ACTS_ERROR("Step failed with " << res.error() << ": "
                                     << res.error().message());
      // pass error to caller
      return res.error();
    }
    // Accumulate the path length
    state.pathLength += *res;
    // Update the position and direction
    state.position = m_stepper.position(state.stepping);
    state.direction =
        state.options.direction * m_stepper.direction(state.stepping);

    // Refresh the navigator's helix-model snapshot from the stepper so the
    // next nextTarget() / Layer::compatibleSurfaces call can use the
    // closed-form helix-plane predicate where applicable. The snapshot is
    // optional — when the stepper exposes qOverP and getField — and only
    // used by the navigator inside the inward-barrel regime.
    if constexpr (requires {
                    state.navigation.currentQOverP;
                    m_stepper.qOverP(state.stepping);
                    m_stepper.getField(state.stepping, state.position);
                  }) {
      state.navigation.currentQOverP = m_stepper.qOverP(state.stepping);
      auto fieldRes = m_stepper.getField(state.stepping, state.position);
      state.navigation.currentBField =
          fieldRes.ok() ? *fieldRes : Vector3::Zero();
    }

    ACTS_VERBOSE("Step with size " << *res << " performed. We are now at "
                                   << state.position.transpose()
                                   << " with direction "
                                   << state.direction.transpose());

    // release actor and aborter constrains after step was performed
    m_stepper.releaseStepSize(state.stepping, ConstrainedStep::Type::Navigator);
    m_stepper.releaseStepSize(state.stepping, ConstrainedStep::Type::Actor);

    // Post-stepping: check target status, call actors, check abort conditions
    state.stage = PropagatorStage::postStep;

    if (!nextTarget.isNone()) {
      IntersectionStatus postStepSurfaceStatus = m_stepper.updateSurfaceStatus(
          state.stepping, *nextTarget.surface,
          nextTarget.surfaceIntersectionIndex, state.options.direction,
          nextTarget.boundaryTolerance, state.options.surfaceTolerance,
          ConstrainedStep::Type::Navigator, state.navigation.isInBarrelVolume,
          logger());
      if (postStepSurfaceStatus == IntersectionStatus::onSurface) {
        m_navigator.handleSurfaceReached(state.navigation, state.position,
                                         state.direction, *nextTarget.surface);
      }
      if (postStepSurfaceStatus != IntersectionStatus::reachable) {
        // Snapshot the surface that just became unreachable for diagnostics
        // (nextTarget itself is reset to None below).
        const Surface* unreachableSurface = nextTarget.surface;
        nextTarget = NavigationTarget::None();
        // Check if target became unreachable due to turning point detection
        // Only if stepper supports turning point detection
        if constexpr (requires { state.stepping.turningPointDetected; }) {
          if (state.stepping.turningPointDetected) {
            if constexpr (requires { state.navigation.navLayers; }) {
              const double rxyTP = std::sqrt(
                  state.position[0] * state.position[0] +
                  state.position[1] * state.position[1]);
              const double zTP = state.position[2];
              const double phiTP = std::atan2(state.position[1],
                                              state.position[0]);
              ACTS_VERBOSE("TP[apex] rxy="
                           << rxyTP << " z=" << zTP << " phi=" << phiTP
                           << " radInward=" << state.navigation.radiallyInward
                           << " -> "
                           << state.stepping.radiallyInward_previous
                           << " stage="
                           << static_cast<int>(
                                  state.navigation.navigationStage));
              ACTS_VERBOSE("TP[apex] currentVolume="
                           << (state.navigation.currentVolume != nullptr
                                   ? state.navigation.currentVolume
                                         ->geometryId()
                                   : Acts::GeometryIdentifier{})
                           << " currentLayer="
                           << (state.navigation.currentLayer != nullptr
                                   ? state.navigation.currentLayer
                                         ->geometryId()
                                   : Acts::GeometryIdentifier{}));
              if (unreachableSurface != nullptr) {
                ACTS_VERBOSE("TP[apex] unreachable target geoId="
                             << unreachableSurface->geometryId()
                             << " sensitive="
                             << (unreachableSurface->geometryId().sensitive() !=
                                 0));
              }
              // Dump the current layer's approach surfaces (r-values) so we
              // can tell whether the apex is geometrically between apr=1
              // and apr=2 of the navigator's current layer.
              if (state.navigation.currentLayer != nullptr) {
                const auto* ad =
                    state.navigation.currentLayer->approachDescriptor();
                if (ad != nullptr) {
                  for (const auto* asurf : ad->containedSurfaces()) {
                    if (asurf == nullptr) continue;
                    // For cylinder approach surfaces the radius is encoded
                    // in surface bounds, not the translation. Try cast.
                    const auto& bounds = asurf->bounds();
                    auto vals = bounds.values();
                    std::ostringstream bs;
                    bs << "[";
                    for (std::size_t k = 0; k < vals.size(); ++k) {
                      if (k != 0) bs << ", ";
                      bs << vals[k];
                    }
                    bs << "]";
                    ACTS_VERBOSE("TP[apex] currentLayer approach geoId="
                                 << asurf->geometryId() << " bounds="
                                 << bs.str());
                  }
                }
              }
              // Push the post-flip direction into the navigator's bookkeeping.
              state.navigation.radiallyInward =
                  state.stepping.radiallyInward_previous;

              // Decide whether the apex sits geometrically INSIDE the
              // navigator's current sensor layer (rxy between apr=1 and
              // apr=2 r-radii). If so, the layerTarget fallback would
              // target apr=1 directly and the trajectory would step
              // through the sensors at the layer's mean r unnoticed (the
              // sensors are only added to navSurfaces AFTER apr=1 is
              // reached, by which point they're behind the trajectory and
              // get culled by negative path length). Force surfaceTarget
              // on the same layer so compatibleSurfaces is called from
              // the apex position with the new inward direction —
              // sensors at r ≈ apr=1 + 5mm are then ahead and reachable.
              bool apexInCurrentLayer = false;
              if (state.navigation.currentLayer != nullptr &&
                  state.navigation.currentLayer->surfaceArray() != nullptr) {
                const auto* ad =
                    state.navigation.currentLayer->approachDescriptor();
                if (ad != nullptr) {
                  double rmin = std::numeric_limits<double>::infinity();
                  double rmax = -std::numeric_limits<double>::infinity();
                  for (const auto* asurf : ad->containedSurfaces()) {
                    if (asurf == nullptr) continue;
                    auto vals = asurf->bounds().values();
                    if (!vals.empty()) {
                      const double r = vals[0];
                      rmin = std::min(rmin, r);
                      rmax = std::max(rmax, r);
                    }
                  }
                  if (rmin < rmax) {
                    apexInCurrentLayer = (rxyTP >= rmin && rxyTP <= rmax);
                  }
                }
              }

              // Env-gate the surfaceTarget fix so we can A/B with vs without.
              // ACTS_DISABLE_TP_INSIDE_LAYER_FIX=1 → always full reset +
              // layerTarget (the natural-baseline behaviour, used to measure
              // the un-fixed failure mode).
              static const bool tpInsideLayerFixDisabled =
                  std::getenv("ACTS_DISABLE_TP_INSIDE_LAYER_FIX") != nullptr;

              if (apexInCurrentLayer && !tpInsideLayerFixDisabled) {
                state.navigation.navSurfaces.clear();
                state.navigation.navSurfaceIndex.reset();
                state.navigation.navLayers.clear();
                state.navigation.navLayerIndex.reset();
                state.navigation.navBoundaries.clear();
                state.navigation.navBoundaryIndex.reset();
                state.navigation.navigationStage = N::Stage::surfaceTarget;
                ACTS_VERBOSE(
                    "TP[apex] apex INSIDE layer "
                    << state.navigation.currentLayer->geometryId()
                    << "; restart surfaceTarget");
              } else {
                // Apex outside currentLayer (or no current layer / fix
                // disabled). Full reset + layerTarget.
                state.navigation.resetAfterVolumeSwitch();
                state.navigation.navigationStage = N::Stage::layerTarget;
                ACTS_VERBOSE(
                    "TP[apex] apex "
                    << (apexInCurrentLayer ? "INSIDE (fix disabled)"
                                           : "OUTSIDE")
                    << " current layer; full reset, restart layerTarget");
              }
            }
            state.stepping.turningPointDetected = false;  // Reset flag
          }
        }
        // Note: sensitive-surface retarget on unreachable is NOT raised here.
        // The navigator handles it inline at navSurfaces-exhaustion time
        // (see Navigator::nextTarget). This lets the existing list-iteration
        // fallback try every original candidate first, and only then re-runs
        // Layer::compatibleSurfaces from the current (closer) position when
        // none of the original candidates produced a sensitive hit.
      }

    }

    state.options.actorList.act(state, m_stepper, m_navigator, logger());

    if (state.options.actorList.checkAbort(state, m_stepper, m_navigator,
                                           logger())) {
      terminatedNormally = true;
      break;
    }

    // Update the position and direction because actors might have changed it
    state.position = m_stepper.position(state.stepping);
    state.direction =
        state.options.direction * m_stepper.direction(state.stepping);

    // Pre-Stepping: target setting
    state.stage = PropagatorStage::preStep;

    if (!nextTarget.isNone() &&
        !m_navigator.checkTargetValid(state.navigation, state.position,
                                      state.direction)) {
      ACTS_VERBOSE("Target is not valid anymore.");
      nextTarget = NavigationTarget::None();
    }

    if (nextTarget.isNone()) {
      // navigator step constraint is not valid anymore
      m_stepper.releaseStepSize(state.stepping,
                                ConstrainedStep::Type::Navigator);

      nextTargetResult = getNextTarget();
      if (!nextTargetResult.ok()) {
        return nextTargetResult.error();
      }
      nextTarget = *nextTargetResult;
    }
  }  // end of stepping loop

  // check if we didn't terminate normally via aborters
  if (!terminatedNormally) {
    ACTS_ERROR("Propagation reached the step count limit of "
               << state.options.maxSteps << " (did " << state.steps
               << " steps)");
    return PropagatorError::StepCountLimitReached;
  }

  ACTS_VERBOSE("Stepping loop done.");

  state.stage = PropagatorStage::postPropagation;

  // Post-stepping call to the actor list
  state.options.actorList.act(state, m_stepper, m_navigator, logger());

  return Result<void>::success();
}

template <typename S, typename N>
template <typename parameters_t, typename propagator_options_t,
          typename path_aborter_t>
auto Acts::Propagator<S, N>::propagate(const parameters_t& start,
                                       const propagator_options_t& options,
                                       bool createFinalParameters) const
    -> Result<
        actor_list_t_result_t<StepperBoundTrackParameters,
                              typename propagator_options_t::actor_list_type>> {
  static_assert(std::copy_constructible<StepperBoundTrackParameters>,
                "return track parameter type must be copy-constructible");

  auto state = makeState<propagator_options_t, path_aborter_t>(options);

  auto initRes =
      initialize<decltype(state), parameters_t, path_aborter_t>(state, start);
  if (!initRes.ok()) {
    return initRes.error();
  }

  // Perform the actual propagation
  auto propagationResult = propagate(state);

  return makeResult(std::move(state), propagationResult, options,
                    createFinalParameters);
}

template <typename S, typename N>
template <typename parameters_t, typename propagator_options_t,
          typename target_aborter_t, typename path_aborter_t>
auto Acts::Propagator<S, N>::propagate(
    const parameters_t& start, const Surface& target,
    const propagator_options_t& options) const
    -> Result<
        actor_list_t_result_t<StepperBoundTrackParameters,
                              typename propagator_options_t::actor_list_type>> {
  static_assert(BoundTrackParametersConcept<parameters_t>,
                "Parameters do not fulfill bound parameters concept.");

  auto state =
      makeState<propagator_options_t, target_aborter_t, path_aborter_t>(
          target, options);

  auto initRes =
      initialize<decltype(state), parameters_t, path_aborter_t>(state, start);
  if (!initRes.ok()) {
    return initRes.error();
  }

  // Perform the actual propagation
  auto propagationResult = propagate(state);

  return makeResult(std::move(state), propagationResult, target, options);
}

template <typename S, typename N>
template <typename propagator_options_t, typename path_aborter_t>
auto Acts::Propagator<S, N>::makeState(
    const propagator_options_t& options) const {
  // Type of track parameters produced by the propagation
  using ReturnParameterType = StepperBoundTrackParameters;

  static_assert(std::copy_constructible<ReturnParameterType>,
                "return track parameter type must be copy-constructible");

  // Expand the actor list with a path aborter
  path_aborter_t pathAborter;
  pathAborter.internalLimit = options.pathLimit;

  auto actorList = options.actorList.append(pathAborter);

  // Create the extended options and declare their type
  auto eOptions = options.extend(actorList);

  using OptionsType = decltype(eOptions);
  using StateType =
      actor_list_t_state_t<OptionsType,
                           typename propagator_options_t::actor_list_type>;

  StateType state{eOptions, m_stepper.makeState(eOptions.stepping),
                  m_navigator.makeState(eOptions.navigation)};

  return state;
}

template <typename S, typename N>
template <typename propagator_options_t, typename target_aborter_t,
          typename path_aborter_t>
auto Acts::Propagator<S, N>::makeState(
    const Surface& target, const propagator_options_t& options) const {
  // Expand the actor list with a target and path aborter
  target_aborter_t targetAborter;
  targetAborter.surface = &target;
  path_aborter_t pathAborter;
  pathAborter.internalLimit = options.pathLimit;

  auto actorList = options.actorList.append(targetAborter, pathAborter);

  // Create the extended options and declare their type
  auto eOptions = options.extend(actorList);
  eOptions.navigation.targetSurface = &target;

  using OptionsType = decltype(eOptions);
  using StateType =
      actor_list_t_state_t<OptionsType,
                           typename propagator_options_t::actor_list_type>;

  StateType state{eOptions, m_stepper.makeState(eOptions.stepping),
                  m_navigator.makeState(eOptions.navigation)};

  return state;
}

template <typename S, typename N>
template <typename propagator_state_t, typename parameters_t,
          typename path_aborter_t>
Acts::Result<void> Acts::Propagator<S, N>::initialize(
    propagator_state_t& state, const parameters_t& start) const {
  static_assert(BoundTrackParametersConcept<parameters_t>,
                "Parameters do not fulfill bound parameters concept.");

  m_stepper.initialize(state.stepping, start.toBound());

  state.position = m_stepper.position(state.stepping);
  state.direction =
      state.options.direction * m_stepper.direction(state.stepping);

  state.navigation.options.startSurface = &start.referenceSurface();

  // Navigator initialize state call
  auto navInitRes =
      m_navigator.initialize(state.navigation, state.position, state.direction,
                             state.options.direction);
  if (!navInitRes.ok()) {
    return navInitRes.error();
  }

  // Apply the loop protection - it resets the internal path limit
  detail::setupLoopProtection(
      state, m_stepper, state.options.actorList.template get<path_aborter_t>(),
      false, logger());

  return Result<void>::success();
}

template <typename S, typename N>
template <typename propagator_state_t, typename propagator_options_t>
auto Acts::Propagator<S, N>::makeResult(propagator_state_t state,
                                        Result<void> propagationResult,
                                        const propagator_options_t& /*options*/,
                                        bool createFinalParameters) const
    -> Result<
        actor_list_t_result_t<StepperBoundTrackParameters,
                              typename propagator_options_t::actor_list_type>> {
  // Type of track parameters produced by the propagation
  using ReturnParameterType = StepperBoundTrackParameters;

  static_assert(std::copy_constructible<ReturnParameterType>,
                "return track parameter type must be copy-constructible");

  // Type of the full propagation result, including output from actors
  using ResultType =
      actor_list_t_result_t<ReturnParameterType,
                            typename propagator_options_t::actor_list_type>;

  if (!propagationResult.ok()) {
    return propagationResult.error();
  }

  ResultType result{};
  moveStateToResult(state, result);

  const Surface* currentSurface = m_navigator.currentSurface(state.navigation);
  if (createFinalParameters && currentSurface == nullptr) {
    if (!m_stepper.prepareCurvilinearState(state.stepping)) {
      // information to compute curvilinearState is incomplete.
      return propagationResult.error();
    }
    /// Convert into return type and fill the result object
    auto curvState = m_stepper.curvilinearState(state.stepping);
    // Fill the end parameters
    result.endParameters = std::get<StepperBoundTrackParameters>(curvState);
    // Only fill the transport jacobian when covariance transport was done
    if (state.stepping.covTransport) {
      result.transportJacobian = std::get<Jacobian>(curvState);
    }
  } else if (createFinalParameters && currentSurface != nullptr) {
    // We are at a surface, so we need to compute the bound state
    auto boundState =
        m_stepper.boundState(state.stepping, *currentSurface).value();
    // Fill the end parameters
    result.endParameters = std::get<StepperBoundTrackParameters>(boundState);
    // Only fill the transport jacobian when covariance transport was done
    if (state.stepping.covTransport) {
      result.transportJacobian = std::get<Jacobian>(boundState);
    }
  }

  return Result<ResultType>::success(std::move(result));
}

template <typename S, typename N>
template <typename propagator_state_t, typename propagator_options_t>
auto Acts::Propagator<S, N>::makeResult(
    propagator_state_t state, Result<void> propagationResult,
    const Surface& target, const propagator_options_t& /*options*/) const
    -> Result<
        actor_list_t_result_t<StepperBoundTrackParameters,
                              typename propagator_options_t::actor_list_type>> {
  // Type of track parameters produced at the end of the propagation
  using ReturnParameterType = StepperBoundTrackParameters;

  static_assert(std::copy_constructible<ReturnParameterType>,
                "return track parameter type must be copy-constructible");

  // Type of the full propagation result, including output from actors
  using ResultType =
      actor_list_t_result_t<ReturnParameterType,
                            typename propagator_options_t::actor_list_type>;

  if (!propagationResult.ok()) {
    return propagationResult.error();
  }

  ResultType result{};
  moveStateToResult(state, result);

  // Compute the final results and mark the propagation as successful
  auto bsRes = m_stepper.boundState(state.stepping, target);
  if (!bsRes.ok()) {
    return bsRes.error();
  }
  const auto& bs = *bsRes;

  // Fill the end parameters
  result.endParameters = std::get<StepperBoundTrackParameters>(bs);
  // Only fill the transport jacobian when covariance transport was done
  if (state.stepping.covTransport) {
    result.transportJacobian = std::get<Jacobian>(bs);
  }
  return Result<ResultType>::success(std::move(result));
}

template <typename S, typename N>
template <typename propagator_state_t, typename result_t>
void Acts::Propagator<S, N>::moveStateToResult(propagator_state_t& state,
                                               result_t& result) const {
  result.tuple() = std::move(state.tuple());

  result.steps = state.steps;
  result.pathLength = state.pathLength;

  result.statistics.stepping = state.stepping.statistics;
  result.statistics.navigation = state.navigation.statistics;
}

template <typename derived_t>
Acts::Result<Acts::BoundTrackParameters>
Acts::detail::BasePropagatorHelper<derived_t>::propagateToSurface(
    const BoundTrackParameters& start, const Surface& target,
    const Options& options) const {
  using ResultType = Result<typename derived_t::template actor_list_t_result_t<
      BoundTrackParameters, ActorList<>>>;
  using DerivedOptions = typename derived_t::template Options<>;

  DerivedOptions derivedOptions(options);

  // dummy initialization
  ResultType res = ResultType::failure(PropagatorError::Failure);

  // Due to the geometry of the perigee surface the overstepping tolerance
  // is sometimes not met.
  if (target.type() == Surface::SurfaceType::Perigee) {
    res = static_cast<const derived_t*>(this)
              ->template propagate<BoundTrackParameters, DerivedOptions,
                                   ForcedSurfaceReached, PathLimitReached>(
                  start, target, derivedOptions);
  } else {
    res = static_cast<const derived_t*>(this)
              ->template propagate<BoundTrackParameters, DerivedOptions,
                                   SurfaceReached, PathLimitReached>(
                  start, target, derivedOptions);
  }

  if (!res.ok()) {
    return res.error();
  }

  // Without errors we can expect a valid endParameters when propagating to a
  // target surface
  assert((*res).endParameters);
  return std::move((*res).endParameters.value());
}
