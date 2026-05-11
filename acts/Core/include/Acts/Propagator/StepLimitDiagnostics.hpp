// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>

namespace Acts::detail {

/// Process-global counters, split by stage (Fatras vs CKF) and by failure
/// type (any propagator error vs the StepCountLimitReached subset).
/// Storage lives in libActsCore (Core/src/Utilities/Intersection.cpp) so every
/// DSO sees the same atomic instances.
///
/// "Failed" counters increment on ANY non-ok propagator result; "StepLimit"
/// counters increment only on PropagatorError::StepCountLimitReached. The
/// step-limit counter is therefore a strict subset of the total-failed one.
std::atomic<std::size_t>& fatrasFailedCounter();
std::atomic<std::size_t>& fatrasStepLimitCounter();
std::atomic<std::size_t>& trackFindingFailedCounter();
std::atomic<std::size_t>& trackFindingStepLimitCounter();

/// Sensitive-surface retarget counters (Navigator-internal).
/// "Invoked" increments each time the navigator clears navSurfaces and re-runs
/// Layer::compatibleSurfaces in response to an unreachable verdict on a
/// sensitive module. "Succeeded" increments at most once per layer entry, the
/// first time the navigator subsequently reports onSurface for any sensitive
/// module on that layer. Split by stage so end-of-run prints can show Fatras
/// vs CKF separately.
std::atomic<std::size_t>& fatrasRetargetInvokedCounter();
std::atomic<std::size_t>& fatrasRetargetSucceededCounter();
std::atomic<std::size_t>& ckfRetargetInvokedCounter();
std::atomic<std::size_t>& ckfRetargetSucceededCounter();

inline std::size_t fatrasFailedCount() {
  return fatrasFailedCounter().load(std::memory_order_relaxed);
}
inline std::size_t fatrasStepLimitCount() {
  return fatrasStepLimitCounter().load(std::memory_order_relaxed);
}
inline std::size_t trackFindingFailedCount() {
  return trackFindingFailedCounter().load(std::memory_order_relaxed);
}
inline std::size_t trackFindingStepLimitCount() {
  return trackFindingStepLimitCounter().load(std::memory_order_relaxed);
}
inline std::size_t fatrasRetargetInvokedCount() {
  return fatrasRetargetInvokedCounter().load(std::memory_order_relaxed);
}
inline std::size_t fatrasRetargetSucceededCount() {
  return fatrasRetargetSucceededCounter().load(std::memory_order_relaxed);
}
inline std::size_t ckfRetargetInvokedCount() {
  return ckfRetargetInvokedCounter().load(std::memory_order_relaxed);
}
inline std::size_t ckfRetargetSucceededCount() {
  return ckfRetargetSucceededCounter().load(std::memory_order_relaxed);
}
/// Per-layer apex-inside-shell counters (temporary diagnostic).
/// Increment when the apex turning point fires geometrically inside a barrel
/// layer's apr=1/apr=2 shell; the per-track latch in Navigator::State then
/// triggers `bumpApexInsideShellRecovered` on the first subsequent sensor
/// hit on the same layer. `not_recovered` per layer = inside - recovered.
void bumpApexInsideShell(std::uint32_t volumeId, std::uint32_t layerId);
void bumpApexInsideShellRecovered(std::uint32_t volumeId,
                                  std::uint32_t layerId);

/// Snapshot of per-layer counters as a vector of (vol, lay, inside,
/// recovered) tuples. Used by Python bindings + end-of-run print.
std::vector<std::tuple<std::uint32_t, std::uint32_t, std::size_t, std::size_t>>
apexInsideShellSnapshot();

/// Zero out the per-layer apex-inside-shell map.
void resetApexInsideShellCounts();

inline void resetStepLimitCounts() {
  fatrasFailedCounter().store(0, std::memory_order_relaxed);
  fatrasStepLimitCounter().store(0, std::memory_order_relaxed);
  trackFindingFailedCounter().store(0, std::memory_order_relaxed);
  trackFindingStepLimitCounter().store(0, std::memory_order_relaxed);
  fatrasRetargetInvokedCounter().store(0, std::memory_order_relaxed);
  fatrasRetargetSucceededCounter().store(0, std::memory_order_relaxed);
  ckfRetargetInvokedCounter().store(0, std::memory_order_relaxed);
  ckfRetargetSucceededCounter().store(0, std::memory_order_relaxed);
  resetApexInsideShellCounts();
}

}  // namespace Acts::detail
