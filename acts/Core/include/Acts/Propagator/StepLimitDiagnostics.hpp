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
inline void resetStepLimitCounts() {
  fatrasFailedCounter().store(0, std::memory_order_relaxed);
  fatrasStepLimitCounter().store(0, std::memory_order_relaxed);
  trackFindingFailedCounter().store(0, std::memory_order_relaxed);
  trackFindingStepLimitCounter().store(0, std::memory_order_relaxed);
}

}  // namespace Acts::detail
