// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Intersection.hpp"

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/Propagator/StepLimitDiagnostics.hpp"

#include <map>
#include <memory>
#include <mutex>
#include <utility>

namespace Acts {

namespace {
struct ApexInsideShellEntry {
  std::atomic<std::size_t> inside{0};
  std::atomic<std::size_t> recovered{0};
};

// Per-(volume, layer) atomic counters. The map is grown under a mutex on
// first use of each key; subsequent increments are atomic and lock-free.
std::mutex& apexInsideShellMutex() {
  static std::mutex m;
  return m;
}
std::map<std::pair<std::uint32_t, std::uint32_t>,
         std::unique_ptr<ApexInsideShellEntry>>&
apexInsideShellStore() {
  static std::map<std::pair<std::uint32_t, std::uint32_t>,
                  std::unique_ptr<ApexInsideShellEntry>>
      store;
  return store;
}
ApexInsideShellEntry& apexInsideShellEntry(std::uint32_t vol,
                                           std::uint32_t lay) {
  auto key = std::make_pair(vol, lay);
  std::lock_guard<std::mutex> g(apexInsideShellMutex());
  auto& slot = apexInsideShellStore()[key];
  if (!slot) {
    slot = std::make_unique<ApexInsideShellEntry>();
  }
  return *slot;
}
}  // namespace

void detail::bumpApexInsideShell(std::uint32_t vol, std::uint32_t lay) {
  apexInsideShellEntry(vol, lay).inside.fetch_add(1, std::memory_order_relaxed);
}
void detail::bumpApexInsideShellRecovered(std::uint32_t vol,
                                          std::uint32_t lay) {
  apexInsideShellEntry(vol, lay).recovered.fetch_add(1,
                                                     std::memory_order_relaxed);
}
std::vector<std::tuple<std::uint32_t, std::uint32_t, std::size_t, std::size_t>>
detail::apexInsideShellSnapshot() {
  std::vector<std::tuple<std::uint32_t, std::uint32_t, std::size_t,
                          std::size_t>>
      out;
  std::lock_guard<std::mutex> g(apexInsideShellMutex());
  for (const auto& [key, ent] : apexInsideShellStore()) {
    out.emplace_back(key.first, key.second,
                     ent->inside.load(std::memory_order_relaxed),
                     ent->recovered.load(std::memory_order_relaxed));
  }
  return out;
}
void detail::resetApexInsideShellCounts() {
  std::lock_guard<std::mutex> g(apexInsideShellMutex());
  for (auto& [key, ent] : apexInsideShellStore()) {
    ent->inside.store(0, std::memory_order_relaxed);
    ent->recovered.store(0, std::memory_order_relaxed);
  }
}


std::atomic<std::size_t>& detail::fatrasFailedCounter() {
  static std::atomic<std::size_t> counter{0};
  return counter;
}
std::atomic<std::size_t>& detail::fatrasStepLimitCounter() {
  static std::atomic<std::size_t> counter{0};
  return counter;
}
std::atomic<std::size_t>& detail::trackFindingFailedCounter() {
  static std::atomic<std::size_t> counter{0};
  return counter;
}
std::atomic<std::size_t>& detail::trackFindingStepLimitCounter() {
  static std::atomic<std::size_t> counter{0};
  return counter;
}
std::atomic<std::size_t>& detail::fatrasRetargetInvokedCounter() {
  static std::atomic<std::size_t> counter{0};
  return counter;
}
std::atomic<std::size_t>& detail::fatrasRetargetSucceededCounter() {
  static std::atomic<std::size_t> counter{0};
  return counter;
}
std::atomic<std::size_t>& detail::ckfRetargetInvokedCounter() {
  static std::atomic<std::size_t> counter{0};
  return counter;
}
std::atomic<std::size_t>& detail::ckfRetargetSucceededCounter() {
  static std::atomic<std::size_t> counter{0};
  return counter;
}

bool detail::checkPathLength(double pathLength, double nearLimit,
                             double farLimit, const Logger& logger) {
  // TODO why?
  const double tolerance = s_onSurfaceTolerance;

  ACTS_VERBOSE(" -> near limit, far limit, distance: "
               << nearLimit << ", " << farLimit << ", " << pathLength);

  const bool coCriterion = pathLength > nearLimit;
  const bool cpCriterion = pathLength < farLimit + tolerance;

  const bool accept = coCriterion && cpCriterion;

  if (accept) {
    ACTS_VERBOSE("Intersection is WITHIN limit");
  } else {
    ACTS_VERBOSE("Intersection is OUTSIDE limit because: ");
    if (!coCriterion) {
      ACTS_VERBOSE("- intersection path length "
                   << pathLength << " <= near limit " << nearLimit);
    }
    if (!cpCriterion) {
      ACTS_VERBOSE("- intersection path length "
                   << pathLength << " is over the far limit "
                   << (farLimit + tolerance) << " (including tolerance of "
                   << tolerance << ")");
    }
  }

  return accept;
}

}  // namespace Acts
