// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/BoundaryTolerance.hpp"

namespace Acts {

class Surface;

/// @brief The navigation target
///
/// This struct represents a navigation target which is communicated from the
/// navigator to the stepper through the propagator.
///
/// @note This incorporates `std::optional` semantics as the next target might
///       not exist.
struct NavigationTarget {
  const Surface* surface = nullptr;
  std::uint8_t surfaceIntersectionIndex = 0;
  BoundaryTolerance boundaryTolerance = BoundaryTolerance::None();
  /// True iff the navigator issued this target while in shell mode
  /// (trajectory descending through a sensor barrel layer's shell with no
  /// sensor hit yet). The stepper/SteppingHelper uses this to bypass the
  /// helix-plane intersect for sensitives and use a pure radial-inward
  /// line predicate instead — "just look down at the next module".
  bool shellMode = false;

  static NavigationTarget None() { return NavigationTarget(); }

  NavigationTarget(const Surface& surface_,
                   std::uint8_t surfaceIntersectionIndex_,
                   BoundaryTolerance boundaryTolerance_,
                   bool shellMode_ = false)
      : surface(&surface_),
        surfaceIntersectionIndex(surfaceIntersectionIndex_),
        boundaryTolerance(boundaryTolerance_),
        shellMode(shellMode_) {}

  bool isNone() const { return surface == nullptr; }

 private:
  NavigationTarget() = default;
};

}  // namespace Acts
