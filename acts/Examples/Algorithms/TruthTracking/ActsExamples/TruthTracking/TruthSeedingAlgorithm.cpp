// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TruthSeedingAlgorithm.hpp"

#include "Acts/Definitions/Common.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace ActsExamples {

TruthSeedingAlgorithm::TruthSeedingAlgorithm(Config cfg,
                                             Acts::Logging::Level lvl)
    : IAlgorithm("TruthSeedingAlgorithm", lvl), m_cfg(std::move(cfg)) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputParticleMeasurementsMap.empty()) {
    throw std::invalid_argument(
        "Missing input particle-measurements map collection");
  }
  if (m_cfg.inputSpacePoints.empty()) {
    throw std::invalid_argument("Missing seeds or space point collection");
  }

  for (const auto& spName : m_cfg.inputSpacePoints) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto& handle = m_inputSpacePoints.emplace_back(
        std::make_unique<ReadDataHandle<SimSpacePointContainer>>(
            this,
            "InputSpacePoints#" + std::to_string(m_inputSpacePoints.size())));
    handle->initialize(spName);
  }

  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument("Missing output particles collection");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing seeds output collections");
  }
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing proto tracks output collections");
  }

  if (m_cfg.inputSimHits.empty()) {
    throw std::invalid_argument("Missing input simulated hits collection");
  }

  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing input simulated hits measurements map");
  }

  m_inputParticles.initialize(m_cfg.inputParticles);
  m_inputParticleMeasurementsMap.initialize(m_cfg.inputParticleMeasurementsMap);
  m_outputParticles.initialize(m_cfg.outputParticles);
  m_outputProtoTracks.initialize(m_cfg.outputProtoTracks);
  m_outputSeeds.initialize(m_cfg.outputSeeds);

  m_inputSimHits.initialize(m_cfg.inputSimHits);
  m_inputMeasurementSimHitsMap.initialize(m_cfg.inputMeasurementSimHitsMap);
}

ProcessCode TruthSeedingAlgorithm::execute(const AlgorithmContext& ctx) const {
  // prepare input collections
  const auto& particles = m_inputParticles(ctx);
  const auto& particleMeasurementsMap = m_inputParticleMeasurementsMap(ctx);
  const auto& simHits = m_inputSimHits(ctx);
  const auto& measurementSimHitsMap = m_inputMeasurementSimHitsMap(ctx);

  // construct the combined input container of space point pointers from all
  // configured input sources.
  // pre-compute the total size required so we only need to allocate once
  std::size_t nSpacePoints = 0;
  for (const auto& isp : m_inputSpacePoints) {
    nSpacePoints += (*isp)(ctx).size();
  }

  std::vector<const SimSpacePoint*> spacePointPtrs;
  spacePointPtrs.reserve(nSpacePoints);
  for (const auto& isp : m_inputSpacePoints) {
    for (const auto& spacePoint : (*isp)(ctx)) {
      // since the event store owns the space points, their pointers should be
      // stable and we do not need to create local copies.
      spacePointPtrs.push_back(&spacePoint);
    }
  }

  SimParticleContainer seededParticles;
  SimSeedContainer seeds;
  ProtoTrackContainer tracks;

  seededParticles.reserve(particles.size());
  seeds.reserve(particles.size());
  tracks.reserve(particles.size());

  std::unordered_map<Index, const SimSpacePoint*> spMap;

  for (const auto& spp : spacePointPtrs) {
    if (spp->sourceLinks().empty()) {
      ACTS_WARNING("Missing source link in space point");
      continue;
    }
    for (const auto& slink : spp->sourceLinks()) {
      const IndexSourceLink& islink = slink.get<IndexSourceLink>();
      spMap.emplace(islink.index(), spp);
    }
  }

  for (const auto& particle : particles) {
    // find the corresponding measurements for this particle
    const auto& measurements =
        makeRange(particleMeasurementsMap.equal_range(particle.particleId()));
    // fill measurement indices to create the proto track
    ProtoTrack track;
    track.reserve(measurements.size());

    std::vector<std::pair<const SimHit*, Index>> hits;
    hits.reserve(measurements.size());

    for (const auto& [barcode, index] : measurements) {
      const auto simHitMapIt = measurementSimHitsMap.find(index);
      if (simHitMapIt == measurementSimHitsMap.end()) {
        ACTS_WARNING("No sim hit found for measurement index " << index);
        continue;
      }

      const auto simHitIt = simHits.nth(simHitMapIt->second);
      if (simHitIt == simHits.end()) {
        ACTS_WARNING("No sim hit found for index " << simHitMapIt->second);
        continue;
      }

      const auto& simHit = *simHitIt;

      hits.emplace_back(&simHit, index);
    }

    std::sort(hits.begin(), hits.end(), [](const auto& a, const auto& b) {
      return a.first->time() < b.first->time();
    });

    for (const auto& [hit, index] : hits) {
      track.push_back(index);
    }

    // The list of measurements and the initial start parameters
    if (track.size() < 3) {
      ACTS_WARNING("Particle " << particle << " has less than 3 measurements");
      continue;
    }

    // Walk the hit-time-ordered measurement list, look up each SP, and
    // classify the SP by the sign of the truth radial momentum at that hit:
    //   p_r = (p_x * x + p_y * y) / r
    //   p_r > 0  -> outgoing arm (particle moving radially outward)
    //   p_r < 0  -> incoming arm (particle moving radially inward)
    //
    // This is a direct, per-hit physical label using SimHit::momentum4Before
    // and avoids the r-monotonicity edge cases (apex, pericenter, equal-r
    // SPs across arms).
    //
    // Low-pT loopers leave SPs on multiple arms (outgoing -> turnaround ->
    // incoming -> outgoing again -> ...). Each maximal run of same-sign p_r
    // is one "arm" and produces its own seed below, so a single looper can
    // contribute multiple seeds (one per arm) instead of being thrown out.
    std::vector<std::vector<const SimSpacePoint*>> arms;
    int armDir = 0;  // +1 outgoing, -1 incoming
    for (const auto& [simHit, measurementIndex] : hits) {
      auto spIt = spMap.find(measurementIndex);
      if (spIt == spMap.end()) {
        continue;
      }
      const SimSpacePoint* sp = spIt->second;

      const auto& pos4 = simHit->fourPosition();
      const auto& mom4 = simHit->momentum4Before();
      const double x = pos4[Acts::ePos0];
      const double y = pos4[Acts::ePos1];
      const double r = std::hypot(x, y);
      if (r == 0) {
        continue;
      }
      const double pr =
          (mom4[Acts::eMom0] * x + mom4[Acts::eMom1] * y) / r;
      if (pr == 0) {
        // exactly tangential motion at this hit -- ambiguous, skip
        continue;
      }
      const int dir = (pr > 0) ? +1 : -1;

      if (armDir == 0 || dir != armDir) {
        // first SP, or direction flipped -> start a new arm
        arms.emplace_back();
        armDir = dir;
      }
      arms.back().push_back(sp);
    }

    bool anySeedFound = false;
    for (auto& arm : arms) {
      if (arm.size() < 3) {
        continue;
      }

      // Re-order each arm by r ascending so bottom = innermost SP and
      // top = outermost SP, regardless of whether the arm was outgoing
      // (time-order already r-ascending) or incoming (time-order
      // r-descending). Both seeds therefore present an outward-going
      // (bottom -> middle -> top) helix to TrackParamsEstimationAlgorithm,
      // matching ACTS' standard seed convention.
      std::ranges::sort(arm, [](const auto* a, const auto* b) {
        return a->r() < b->r();
      });

      // Loop over the SPs on this arm to find the triplet with maximum
      // deltaR between the bottom-middle and middle-top pairs.
      // @todo add the check of deltaZ
      bool seedFound = false;
      std::array<std::size_t, 3> bestSPIndices{};
      double maxDeltaR = std::numeric_limits<double>::min();
      for (std::size_t ib = 0; ib < arm.size() - 2; ++ib) {
        for (std::size_t im = ib + 1; im < arm.size() - 1; ++im) {
          for (std::size_t it = im + 1; it < arm.size(); ++it) {
            double bmDeltaR = arm[im]->r() - arm[ib]->r();
            double mtDeltaR = arm[it]->r() - arm[im]->r();
            if (bmDeltaR >= m_cfg.deltaRMin && bmDeltaR <= m_cfg.deltaRMax &&
                mtDeltaR >= m_cfg.deltaRMin && mtDeltaR <= m_cfg.deltaRMax &&
                (bmDeltaR + mtDeltaR) > maxDeltaR) {
              maxDeltaR = bmDeltaR + mtDeltaR;
              bestSPIndices = {ib, im, it};
              seedFound = true;
            }
          }
        }
      }

      if (seedFound) {
        SimSeed seed{*arm[bestSPIndices[0]], *arm[bestSPIndices[1]],
                     *arm[bestSPIndices[2]]};
        seed.setVertexZ(static_cast<float>(arm[bestSPIndices[1]]->z()));

        seeds.emplace_back(seed);
        // Mirror the original 1:1 seed/proto-track output by emitting a
        // copy of the particle's full hit list per seed. Downstream
        // (SeedsToPrototracks) only consumes seeds, so this is mostly a
        // bookkeeping copy.
        tracks.emplace_back(track);
        anySeedFound = true;
      }
    }

    if (anySeedFound) {
      seededParticles.insert(particle);
    }
  }

  ACTS_VERBOSE("Found " << seeds.size() << " seeds");

  m_outputParticles(ctx, std::move(seededParticles));
  m_outputProtoTracks(ctx, std::move(tracks));
  m_outputSeeds(ctx, std::move(seeds));

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
