/*
 * Copyright 2018 The Cartographer Authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef CARTOGRAPHER_MAPPING_TSDF_VALUES_H_
#define CARTOGRAPHER_MAPPING_TSDF_VALUES_H_

#include <cmath>
#include <vector>

#include "cartographer/common/math.h"
#include "cartographer/common/port.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping {

constexpr float kMinTSDF = -0.3f;
constexpr float kMaxTSDF = -kMinTSDF;

constexpr float kMinWeight = 0.0f;
constexpr float kMaxWeight = 100.0f;

// Clamps probability to be in the range [kMinTSDF, kMaxTSDF].
inline float ClampTSDF(const float tsdf) {
  return common::Clamp(tsdf, kMinTSDF, kMaxTSDF);
}
// Clamps probability to be in the range [kMinTSDF, kMaxTSDF].
inline float ClampWeight(const float tsdf) {
  return common::Clamp(tsdf, kMinWeight, kMaxWeight);
}

constexpr uint16 kUnknownTSDFValue = 0;
constexpr uint16 kUnknownWeightValue = 0;
constexpr uint16 kUpdateMarker = 1u << 15;

// Converts a tsdf to a uint16 in the [1, 32767] range.
inline uint16 TSDFToValue(const float tsdf) {
  const int value = common::RoundToInt((ClampTSDF(tsdf) - kMinTSDF) *
                                       (32766.f / (kMaxTSDF - kMinTSDF))) +
                    1;
  // DCHECK for performance.
  DCHECK_GE(value, 1);
  DCHECK_LE(value, 32767);
  return value;
}

// Converts a weight to a uint16 in the [1, 32767] range.
inline uint16 WeightToValue(const float weight) {
  const int value = common::RoundToInt((ClampWeight(weight) - kMinWeight) *
                                       (32766.f / (kMaxWeight - kMinWeight))) +
                    1;
  // DCHECK for performance.
  DCHECK_GE(value, 1);
  DCHECK_LE(value, 32767);
  return value;
}

extern const std::vector<float>* const kValueToTSDF;
extern const std::vector<float>* const kValueToWeight;

// Converts a uint16 (which may or may not have the update marker set) to a
// probability in the range [kMinProbability, kMaxProbability].
inline float ValueToTSDF(const uint16 value) {
  return (*kValueToTSDF)[value];

}  // Converts a uint16 (which may or may not have the update marker set) to a
// probability in the range [kMinProbability, kMaxProbability].
inline float ValueToWeight(const uint16 value) {
  return (*kValueToWeight)[value];
}

}  // namespace mapping
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_TSDF_VALUES_H_
