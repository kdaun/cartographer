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

#include "cartographer/mapping/2d/range_data_inserter_2d_tsdf.h"

#include <cstdlib>

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "cartographer/mapping/2d/xy_index.h"
#include "cartographer/mapping/internal/2d/ray_casting.h"
#include "cartographer/mapping/tsdf_values.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping {

RangeDataInserter2DTSDF::RangeDataInserter2DTSDF(
    const proto::RangeDataInserterOptions2D& options)
    : options_(options) {}

void RangeDataInserter2DTSDF::Insert(const sensor::RangeData& range_data,
                                     TSDF2D* const tsdf) const {
  const float truncation_distance = options_.tsdf().truncation_distance();
  Eigen::AlignedBox2f bounding_box(range_data.origin.head<2>());
  for (const Eigen::Vector3f& hit : range_data.returns) {
    const Eigen::Vector3f direction = (hit - range_data.origin).normalized();
    const Eigen::Vector3f end_position = hit + truncation_distance * direction;
    bounding_box.extend(end_position.head<2>());
  }
  constexpr float kPadding = 1e-6f;
  tsdf->GrowLimits(bounding_box.min() - kPadding * Eigen::Vector2f::Ones());
  tsdf->GrowLimits(bounding_box.max() + kPadding * Eigen::Vector2f::Ones());

  const Eigen::Vector2f origin = range_data.origin.head<2>();
  const Eigen::Array2i origin_cell = tsdf->limits().GetCellIndex(origin);
  for (const Eigen::Vector3f& hit_3d : range_data.returns) {
    const Eigen::Vector2f hit = hit_3d.head<2>();
    const Eigen::Vector2f direction = (hit - origin).normalized();
    float ray_length = (hit - origin).norm();
    const Eigen::Vector2f end_position = hit + truncation_distance * direction;
    const Eigen::Array2i end_cell = tsdf->limits().GetCellIndex(end_position);
    const Eigen::Array2i delta = end_cell - origin_cell;
    const int num_samples = delta.cwiseAbs().maxCoeff();
    CHECK_LT(num_samples, 1 << 15);
    // 'num_samples' is the number of samples we equi-distantly place on the
    // line between 'origin' and 'hit'. (including a fractional part for sub-
    // voxels) It is chosen so that between two samples we change from one voxel
    // to the next on the fastest changing dimension.
    for (int position = 0; position < num_samples; ++position) {
      const Eigen::Array2i cell = origin_cell + delta * position / num_samples;
      // float resolution = tsdf->limits().resolution();
      const Eigen::Vector2f cell_position = tsdf->limits().GetCellCenter(cell);
      float distance = (hit - cell_position).dot(direction);
      if (distance > truncation_distance) {
        distance = truncation_distance;
      } else if (distance < -truncation_distance) {
        distance = -truncation_distance;
      }
      UpdateCell(tsdf, cell, distance, ray_length);
    }
  }
  tsdf->FinishUpdate();
}

void RangeDataInserter2DTSDF::UpdateCell(TSDF2D* const tsdf,
                                         const Eigen::Array2i& cell,
                                         float update_sdf,
                                         float ray_length) const {
  float update_weight = 0.f;
  if (proto::CONSTANT_WEIGHT == options_.tsdf().range_data_inserter_type()) {
    update_weight = ComputeWeightConstant(update_sdf);
  } else if (proto::LINEAR_WEIGHT ==
      options_.tsdf().range_data_inserter_type()) {
    update_weight = ComputeWeightLinear(update_sdf, ray_length);
  } else if (proto::QUADRATIC_WEIGHT ==
      options_.tsdf().range_data_inserter_type()) {
    update_weight = ComputeWeightQuadratic(update_sdf, ray_length);
  }

  float updated_weight = tsdf->GetWeight(cell) + update_weight;
  float updated_sdf = updated_weight > 0.f
                          ? (tsdf->GetTSDF(cell) * tsdf->GetWeight(cell) +
                             update_sdf * update_weight) /
                                updated_weight
                          : tsdf->GetTSDF(cell);
  tsdf->UpdateCell(cell, updated_sdf, updated_weight);
}

float RangeDataInserter2DTSDF::ComputeWeightConstant(float sdf) const {
  float behind_surface_factor = 1.0;
  if (-options_.tsdf().behind_surface_distance() > sdf) {
    behind_surface_factor = (sdf + options_.tsdf().truncation_distance()) /
                            (options_.tsdf().truncation_distance() -
                             options_.tsdf().behind_surface_distance());
  }
  return behind_surface_factor * options_.tsdf().update_weight();
}

float RangeDataInserter2DTSDF::ComputeWeightLinear(float sdf,
                                                      float ray_length) const {
  float weight = 0.f;
  if (-options_.tsdf().behind_surface_distance() < sdf) {
    weight = options_.tsdf().update_weight() / ray_length;
  } else if (-options_.tsdf().truncation_distance() < sdf) {
    float behind_surface_factor =
        (sdf + options_.tsdf().truncation_distance()) /
            (options_.tsdf().truncation_distance() -
                options_.tsdf().behind_surface_distance());
    weight = behind_surface_factor * options_.tsdf().update_weight() /
        ray_length;
  }
  return weight;
}

float RangeDataInserter2DTSDF::ComputeWeightQuadratic(float sdf,
                                                      float ray_length) const {
  float weight = 0.f;
  if (-options_.tsdf().behind_surface_distance() < sdf) {
    weight = options_.tsdf().update_weight() / (std::pow(ray_length, 2));
  } else if (-options_.tsdf().truncation_distance() < sdf) {
    float behind_surface_factor =
        (sdf + options_.tsdf().truncation_distance()) /
        (options_.tsdf().truncation_distance() -
         options_.tsdf().behind_surface_distance());
    weight = behind_surface_factor * options_.tsdf().update_weight() /
             (std::pow(ray_length, 2));
  }
  return weight;
}

}  // namespace mapping
}  // namespace cartographer
