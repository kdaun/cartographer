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

proto::RangeDataInserterOptions2D CreateRangeDataInserterOptions2D(
    common::LuaParameterDictionary* const parameter_dictionary) {
  proto::RangeDataInserterOptions2D options;
  options.set_hit_probability(
      parameter_dictionary->GetDouble("hit_probability"));
  options.set_miss_probability(
      parameter_dictionary->GetDouble("miss_probability"));
  options.set_insert_free_space(
      parameter_dictionary->HasKey("insert_free_space")
          ? parameter_dictionary->GetBool("insert_free_space")
          : true);
  CHECK_GT(options.hit_probability(), 0.5);
  CHECK_LT(options.miss_probability(), 0.5);
  return options;
}

RangeDataInserter2DTSDF::RangeDataInserter2DTSDF(
    const proto::RangeDataInserterOptions2D& options)
    : options_(options) {}

void RangeDataInserter2DTSDF::Insert(const sensor::RangeData& range_data,
                                     TSDF2D* const tsdf) const {
  const float tau = 0.3;
  const float update_weight = 1.f / (2.f * tau);
  const Eigen::Array2i origin_cell =
      tsdf->limits().GetCellIndex(range_data.origin.head<2>());

  for (const Eigen::Vector3f& hit : range_data.returns) {
    const Eigen::Vector3f direction = (hit - range_data.origin).normalized();
    const Eigen::Vector3f end_position = hit + tau * direction;
    const Eigen::Array2i end_cell =
        tsdf->limits().GetCellIndex(end_position.head<2>());

    const Eigen::Array2i delta = end_cell - origin_cell;
    const int num_samples = delta.cwiseAbs().maxCoeff();
    CHECK_LT(num_samples, 1 << 15);
    // 'num_samples' is the number of samples we equi-distantly place on the
    // line between 'origin' and 'hit'. (including a fractional part for sub-
    // voxels) It is chosen so that between two samples we change from one voxel
    // to the next on the fastest changing dimension.
    //
    // Only the last 'num_free_space_voxels' are updated for performance.
    for (int position = 0; position < num_samples; ++position) {
      const Eigen::Array2i cell = origin_cell + delta * position / num_samples;
      float resolution = tsdf->limits().resolution();
      const Eigen::Vector3f cell_position = {resolution * cell[0],
                                             resolution * cell[1], 0.f};
      float distance = (hit - cell_position).dot(direction);
      if (distance > tau) {
        distance = tau;
      } else if (distance < -tau) {
        LOG(INFO) << "distance " << distance;
        distance = -tau;
      }
      tsdf->UpdateCell(cell, distance, update_weight);
    }
  }

  tsdf->FinishUpdate();
}

}  // namespace mapping
}  // namespace cartographer
