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
  const float tau = options_.tsdf().truncation_distance();
  const float update_weight = options_.tsdf().update_weight();

  for (const Eigen::Vector3f& hit : range_data.returns) {
    const Eigen::Vector3f direction = (hit - range_data.origin).normalized();
    const Eigen::Vector3f end_position = hit + tau * direction;
    const Eigen::Array2i end_cell =
        tsdf->limits().GetCellIndex(end_position.head<2>());
    Eigen::AlignedBox2f bounding_box(range_data.origin.head<2>());
    bounding_box.extend(end_position.head<2>());
    constexpr float kPadding = 1e-6f;
    tsdf->GrowLimits(bounding_box.min() - kPadding * Eigen::Vector2f::Ones());
    tsdf->GrowLimits(bounding_box.max() + kPadding * Eigen::Vector2f::Ones());
  }

  const Eigen::Vector2f origin = range_data.origin.head<2>();
  const Eigen::Array2i origin_cell = tsdf->limits().GetCellIndex(origin);
  for (const Eigen::Vector3f& hit_3d : range_data.returns) {
    const Eigen::Vector2f hit = hit_3d.head<2>();
    const Eigen::Vector2f direction = (hit - origin).normalized();
    const Eigen::Vector2f end_position = hit + tau * direction;
    const Eigen::Array2i end_cell = tsdf->limits().GetCellIndex(end_position);
    const Eigen::Array2i delta =
        end_cell - origin_cell;  // todo(kdaun) remove code duplication
    const int num_samples = delta.cwiseAbs().maxCoeff();
    CHECK_LT(num_samples, 1 << 15);
    // 'num_samples' is the number of samples we equi-distantly place on the
    // line between 'origin' and 'hit'. (including a fractional part for sub-
    // voxels) It is chosen so that between two samples we change from one voxel
    // to the next on the fastest changing dimension.

    for (int position = 0; position < num_samples; ++position) {
      const Eigen::Array2i cell = origin_cell + delta * position / num_samples;
      float resolution = tsdf->limits().resolution();
      const Eigen::Vector2f cell_position = tsdf->limits().GetCellCenter(cell);
      float distance = (hit - cell_position).dot(direction);
      if (distance > tau) {
        distance = tau;
      } else if (distance < -tau) {
        LOG(INFO) << "\n\ndirection:\n " << direction
                  << "\nrange_data.origin:\n " << range_data.origin
                  << "\nhit:\n " << hit << "\nend_position:\n " << end_position
                  << "\ncell_position:\n " << cell_position << "\ndistance:\n "
                  << distance << "\norigin_cell:\n " << origin_cell
                  << "\nend_cell:\n " << end_cell << "\ncurrent_cell:\n "
                  << cell;
        distance = -tau;
      }
      tsdf->UpdateCell(cell, distance, update_weight);
    }
  }

  tsdf->FinishUpdate();
}

}  // namespace mapping
}  // namespace cartographer
