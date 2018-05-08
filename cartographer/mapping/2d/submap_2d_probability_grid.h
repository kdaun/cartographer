/*
 * Copyright 2016 The Cartographer Authors
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

#ifndef CARTOGRAPHER_MAPPING_2D_SUBMAP_2D_PROBABILITY_GRID_H_
#define CARTOGRAPHER_MAPPING_2D_SUBMAP_2D_PROBABILITY_GRID_H_

#include <memory>
#include <vector>

#include "Eigen/Core"
#include "cartographer/common/lua_parameter_dictionary.h"
#include "cartographer/mapping/2d/map_limits.h"
#include "cartographer/mapping/2d/probability_grid.h"
#include "cartographer/mapping/2d/proto/submaps_options_2d.pb.h"
#include "cartographer/mapping/2d/range_data_inserter_2d_probability_grid.h"
#include "cartographer/mapping/2d/range_data_inserter_2d_tsdf.h"
#include "cartographer/mapping/2d/submap_2d.h"
#include "cartographer/mapping/2d/tsdf_2d.h"
#include "cartographer/mapping/proto/serialization.pb.h"
#include "cartographer/mapping/proto/submap_visualization.pb.h"
#include "cartographer/mapping/submaps.h"
#include "cartographer/mapping/trajectory_node.h"
#include "cartographer/sensor/range_data.h"
#include "cartographer/transform/rigid_transform.h"

namespace cartographer {
namespace mapping {

ProbabilityGrid ComputeCroppedProbabilityGrid(
    const ProbabilityGrid& probability_grid);

class Submap2DProbabilityGrid : public Submap2D {
 public:
  Submap2DProbabilityGrid(const MapLimits& limits,
                          const Eigen::Vector2f& origin);
  explicit Submap2DProbabilityGrid(const proto::Submap& proto);

  void ToProto(proto::Submap* proto,
               bool include_probability_grid_data) const override;
  void UpdateFromProto(const proto::Submap& proto) override;

  void ToResponseProto(const transform::Rigid3d& global_submap_pose,
                       proto::SubmapQuery::Response* response) const override;

  const ProbabilityGrid& probability_grid() const override {
    return probability_grid_;
  }
  const TSDF2D& tsdf() const override { return tsdf_; }

  // Insert 'range_data' into this submap using 'range_data_inserter'. The
  // submap must not be finished yet.
  void InsertRangeData(
      const sensor::RangeData& range_data,
      const RangeDataInserter2DProbabilityGrid& range_data_inserter);
  void InsertRangeData(const sensor::RangeData& range_data,
                       const RangeDataInserter2DTSDF& range_data_inserter);
  void Finish() override;

 private:
  ProbabilityGrid probability_grid_;
  TSDF2D tsdf_;
};

}  // namespace mapping
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_2D_SUBMAP_2D_PROBABILITY_GRID_H_
