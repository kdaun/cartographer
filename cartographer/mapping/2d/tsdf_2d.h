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

#ifndef CARTOGRAPHER_MAPPING_2D_TSDF_2D_H_
#define CARTOGRAPHER_MAPPING_2D_TSDF_2D_H_

#include <vector>

#include "cartographer/common/port.h"
#include "cartographer/mapping/2d/grid_2d.h"
#include "cartographer/mapping/2d/map_limits.h"
#include "cartographer/mapping/2d/proto/submap_2d.pb.h"
#include "cartographer/mapping/2d/xy_index.h"
#include "cartographer/mapping/tsdf_values.h"

namespace cartographer {
namespace mapping {

// Represents a 2D grid of probabilities.
class TSDF2D : public Grid2D {
 public:
  explicit TSDF2D(const MapLimits& limits, float truncation_distance,
                  float max_weight);
  explicit TSDF2D(const proto::Submap2D& proto);

  void FinishUpdate() override;

  void SetCell(const Eigen::Array2i& cell_index, const float tsdf,
               const float weight);

  bool UpdateCell(const Eigen::Array2i& cell_index, const float updated_sdf,
                  const float updated_weight = 0.1);

  float GetTSDF(const Eigen::Array2i& cell_index) const;

  float GetWeight(const Eigen::Array2i& cell_index) const;

  float GetCorrespondenceCost(const Eigen::Array2i& cell_index) const override;
  float GetMinCorrespondenceCost() const override;
  float GetMinAbsCorrespondenceCost() const override;
  float GetMaxCorrespondenceCost() const override;
  float GetMaxTSDF() const;
  float GetMinTSDF() const;
  float GetMaxWeight() const;
  float GetMinWeight() const;

  bool IsKnown(const Eigen::Array2i& cell_index) const override;

  virtual void GrowLimits(const Eigen::Vector2f& point) override;

  proto::Submap2D ToProto() const;

 private:
  float truncation_distance_;
  float max_weight_;
  TSDFValueHelper value_helper;
  std::vector<uint16> tsdf_cells_;    // Highest bit is update marker.
  std::vector<uint16> weight_cells_;  // Highest bit is update marker.
};

}  // namespace mapping
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_2D_TSDF_2D_H_
