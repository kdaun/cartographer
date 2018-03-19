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

#include "cartographer/mapping/2d/submap_2d_probability_grid.h"

#include <cinttypes>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>

#include "Eigen/Geometry"
#include "cartographer/common/make_unique.h"
#include "cartographer/common/port.h"
#include "cartographer/mapping/2d/probability_grid.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping {

ProbabilityGrid ComputeCroppedProbabilityGrid(
    const ProbabilityGrid& probability_grid) {
  Eigen::Array2i offset;
  CellLimits limits;
  probability_grid.ComputeCroppedLimits(&offset, &limits);
  const double resolution = probability_grid.limits().resolution();
  const Eigen::Vector2d max =
      probability_grid.limits().max() -
      resolution * Eigen::Vector2d(offset.y(), offset.x());
  ProbabilityGrid cropped_grid(MapLimits(resolution, max, limits));
  for (const Eigen::Array2i& xy_index : XYIndexRangeIterator(limits)) {
    if (probability_grid.IsKnown(xy_index + offset)) {
      cropped_grid.SetProbability(
          xy_index, probability_grid.GetProbability(xy_index + offset));
    }
  }
  return cropped_grid;
}

Submap2DProbabilityGrid::Submap2DProbabilityGrid(
    const MapLimits& limits, const Eigen::Vector2f& origin,
    std::shared_ptr<RangeDataInserter2DProbabilityGrid> range_data_inserter)
    : Submap2D(limits, origin),
      probability_grid_(limits),
      range_data_inserter_(range_data_inserter) {}

Submap2DProbabilityGrid::Submap2DProbabilityGrid(
    const proto::Submap&
        proto)  // todo(kdaun) discuss how to handle range_data_inserter here
    : Submap2D(proto), probability_grid_(ProbabilityGrid(proto.submap_2d())) {}

void Submap2DProbabilityGrid::ToProto(
    proto::Submap* const proto, bool include_probability_grid_data) const {
  auto* const local_pose = proto->mutable_local_pose();
  *local_pose = transform::ToProto(this->local_pose());
  proto->set_num_range_data(num_range_data());
  proto->set_finished(finished());
  auto* const submap_2d = proto->mutable_submap_2d();
  if (include_probability_grid_data) {
    *submap_2d = probability_grid_.ToProto();
  }
}

void Submap2DProbabilityGrid::UpdateFromProto(const proto::Submap& proto) {
  CHECK(proto.has_submap_2d());
  const auto& submap_2d = proto.submap_2d();
  SetNumRangeData(proto.num_range_data());
  SetFinished(proto.finished());
  if (submap_2d.has_details()) {
    CHECK(submap_2d.details().Is<proto::Submap2DProbabilityGridDetails>());
    probability_grid_ = ProbabilityGrid(submap_2d);
  }
}

void Submap2DProbabilityGrid::ToResponseProto(
    const transform::Rigid3d&,
    proto::SubmapQuery::Response* const response) const {
  response->set_submap_version(num_range_data());

  Eigen::Array2i offset;
  CellLimits limits;
  probability_grid_.ComputeCroppedLimits(&offset, &limits);

  std::string cells;
  for (const Eigen::Array2i& xy_index : XYIndexRangeIterator(limits)) {
    if (probability_grid_.IsKnown(xy_index + offset)) {
      // We would like to add 'delta' but this is not possible using a value and
      // alpha. We use premultiplied alpha, so when 'delta' is positive we can
      // add it by setting 'alpha' to zero. If it is negative, we set 'value' to
      // zero, and use 'alpha' to subtract. This is only correct when the pixel
      // is currently white, so walls will look too gray. This should be hard to
      // detect visually for the user, though.
      const int delta =
          128 - ProbabilityToLogOddsInteger(
                    probability_grid_.GetProbability(xy_index + offset));
      const uint8 alpha = delta > 0 ? 0 : -delta;
      const uint8 value = delta > 0 ? delta : 0;
      cells.push_back(value);
      cells.push_back((value || alpha) ? alpha : 1);
    } else {
      constexpr uint8 kUnknownLogOdds = 0;
      cells.push_back(static_cast<uint8>(kUnknownLogOdds));  // value
      cells.push_back(0);                                    // alpha
    }
  }
  proto::SubmapQuery::Response::SubmapTexture* const texture =
      response->add_textures();
  common::FastGzipString(cells, texture->mutable_cells());

  texture->set_width(limits.num_x_cells);
  texture->set_height(limits.num_y_cells);
  const double resolution = probability_grid_.limits().resolution();
  texture->set_resolution(resolution);
  const double max_x =
      probability_grid_.limits().max().x() - resolution * offset.y();
  const double max_y =
      probability_grid_.limits().max().y() - resolution * offset.x();
  *texture->mutable_slice_pose() = transform::ToProto(
      local_pose().inverse() *
      transform::Rigid3d::Translation(Eigen::Vector3d(max_x, max_y, 0.)));
}

void Submap2DProbabilityGrid::InsertRangeData(
    const sensor::RangeData& range_data) {
  CHECK(!finished());
  range_data_inserter_->Insert(range_data, &probability_grid_);
  SetNumRangeData(num_range_data() + 1);
}

void Submap2DProbabilityGrid::Finish() {
  CHECK(!finished());
  probability_grid_ = ComputeCroppedProbabilityGrid(probability_grid_);
  SetFinished(true);
}
}  // namespace mapping
}  // namespace cartographer