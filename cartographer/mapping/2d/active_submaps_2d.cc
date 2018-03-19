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

#include "cartographer/mapping/2d/active_submaps_2d.h"

#include <cinttypes>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>

#include "Eigen/Geometry"
#include "cartographer/common/make_unique.h"
#include "cartographer/common/port.h"
#include "glog/logging.h"
#include "submap_2d_probability_grid.h"

namespace cartographer {
namespace mapping {

ActiveSubmaps2D::ActiveSubmaps2D(const proto::SubmapsOptions2D& options)
    : options_(options),
      range_data_inserter_probability_grid_(
          std::make_shared<RangeDataInserter2DProbabilityGrid>(
              options.range_data_inserter_options())) {
  // We always want to have at least one likelihood field which we can return,
  // and will create it at the origin in absence of a better choice.
  AddSubmap(Eigen::Vector2f::Zero());
}

void ActiveSubmaps2D::InsertRangeData(const sensor::RangeData& range_data) {
  for (auto& submap : submaps_) {
    submap->InsertRangeData(range_data);
  }
  if (submaps_.back()->num_range_data() == options_.num_range_data()) {
    AddSubmap(range_data.origin.head<2>());
  }
}

std::vector<std::shared_ptr<Submap2D>> ActiveSubmaps2D::submaps() const {
  return submaps_;
}

int ActiveSubmaps2D::matching_index() const { return matching_submap_index_; }

void ActiveSubmaps2D::FinishSubmap() {
  Submap2D* submap = submaps_.front().get();
  submap->Finish();
  ++matching_submap_index_;
  submaps_.erase(submaps_.begin());
}

void ActiveSubmaps2D::AddSubmap(const Eigen::Vector2f& origin) {
  if (submaps_.size() > 1) {
    // This will crop the finished Submap before inserting a new Submap to
    // reduce peak memory usage a bit.
    FinishSubmap();
  }
  constexpr int kInitialSubmapSize = 100;
  submaps_.push_back(common::make_unique<Submap2DProbabilityGrid>(
      MapLimits(options_.resolution(),
                origin.cast<double>() + 0.5 * kInitialSubmapSize *
                                            options_.resolution() *
                                            Eigen::Vector2d::Ones(),
                CellLimits(kInitialSubmapSize, kInitialSubmapSize)),
      origin,
      range_data_inserter_probability_grid_));  // TODO(kdaun) check
                                                // submap type

  LOG(INFO) << "Added submap " << matching_submap_index_ + submaps_.size();
}

}  // namespace mapping
}  // namespace cartographer
