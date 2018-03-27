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

#ifndef CARTOGRAPHER_MAPPING_2D_ACTIVE_SUBMAPS_2D_H_
#define CARTOGRAPHER_MAPPING_2D_ACTIVE_SUBMAPS_2D_H_

#include "Eigen/Core"
#include "cartographer/common/lua_parameter_dictionary.h"
#include "cartographer/mapping/2d/grid_2d.h"
#include "cartographer/mapping/2d/map_limits.h"
#include "cartographer/mapping/2d/proto/submaps_options_2d.pb.h"
#include "cartographer/mapping/2d/range_data_inserter_2d_probability_grid.h"
#include "cartographer/mapping/2d/range_data_inserter_2d_tsdf.h"
#include "cartographer/mapping/2d/submap_2d.h"
#include "cartographer/mapping/proto/serialization.pb.h"
#include "cartographer/mapping/proto/submap_visualization.pb.h"
#include "cartographer/mapping/submaps.h"
#include "cartographer/mapping/trajectory_node.h"
#include "cartographer/sensor/range_data.h"
#include "cartographer/transform/rigid_transform.h"

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
#include "submap_2d_tsdf.h"

namespace cartographer {
namespace mapping {

// Except during initialization when only a single submap exists, there are
// always two submaps into which range data is inserted: an old submap that is
// used for matching, and a new one, which will be used for matching next, that
// is being initialized.
//
// Once a certain number of range data have been inserted, the new submap is
// considered initialized: the old submap is no longer changed, the "new" submap
// is now the "old" submap and is used for scan-to-map matching. Moreover, a
// "new" submap gets created. The "old" submap is forgotten by this object.
class ActiveSubmaps2D {
 public:
  // explicit ActiveSubmaps2D(const proto::SubmapsOptions2D& options);

  // ActiveSubmaps2D(const ActiveSubmaps2D&) = delete;
  ActiveSubmaps2D& operator=(const ActiveSubmaps2D&) = delete;

  // Returns the index of the newest initialized Submap which can be
  // used for scan-to-map matching.
  virtual int matching_index() const = 0;

  // Inserts 'range_data' into the Submap collection.
  virtual void InsertRangeData(const sensor::RangeData& range_data) = 0;

  virtual std::vector<std::shared_ptr<Submap2D>> submaps() const = 0;
};

template <typename Submap2DType, typename RangeDataInserter2DType>
class ActiveSubmaps2DI : public ActiveSubmaps2D {
 public:
  explicit ActiveSubmaps2DI(const proto::SubmapsOptions2D& options)
      : options_(options),
        range_data_inserter_(options.range_data_inserter_options()) {
    // We always want to have at least one likelihood field which we can return,
    // and will create it at the origin in absence of a better choice.
    AddSubmap(Eigen::Vector2f::Zero());
  }

  ActiveSubmaps2DI(const ActiveSubmaps2DI&) = delete;
  ActiveSubmaps2DI& operator=(const ActiveSubmaps2DI&) = delete;
  // Returns the index of the newest initialized Submap which can be
  // used for scan-to-map matching.
  int matching_index() const { return matching_submap_index_; };

  // Inserts 'range_data' into the Submap collection.
  void InsertRangeData(const sensor::RangeData& range_data) {
    for (auto& submap : submaps_) {
      submap->InsertRangeData(range_data, range_data_inserter_);
    }
    if (submaps_.back()->num_range_data() == options_.num_range_data()) {
      AddSubmap(range_data.origin.head<2>());
    }
  };

  std::vector<std::shared_ptr<Submap2D>> submaps() const {
    std::vector<std::shared_ptr<Submap2D>> result(submaps_.begin(),
                                                  submaps_.end());
    return result;
  }

 private:
  void FinishSubmap() {
    Submap2D* submap = submaps_.front().get();
    submap->Finish();
    ++matching_submap_index_;
    submaps_.erase(submaps_.begin());
  }

  void AddSubmap(const Eigen::Vector2f& origin) {
    if (submaps_.size() > 1) {
      // This will crop the finished Submap before inserting a new Submap to
      // reduce peak memory usage a bit.
      FinishSubmap();
    }
    constexpr int kInitialSubmapSize = 100;

    submaps_.push_back(common::make_unique<Submap2DType>(
        MapLimits(options_.resolution(),
                  origin.cast<double>() + 0.5 * kInitialSubmapSize *
                                              options_.resolution() *
                                              Eigen::Vector2d::Ones(),
                  CellLimits(kInitialSubmapSize, kInitialSubmapSize)),
        origin));

    LOG(INFO) << "Added submap " << matching_submap_index_ + submaps_.size();
  }

  std::vector<std::shared_ptr<Submap2DType>> submaps_;
  const proto::SubmapsOptions2D options_;
  RangeDataInserter2DType range_data_inserter_;
  int matching_submap_index_ = 0;
};

template <>
inline void ActiveSubmaps2DI<Submap2DTSDF, RangeDataInserter2DTSDF>::AddSubmap(
    const Eigen::Vector2f& origin) {
  if (submaps_.size() > 1) {
    // This will crop the finished Submap before inserting a new Submap to
    // reduce peak memory usage a bit.
    FinishSubmap();
  }
  constexpr int kInitialSubmapSize = 100;

  const proto::RangeDataInserterOptions2DTSDF& tsdf_options =
      options_.range_data_inserter_options().tsdf();
  submaps_.push_back(common::make_unique<Submap2DTSDF>(
      MapLimits(options_.resolution(),
                origin.cast<double>() + 0.5 * kInitialSubmapSize *
                                            options_.resolution() *
                                            Eigen::Vector2d::Ones(),
                CellLimits(kInitialSubmapSize, kInitialSubmapSize)),
      origin, tsdf_options.truncation_distance(),
      tsdf_options.maximum_weight()));

  LOG(INFO) << "Added submap " << matching_submap_index_ + submaps_.size();
}

}  // namespace mapping
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_2D_ACTIVE_SUBMAPS_2D_H_
