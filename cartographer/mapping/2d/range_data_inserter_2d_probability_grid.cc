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

#include "cartographer/mapping/2d/range_data_inserter_2d_probability_grid.h"

#include <cstdlib>

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "cartographer/mapping/2d/xy_index.h"
#include "cartographer/mapping/internal/2d/ray_casting.h"
#include "cartographer/mapping/probability_values.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping {

RangeDataInserter2DProbabilityGrid::RangeDataInserter2DProbabilityGrid(
    const proto::RangeDataInserterOptions2D& options)
    : options_(options),
      hit_table_(ComputeLookupTableToApplyOdds(
          Odds(options_.probability_grid().hit_probability()))),
      miss_table_(ComputeLookupTableToApplyOdds(
          Odds(options_.probability_grid().miss_probability()))) {}

void RangeDataInserter2DProbabilityGrid::Insert(
    const sensor::RangeData& range_data,
    ProbabilityGrid* const probability_grid) const {
  // By not finishing the update after hits are inserted, we give hits priority
  // (i.e. no hits will be ignored because of a miss in the same cell).
  CastRays(range_data, hit_table_, miss_table_,
           options_.probability_grid().insert_free_space(),
           CHECK_NOTNULL(probability_grid));
  probability_grid->FinishUpdate();
}

}  // namespace mapping
}  // namespace cartographer
