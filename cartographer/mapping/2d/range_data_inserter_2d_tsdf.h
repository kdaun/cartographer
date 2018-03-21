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

#ifndef CARTOGRAPHER_MAPPING_2D_RANGE_DATA_INSERTER_2D_TSDF_H_
#define CARTOGRAPHER_MAPPING_2D_RANGE_DATA_INSERTER_2D_TSDF_H_

#include <utility>
#include <vector>

#include "cartographer/common/lua_parameter_dictionary.h"
#include "cartographer/common/port.h"
#include "cartographer/mapping/2d/proto/range_data_inserter_options_2d.pb.h"
#include "cartographer/mapping/2d/range_data_inserter_2d_probability_grid.h" //TODO(kdaun) remove after options are split up
#include "cartographer/mapping/2d/tsdf_2d.h"
#include "cartographer/mapping/2d/xy_index.h"
#include "cartographer/sensor/point_cloud.h"
#include "cartographer/sensor/range_data.h"

namespace cartographer {
namespace mapping {

class RangeDataInserter2DTSDF {
 public:
  explicit RangeDataInserter2DTSDF(
      const proto::RangeDataInserterOptions2D& options);

  RangeDataInserter2DTSDF(const RangeDataInserter2DTSDF&) = delete;
  RangeDataInserter2DTSDF& operator=(const RangeDataInserter2DTSDF&) = delete;

  // Inserts 'range_data' into 'probability_grid'.
  void Insert(const sensor::RangeData& range_data, TSDF2D* tsdf) const;

 private:
  const proto::RangeDataInserterOptions2D options_;
};

}  // namespace mapping
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_2D_RANGE_DATA_INSERTER_2D_TSDF_H_
