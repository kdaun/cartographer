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

#include "cartographer/mapping/2d/submap_2d.h"

#include <cinttypes>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>

#include "Eigen/Geometry"
#include "cartographer/common/make_unique.h"
#include "cartographer/common/port.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping {

proto::SubmapsOptions2D CreateSubmapsOptions2D(
    common::LuaParameterDictionary* const parameter_dictionary) {
  proto::SubmapsOptions2D options;
  options.set_resolution(parameter_dictionary->GetDouble("resolution"));
  const std::string map_type_string = parameter_dictionary->GetString("map_type");
  proto::MapType map_type;
  CHECK(proto::MapType_Parse(map_type_string, &map_type))
  << "Unknown MapType kind: " << map_type_string;
  options.set_map_type(map_type);
  options.set_num_range_data(
      parameter_dictionary->GetNonNegativeInt("num_range_data"));
  *options.mutable_range_data_inserter_options() =
      CreateRangeDataInserterOptions2D(
          parameter_dictionary->GetDictionary("range_data_inserter").get());
  CHECK_GT(options.num_range_data(), 0);
  return options;
}

Submap2D::Submap2D(const MapLimits& limits, const Eigen::Vector2f& origin)
    : Submap(transform::Rigid3d::Translation(
          Eigen::Vector3d(origin.x(), origin.y(), 0.))) {}

Submap2D::Submap2D(const proto::Submap& proto)
    : Submap(transform::ToRigid3(proto.local_pose())) {
  SetNumRangeData(proto.num_range_data());
  SetFinished(proto.finished());
}

}  // namespace mapping
}  // namespace cartographer
