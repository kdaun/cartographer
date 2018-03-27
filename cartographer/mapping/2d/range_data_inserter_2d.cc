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

#include "cartographer/mapping/2d/range_data_inserter_2d.h"

namespace cartographer {
namespace mapping {

proto::RangeDataInserterOptions2DTSDF CreateRangeDataInserterOptions2DTSDF(
    common::LuaParameterDictionary* const parameter_dictionary) {
  proto::RangeDataInserterOptions2DTSDF options;
  const std::string range_data_inserter_type_string =
      parameter_dictionary->GetString("range_data_inserter_type");
  proto::RangeDataInserter2DTSDFType range_data_inserter_type;
  CHECK(proto::RangeDataInserter2DTSDFType_Parse(
      range_data_inserter_type_string, &range_data_inserter_type))
      << "Unknown RangeDataInserter2DTSDFType kind: "
      << range_data_inserter_type_string;
  options.set_range_data_inserter_type(range_data_inserter_type);
  options.set_truncation_distance(
      parameter_dictionary->GetDouble("truncation_distance"));
  options.set_behind_surface_distance(
      parameter_dictionary->GetDouble("behind_surface_distance"));
  options.set_update_weight(parameter_dictionary->GetDouble("update_weight"));
  options.set_maximum_weight(parameter_dictionary->GetDouble("maximum_weight"));
  CHECK_GT(options.truncation_distance(), 0.0);
  CHECK_GT(options.behind_surface_distance(), 0.0);
  CHECK_GT(options.maximum_weight(), 0.0);
  CHECK_GT(options.maximum_weight(), 0.0);
  CHECK_GT(options.maximum_weight(), options.update_weight());
  return options;
}

proto::RangeDataInserterOptions2DProbabilityGrid
CreateRangeDataInserterOptions2DProbabilityGrid(
    common::LuaParameterDictionary* const parameter_dictionary) {
  proto::RangeDataInserterOptions2DProbabilityGrid options;
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

proto::RangeDataInserterOptions2D CreateRangeDataInserterOptions2D(
    common::LuaParameterDictionary* const parameter_dictionary) {
  proto::RangeDataInserterOptions2D options;
  *options.mutable_probability_grid() =
      CreateRangeDataInserterOptions2DProbabilityGrid(
          parameter_dictionary->GetDictionary("probability_grid").get());
  *options.mutable_tsdf() = CreateRangeDataInserterOptions2DTSDF(
      parameter_dictionary->GetDictionary("tsdf").get());
  return options;
}

}  // namespace mapping
}  // namespace cartographer