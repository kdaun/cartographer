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

#include <memory>

#include "cartographer/common/lua_parameter_dictionary.h"
#include "cartographer/common/lua_parameter_dictionary_test_helpers.h"
#include "cartographer/common/make_unique.h"
#include "cartographer/mapping/2d/probability_grid.h"
#include "cartographer/mapping/2d/range_data_inserter_2d.h"
#include "cartographer/mapping/probability_values.h"
#include "gmock/gmock.h"

namespace cartographer {
namespace mapping {
namespace {

class RangeDataInserterTest2DTSDF : public ::testing::Test {
 protected:
  RangeDataInserterTest2DTSDF()
      : tsdf_(MapLimits(1., Eigen::Vector2d(0., 7.), CellLimits(8, 1)), 2.0,
              10.0) {
    auto parameter_dictionary = common::MakeDictionary(
        "return { "
        "probability_grid = {"
        "insert_free_space = true, "
        "hit_probability = 0.7, "
        "miss_probability = 0.4, "
        "},"
        "tsdf = {"
        "range_data_inserter_type = \"CONSTANT_WEIGHT\","
        "truncation_distance = 2.0,"
        "behind_surface_distance = 2.0,"
        "update_weight = 1.0,"
        "maximum_weight = 10.,"
        "},"
        "}");
    options_ = CreateRangeDataInserterOptions2D(parameter_dictionary.get());
    range_data_inserter_ =
        common::make_unique<RangeDataInserter2DTSDF>(options_);
  }

  void InsertPoint() {
    sensor::RangeData range_data;
    range_data.returns.emplace_back(-0.5f, 3.5f, 0.f);
    range_data.origin.x() = -0.5f;
    range_data.origin.y() = -0.5f;
    range_data_inserter_->Insert(range_data, &tsdf_);
    tsdf_.FinishUpdate();
  }

  proto::RangeDataInserterOptions2D options_;
  TSDF2D tsdf_;
  std::unique_ptr<RangeDataInserter2DTSDF> range_data_inserter_;
};

TEST_F(RangeDataInserterTest2DTSDF, InsertPointConstantIntegrator) {
  InsertPoint();

  for (float y = -0.5; y < 5.; ++y) {
    float x = -0.5f;
    Eigen::Array2i cell_index =
        tsdf_.limits().GetCellIndex(Eigen::Vector2f(x, y));
    float expected_tsdf =
        std::max(std::min(3.5f - y, tsdf_.GetMaxTSDF()), tsdf_.GetMinTSDF());
    float expected_weight = options_.tsdf().update_weight();
    EXPECT_TRUE(tsdf_.IsKnown(cell_index));
    EXPECT_NEAR(expected_tsdf, tsdf_.GetTSDF(cell_index), 1e-4);
    EXPECT_NEAR(expected_weight, tsdf_.GetWeight(cell_index), 1e-3);
    x = 0.5f;
    cell_index = tsdf_.limits().GetCellIndex(Eigen::Vector2f(x, y));
    expected_tsdf = tsdf_.GetMaxTSDF();
    expected_weight = 0.;
    EXPECT_FALSE(tsdf_.IsKnown(cell_index));
    EXPECT_NEAR(expected_tsdf, tsdf_.GetTSDF(cell_index), 1e-4);
    EXPECT_NEAR(expected_weight, tsdf_.GetWeight(cell_index), 1e-3);
    x = 1.5f;
    cell_index = tsdf_.limits().GetCellIndex(Eigen::Vector2f(x, y));
    EXPECT_FALSE(tsdf_.IsKnown(cell_index));
    EXPECT_NEAR(expected_tsdf, tsdf_.GetTSDF(cell_index), 1e-4);
    EXPECT_NEAR(expected_weight, tsdf_.GetWeight(cell_index), 1e-3);
  }

  Eigen::Array2i cell_index =
      tsdf_.limits().GetCellIndex(Eigen::Vector2f(-0.5, 5.5));
  EXPECT_FALSE(tsdf_.IsKnown(cell_index));
  EXPECT_NEAR(tsdf_.GetMaxTSDF(), tsdf_.GetTSDF(cell_index), 1e-4);
  EXPECT_NEAR(0., tsdf_.GetWeight(cell_index), 1e-3);
  cell_index = tsdf_.limits().GetCellIndex(Eigen::Vector2f(-0.5, -1.5));
  EXPECT_FALSE(tsdf_.IsKnown(cell_index));
  EXPECT_NEAR(tsdf_.GetMaxTSDF(), tsdf_.GetTSDF(cell_index), 1e-4);
  EXPECT_NEAR(0., tsdf_.GetWeight(cell_index), 1e-3);

  for (int i = 0; i < 1000; ++i) {
    InsertPoint();
  }
  for (float y = -0.5; y < 5.; ++y) {
    float x = -0.5f;
    Eigen::Array2i cell_index =
        tsdf_.limits().GetCellIndex(Eigen::Vector2f(x, y));
    float expected_tsdf =
        std::max(std::min(3.5f - y, tsdf_.GetMaxTSDF()), tsdf_.GetMinTSDF());
    EXPECT_TRUE(tsdf_.IsKnown(cell_index));
    EXPECT_NEAR(expected_tsdf, tsdf_.GetTSDF(cell_index), 1e-4);
    EXPECT_NEAR(tsdf_.GetMaxWeight(), tsdf_.GetWeight(cell_index), 1e-3);
  }
}

TEST_F(RangeDataInserterTest2DTSDF, InsertPointLinearIntegrator) {
  options_.mutable_tsdf()->set_range_data_inserter_type(proto::LINEAR_WEIGHT);
  range_data_inserter_ = common::make_unique<RangeDataInserter2DTSDF>(options_);
  InsertPoint();
  for (float y = -0.5; y < 5.; ++y) {
    float x = -0.5f;
    Eigen::Array2i cell_index =
        tsdf_.limits().GetCellIndex(Eigen::Vector2f(x, y));
    float expected_tsdf =
        std::max(std::min(3.5f - y, tsdf_.GetMaxTSDF()), tsdf_.GetMinTSDF());
    float expected_weight = options_.tsdf().update_weight() / 4.f;
    EXPECT_TRUE(tsdf_.IsKnown(cell_index));
    EXPECT_NEAR(expected_tsdf, tsdf_.GetTSDF(cell_index), 1e-4);
    EXPECT_NEAR(expected_weight, tsdf_.GetWeight(cell_index), 1e-3);
  }

  for (int i = 0; i < 1000; ++i) {
    InsertPoint();
  }
  for (float y = -0.5; y < 5.; ++y) {
    float x = -0.5f;
    Eigen::Array2i cell_index =
        tsdf_.limits().GetCellIndex(Eigen::Vector2f(x, y));
    float expected_tsdf =
        std::max(std::min(3.5f - y, tsdf_.GetMaxTSDF()), tsdf_.GetMinTSDF());
    EXPECT_TRUE(tsdf_.IsKnown(cell_index));
    EXPECT_NEAR(expected_tsdf, tsdf_.GetTSDF(cell_index), 1e-4);
    EXPECT_NEAR(tsdf_.GetMaxWeight(), tsdf_.GetWeight(cell_index), 1e-3);
  }
}

TEST_F(RangeDataInserterTest2DTSDF, InsertPointQuadraticIntegrator) {
  options_.mutable_tsdf()->set_range_data_inserter_type(
      proto::QUADRATIC_WEIGHT);
  range_data_inserter_ = common::make_unique<RangeDataInserter2DTSDF>(options_);
  InsertPoint();
  for (float y = -0.5; y < 5.; ++y) {
    float x = -0.5f;
    Eigen::Array2i cell_index =
        tsdf_.limits().GetCellIndex(Eigen::Vector2f(x, y));
    float expected_tsdf =
        std::max(std::min(3.5f - y, tsdf_.GetMaxTSDF()), tsdf_.GetMinTSDF());
    float expected_weight = options_.tsdf().update_weight() / 16.f;
    EXPECT_TRUE(tsdf_.IsKnown(cell_index));
    EXPECT_NEAR(expected_tsdf, tsdf_.GetTSDF(cell_index), 1e-4);
    EXPECT_NEAR(expected_weight, tsdf_.GetWeight(cell_index), 1e-3);
  }

  for (int i = 0; i < 1000; ++i) {
    InsertPoint();
  }
  for (float y = -0.5; y < 5.; ++y) {
    float x = -0.5f;
    Eigen::Array2i cell_index =
        tsdf_.limits().GetCellIndex(Eigen::Vector2f(x, y));
    float expected_tsdf =
        std::max(std::min(3.5f - y, tsdf_.GetMaxTSDF()), tsdf_.GetMinTSDF());
    EXPECT_TRUE(tsdf_.IsKnown(cell_index));
    EXPECT_NEAR(expected_tsdf, tsdf_.GetTSDF(cell_index), 1e-4);
    EXPECT_NEAR(tsdf_.GetMaxWeight(), tsdf_.GetWeight(cell_index), 1e-3);
  }
}

}  // namespace
}  // namespace mapping
}  // namespace cartographer
