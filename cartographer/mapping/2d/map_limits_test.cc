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

#include "cartographer/mapping/2d/map_limits.h"

#include "gtest/gtest.h"

namespace cartographer {
namespace mapping {
namespace {

TEST(MapLimitsTest, ToProto) {
  const MapLimits limits(42., Eigen::Vector2d(3., 0.), CellLimits(2, 3));
  const auto proto = ToProto(limits);
  EXPECT_EQ(limits.resolution(), proto.resolution());
  EXPECT_EQ(limits.max().x(), proto.max().x());
  EXPECT_EQ(limits.max().y(), proto.max().y());
  EXPECT_EQ(ToProto(limits.cell_limits()).DebugString(),
            proto.cell_limits().DebugString());
}

TEST(MapLimitsTest, ProtoConstructor) {
  proto::MapLimits limits;
  limits.set_resolution(1.);
  limits.mutable_max()->set_x(2.);
  limits.mutable_max()->set_y(3.);
  limits.mutable_cell_limits()->set_num_x_cells(4);
  limits.mutable_cell_limits()->set_num_y_cells(5);

  const MapLimits native(limits);
  EXPECT_EQ(limits.resolution(), native.resolution());
  EXPECT_EQ(limits.max().x(), native.max().x());
  EXPECT_EQ(limits.max().y(), native.max().y());
  EXPECT_EQ(limits.cell_limits().DebugString(),
            ToProto(native.cell_limits()).DebugString());
}

TEST(MapLimitsTest, ConstructAndGet) {
  const MapLimits limits(42., Eigen::Vector2d(3., 0.), CellLimits(2, 3));

  Eigen::Array2i cell_idx = limits.GetCellIndex({-1.f, -1.f});
  EXPECT_NEAR(limits.GetCellCenter(cell_idx)[0], -18.f, 1e-5);
  EXPECT_NEAR(limits.GetCellCenter(cell_idx)[1], -21.f, 1e-5);

  cell_idx = limits.GetCellIndex({2.9f, -0.1f});
  EXPECT_NEAR(limits.GetCellCenter(cell_idx)[0], -18.f, 1e-5);
  EXPECT_NEAR(limits.GetCellCenter(cell_idx)[1], -21.f, 1e-5);

  cell_idx = limits.GetCellIndex({-38.9f, -41.9f});
  EXPECT_NEAR(limits.GetCellCenter(cell_idx)[0], -18.f, 1e-5);
  EXPECT_NEAR(limits.GetCellCenter(cell_idx)[1], -21.f, 1e-5);

  cell_idx = limits.GetCellIndex({-50.f, -0.1f});
  EXPECT_NEAR(limits.GetCellCenter(cell_idx)[0], -60.f, 1e-5);
  EXPECT_NEAR(limits.GetCellCenter(cell_idx)[1], -21.f, 1e-5);

  cell_idx = limits.GetCellIndex({2.9f, -50.f});
  EXPECT_NEAR(limits.GetCellCenter(cell_idx)[0], -18.f, 1e-5);
  EXPECT_NEAR(limits.GetCellCenter(cell_idx)[1], -63.f, 1e-5);
}

}  // namespace
}  // namespace mapping
}  // namespace cartographer
