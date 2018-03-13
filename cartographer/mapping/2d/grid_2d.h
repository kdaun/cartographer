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

#ifndef CARTOGRAPHER_MAPPING_2D_GRID_2D_H_
#define CARTOGRAPHER_MAPPING_2D_GRID_2D_H_

#include <vector>

#include "cartographer/common/port.h"
#include "cartographer/mapping/2d/map_limits.h"
#include "cartographer/mapping/2d/proto/probability_grid.pb.h"
#include "cartographer/mapping/2d/xy_index.h"

namespace cartographer {
namespace mapping {

// Represents a 2D grid of probabilities.
class Grid2D {
 public:
  explicit Grid2D(const MapLimits& limits) : limits_(limits) {}
  virtual ~Grid2D() {}

  // Returns the limits of this Grid2D.
  const MapLimits& limits() const { return limits_; }

  // Finishes the update sequence.
  virtual void FinishUpdate() = 0;

  // Returns the probability of the cell with 'cell_index'.
  virtual float GetCorrespondance(const Eigen::Array2i& cell_index) const = 0;

  // Returns true if the probability at the specified index is known.
  virtual bool IsKnown(const Eigen::Array2i& cell_index) const = 0;

  // Fills in 'offset' and 'limits' to define a subregion of that contains all
  // known cells.
  void ComputeCroppedLimits(Eigen::Array2i* const offset,
                            CellLimits* const limits) const;

  // Grows the map as necessary to include 'point'. This changes the meaning of
  // these coordinates going forward. This method must be called immediately
  // after 'FinishUpdate', before any calls to 'ApplyLookupTable'.
  virtual void GrowLimits(const Eigen::Vector2f& point) = 0;

 protected:
  MapLimits limits_;
  std::vector<int> update_indices_;

  // Bounding box of known cells to efficiently compute cropping limits.
  Eigen::AlignedBox2i known_cells_box_;
};

}  // namespace mapping
}  // namespace cartographer

#endif  // CARTOGRAPHER_MAPPING_2D_GRID_2D_H_
