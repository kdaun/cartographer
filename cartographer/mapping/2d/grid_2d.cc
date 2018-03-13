#include "cartographer/mapping/2d/grid_2d.h"

#include <limits>

#include "cartographer/mapping/probability_values.h"

namespace cartographer {
namespace mapping {

// Fills in 'offset' and 'limits' to define a subregion of that contains all
// known cells.
void Grid2D::ComputeCroppedLimits(Eigen::Array2i* const offset,
                                  CellLimits* const limits) const {
  if (known_cells_box_.isEmpty()) {
    *offset = Eigen::Array2i::Zero();
    *limits = CellLimits(1, 1);
    return;
  }
  *offset = known_cells_box_.min().array();
  *limits = CellLimits(known_cells_box_.sizes().x() + 1,
                       known_cells_box_.sizes().y() + 1);
}

}  // namespace mapping
}  // namespace cartographer
