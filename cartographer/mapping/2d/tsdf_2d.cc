#include "cartographer/mapping/2d/tsdf_2d.h"

#include <limits>

#include "cartographer/mapping/2d/proto/submap_2d.pb.h"
#include "cartographer/mapping/tsdf_values.h"

namespace cartographer {
namespace mapping {
namespace {

// Converts a 'cell_index' into an index into 'cells_'.
int ToFlatIndex(const Eigen::Array2i& cell_index, const MapLimits& limits) {
  CHECK(limits.Contains(cell_index)) << cell_index;
  return limits.cell_limits().num_x_cells * cell_index.y() + cell_index.x();
}

}  // namespace

TSDF2D::TSDF2D(const MapLimits& limits)
    : Grid2D(limits),
      tsdf_cells_(
          limits_.cell_limits().num_x_cells * limits_.cell_limits().num_y_cells,
          kUnknownTSDFValue),
      weight_cells_(
          limits_.cell_limits().num_x_cells * limits_.cell_limits().num_y_cells,
          kUnknownWeightValue),
      truncation_distance_(0.3f) {}

TSDF2D::TSDF2D(const proto::Submap2D& proto)
    : Grid2D(MapLimits(proto.limits())),
      tsdf_cells_(),
      weight_cells_(),
      truncation_distance_(0.3f) {
  LOG(ERROR) << "TSDF2D(const proto::Submap2D& proto) not implemented";
}

// Finishes the update sequence.
void TSDF2D::FinishUpdate() {
  while (!update_indices_.empty()) {
    DCHECK_GE(tsdf_cells_[update_indices_.back()], kUpdateMarker);
    tsdf_cells_[update_indices_.back()] -= kUpdateMarker;
    DCHECK_GE(weight_cells_[update_indices_.back()], kUpdateMarker);
    weight_cells_[update_indices_.back()] -= kUpdateMarker;
    update_indices_.pop_back();
  }
}

void TSDF2D::SetCell(const Eigen::Array2i& cell_index, const float tsdf,
                     const float weight) {
  uint16& cell_tsdf = tsdf_cells_[ToFlatIndex(cell_index, limits_)];
  // CHECK_EQ(cell, kUnknownProbabilityValue);
  cell_tsdf = TSDFToValue(tsdf);
  uint16& cell_weight = weight_cells_[ToFlatIndex(cell_index, limits_)];
  cell_weight = WeightToValue(weight);
  known_cells_box_.extend(cell_index.matrix());
}

bool TSDF2D::UpdateCell(const Eigen::Array2i& cell_index,
                        const float update_tsdf, const float update_weight) {
  const int flat_index = ToFlatIndex(cell_index, limits_);
  uint16* tsdf_cell = &tsdf_cells_[flat_index];
  uint16* weight_cell = &weight_cells_[flat_index];
  if (*tsdf_cell >= kUpdateMarker) {
    float current_tsdf = ValueToTSDF(*tsdf_cell - kUpdateMarker);
    float current_weight = ValueToWeight(*weight_cell);
    current_tsdf =
        (current_weight * current_tsdf + update_tsdf * update_weight) /
        (current_weight + update_weight);
    current_weight += update_weight;
    *tsdf_cell = TSDFToValue(current_tsdf) + kUpdateMarker;
    *weight_cell = WeightToValue(current_weight);
    return false;
  } else {
    update_indices_.push_back(flat_index);
    float current_tsdf = ValueToTSDF(*tsdf_cell);
    float current_weight =
        ValueToWeight(*weight_cell);  // TODO(kdaun) remove code duplication
    current_tsdf =
        (current_weight * current_tsdf + update_tsdf * update_weight) /
        (current_weight + update_weight);
    current_weight += update_weight;
    *tsdf_cell = TSDFToValue(current_tsdf) + kUpdateMarker;
    *weight_cell = WeightToValue(current_weight);
    DCHECK_GE(*tsdf_cell, kUpdateMarker);
    known_cells_box_.extend(cell_index.matrix());
    return true;
  }
}

float TSDF2D::GetTSDF(const Eigen::Array2i& cell_index) const {
  if (limits_.Contains(cell_index)) {
    return ValueToTSDF(tsdf_cells_[ToFlatIndex(cell_index, limits_)]);
  }
  return kMinTSDF;
}

float TSDF2D::GetWeight(const Eigen::Array2i& cell_index) const {
  if (limits_.Contains(cell_index)) {
    return ValueToWeight(weight_cells_[ToFlatIndex(cell_index, limits_)]);
  }
  return kMaxTSDF;
}

float TSDF2D::GetCorrespondenceCost(const Eigen::Array2i& cell_index) const {
  float correspondence_cost = GetTSDF(cell_index);
  DCHECK_GE(correspondence_cost, GetMinCorrespondenceCost());
  DCHECK_GE(std::abs(correspondence_cost), GetMinAbsCorrespondenceCost());
  DCHECK_LE(correspondence_cost, GetMaxCorrespondenceCost());
  return correspondence_cost;
}

float TSDF2D::GetMinCorrespondenceCost() const { return -truncation_distance_; }

float TSDF2D::GetMinAbsCorrespondenceCost() const { return 0.f; }

float TSDF2D::GetMaxCorrespondenceCost() const { return truncation_distance_; }

// Returns true if the probability at the specified index is known.
bool TSDF2D::IsKnown(const Eigen::Array2i& cell_index) const {
  return limits_.Contains(cell_index) &&
         tsdf_cells_[ToFlatIndex(cell_index, limits_)] != kUnknownTSDFValue;
}

proto::Submap2D TSDF2D::ToProto() const {
  proto::Submap2D result;
  *result.mutable_limits() = mapping::ToProto(limits_);
  CHECK(update_indices_.empty()) << "Serializing a grid during an update is "
      "not supported. Finish the update first.";
  if (!known_cells_box_.isEmpty()) {
    auto* const box = result.mutable_known_cells_box();
    box->set_max_x(known_cells_box_.max().x());
    box->set_max_y(known_cells_box_.max().y());
    box->set_min_x(known_cells_box_.min().x());
    box->set_min_y(known_cells_box_.min().y());
  }
  proto::Submap2DTSDFDetails details;
  details.mutable_tsdf_cells()->Reserve(tsdf_cells_.size());
  for (const auto& cell : tsdf_cells_) {
    details.mutable_tsdf_cells()->Add(cell);
  }
  details.mutable_weight_cells()->Reserve(weight_cells_.size());
  for (const auto& cell : weight_cells_) {
    details.mutable_weight_cells()->Add(cell);
  }
  result.mutable_details()->PackFrom(details);
  return result;
}

// Grows the map as necessary to include 'point'. This changes the meaning of
// these coordinates going forward. This method must be called immediately
// after 'FinishUpdate', before any calls to 'ApplyLookupTable'.
void TSDF2D::GrowLimits(const Eigen::Vector2f& point) {
  CHECK(update_indices_.empty());
  while (!limits_.Contains(limits_.GetCellIndex(point))) {
    const int x_offset = limits_.cell_limits().num_x_cells / 2;
    const int y_offset = limits_.cell_limits().num_y_cells / 2;
    const MapLimits new_limits(
        limits_.resolution(),
        limits_.max() +
            limits_.resolution() * Eigen::Vector2d(y_offset, x_offset),
        CellLimits(2 * limits_.cell_limits().num_x_cells,
                   2 * limits_.cell_limits().num_y_cells));
    const int stride = new_limits.cell_limits().num_x_cells;
    const int offset = x_offset + stride * y_offset;
    const int new_size = new_limits.cell_limits().num_x_cells *
        new_limits.cell_limits().num_y_cells;
    std::vector<uint16> new_tsdf_cells(new_size, kUnknownTSDFValue);
    std::vector<uint16> new_weight_cells(new_size, kUnknownWeightValue);
    for (int i = 0; i < limits_.cell_limits().num_y_cells; ++i) {
      for (int j = 0; j < limits_.cell_limits().num_x_cells; ++j) {
        new_tsdf_cells[offset + j + i * stride] =
            tsdf_cells_[j + i * limits_.cell_limits().num_x_cells];
        new_weight_cells[offset + j + i * stride] =
            weight_cells_[j + i * limits_.cell_limits().num_x_cells];
      }
    }
    tsdf_cells_ = new_tsdf_cells;
    weight_cells_ = new_weight_cells;
    limits_ = new_limits;
    if (!known_cells_box_.isEmpty()) {
      known_cells_box_.translate(Eigen::Vector2i(x_offset, y_offset));
    }
  }
}

}  // namespace mapping
}  // namespace cartographer
