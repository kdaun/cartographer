#include "cartographer/mapping/2d/tsdf_2d.h"

#include <limits>

#include "cartographer/common/make_unique.h"
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

TSDF2D::TSDF2D(const MapLimits& limits, float truncation_distance,
               float max_weight)
    : Grid2D(limits),
      truncation_distance_(truncation_distance),
      max_weight_(max_weight),
      value_helper(common::make_unique<TSDFValueHelper>(truncation_distance,
                                                        max_weight)),
      tsdf_cells_(
          limits_.cell_limits().num_x_cells * limits_.cell_limits().num_y_cells,
          value_helper->getUnknownTSDFValue()),
      weight_cells_(
          limits_.cell_limits().num_x_cells * limits_.cell_limits().num_y_cells,
          value_helper->getUnknownWeightValue()) {}

TSDF2D::TSDF2D(const proto::Submap2D& proto)
    : Grid2D(MapLimits(proto.limits())) {
  if (proto.has_known_cells_box()) {
    const auto& box = proto.known_cells_box();
    known_cells_box_ =
        Eigen::AlignedBox2i(Eigen::Vector2i(box.min_x(), box.min_y()),
                            Eigen::Vector2i(box.max_x(), box.max_y()));
  }
  proto::Submap2DTSDFDetails details;
  CHECK(proto.details().Is<proto::Submap2DTSDFDetails>());
  proto.details().UnpackTo(&details);
  truncation_distance_ = details.truncation_distance();
  max_weight_ = details.maximum_weight();
  value_helper =
      common::make_unique<TSDFValueHelper>(truncation_distance_, max_weight_);
  tsdf_cells_.reserve(details.tsdf_cells_size());
  for (const auto& cell : details.tsdf_cells()) {
    CHECK_LE(cell, std::numeric_limits<uint16>::max());
    tsdf_cells_.push_back(cell);
  }
  weight_cells_.reserve(details.tsdf_cells_size());
  for (const auto& cell : details.weight_cells()) {
    CHECK_LE(cell, std::numeric_limits<uint16>::max());
    weight_cells_.push_back(cell);
  }
}

// Finishes the update sequence.
void TSDF2D::FinishUpdate() {
  while (!update_indices_.empty()) {
    DCHECK_GE(tsdf_cells_[update_indices_.back()],
              value_helper->getUpdateMarker());
    tsdf_cells_[update_indices_.back()] -= value_helper->getUpdateMarker();
    update_indices_.pop_back();
  }
}

void TSDF2D::SetCell(const Eigen::Array2i& cell_index, const float tsdf,
                     const float weight) {
  uint16& cell_tsdf = tsdf_cells_[ToFlatIndex(cell_index, limits_)];
  CHECK_EQ(cell_tsdf, value_helper->getUnknownTSDFValue());
  cell_tsdf = value_helper->TSDFToValue(tsdf);
  uint16& cell_weight = weight_cells_[ToFlatIndex(cell_index, limits_)];
  CHECK_EQ(cell_weight, value_helper->getUnknownWeightValue());
  cell_weight = value_helper->WeightToValue(weight);
  known_cells_box_.extend(cell_index.matrix());
}

bool TSDF2D::UpdateCell(const Eigen::Array2i& cell_index,
                        const float updated_sdf, const float updated_weight) {
  const int flat_index = ToFlatIndex(cell_index, limits_);
  uint16* tsdf_cell = &tsdf_cells_[flat_index];
  uint16* weight_cell = &weight_cells_[flat_index];
  if (*tsdf_cell < value_helper->getUpdateMarker()) {
    update_indices_.push_back(flat_index);
    known_cells_box_.extend(cell_index.matrix());
    *tsdf_cell =
        value_helper->TSDFToValue(updated_sdf) + value_helper->getUpdateMarker();
    *weight_cell = value_helper->WeightToValue(updated_weight);
  }
  return true;
}

float TSDF2D::GetTSDF(const Eigen::Array2i& cell_index) const {
  if (limits_.Contains(cell_index)) {
    return value_helper->ValueToTSDF(
        tsdf_cells_[ToFlatIndex(cell_index, limits_)]);
  }
  return value_helper->getMinTSDF();
}

float TSDF2D::GetWeight(const Eigen::Array2i& cell_index) const {
  if (limits_.Contains(cell_index)) {
    return value_helper->ValueToWeight(
        weight_cells_[ToFlatIndex(cell_index, limits_)]);
  }
  return value_helper->getMinWeight();
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
float TSDF2D::GetMaxTSDF() const { return value_helper->getMaxTSDF(); }
float TSDF2D::GetMinTSDF() const { return value_helper->getMinTSDF(); }
float TSDF2D::GetMaxWeight() const { return value_helper->getMaxWeight(); }
float TSDF2D::GetMinWeight() const { return value_helper->getMinWeight(); }

// Returns true if the probability at the specified index is known.
bool TSDF2D::IsKnown(const Eigen::Array2i& cell_index) const {
  return limits_.Contains(cell_index) &&
         weight_cells_[ToFlatIndex(cell_index, limits_)] !=
             value_helper->getUnknownWeightValue();
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
  details.set_truncation_distance(truncation_distance_);
  details.set_maximum_weight(max_weight_);
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
    std::vector<uint16> new_tsdf_cells(new_size,
                                       value_helper->getUnknownTSDFValue());
    std::vector<uint16> new_weight_cells(new_size,
                                         value_helper->getUnknownTSDFValue());
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
