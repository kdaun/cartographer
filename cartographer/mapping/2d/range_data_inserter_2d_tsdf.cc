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

#include <cstdlib>

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "cartographer/mapping/2d/xy_index.h"
#include "cartographer/mapping/internal/2d/ray_casting.h"
#include "cartographer/mapping/tsdf_values.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping {

RangeDataInserter2DTSDF::RangeDataInserter2DTSDF(
    const proto::RangeDataInserterOptions2D& options)
    : options_(options) {}

// Utility function for obtaining helper variables for the raytracing step.
inline void GetRaytracingHelperVariables(
    const Eigen::Vector2f& observation_origin,
    const Eigen::Vector2f& observation_ray, float t_start, float t_end,
    double grid_size_inv, Eigen::Vector2f* grid_index,
    Eigen::Vector2f* grid_step, Eigen::Vector2f* t_max,
    Eigen::Vector2f* t_delta) {
  CHECK(grid_index != nullptr);
  CHECK(grid_step != nullptr);
  CHECK(t_max != nullptr);
  CHECK(t_delta != nullptr);

  // Start and end of voxel traversal region.
  const Eigen::Vector2f traversal_start =
      observation_origin + t_start * observation_ray;
  const Eigen::Vector2f traversal_end =
      observation_origin + t_end * observation_ray;
  const Eigen::Vector2f traversal_start_scaled =
      traversal_start * grid_size_inv;
  const Eigen::Vector2f traversal_end_scaled = traversal_end * grid_size_inv;
  const Eigen::Vector2f traversal_ray_scaled =
      traversal_end_scaled - traversal_start_scaled;

  const Eigen::Vector2f traversal_ray_scaled_inv(
      1.f / traversal_ray_scaled.x(), 1.f / traversal_ray_scaled.y());

  // Calculate starting voxel index.
  *grid_index << std::floor(traversal_start_scaled.x()),
      std::floor(traversal_start_scaled.y());

  // Calculate the direction of the traversal step.
  *grid_step << (traversal_ray_scaled.x() < 0.0f ? -1 : 1),
      (traversal_ray_scaled.y() < 0.0f ? -1 : 1);

  const Eigen::Vector2f adjustment(grid_step->x() > 0 ? 1 : 0,
                                   grid_step->y() > 0 ? 1 : 0);

  const Eigen::Vector2f grid_index_adjusted = *grid_index + adjustment;

  *t_max = (grid_index_adjusted.cast<float>() - traversal_start_scaled)
               .cwiseProduct(traversal_ray_scaled_inv);

  // Determine how far we must travel along the ray before we have crossed a
  // grid cell.
  *t_delta = grid_step->cast<float>().cwiseProduct(traversal_ray_scaled_inv);
}

void RangeDataInserter2DTSDF::Insert(const sensor::RangeData& range_data,
                                     TSDF2D* const tsdf) const {
  const float truncation_distance = options_.tsdf().truncation_distance();
  Eigen::AlignedBox2f bounding_box(range_data.origin.head<2>());
  for (const Eigen::Vector3f& hit : range_data.returns) {
    const Eigen::Vector3f direction = (hit - range_data.origin).normalized();
    const Eigen::Vector3f end_position = hit + truncation_distance * direction;
    bounding_box.extend(end_position.head<2>());
  }
  constexpr float kPadding = 1e-6f;
  tsdf->GrowLimits(bounding_box.min() - kPadding * Eigen::Vector2f::Ones());
  tsdf->GrowLimits(bounding_box.max() + kPadding * Eigen::Vector2f::Ones());
  bool use_experimental_traversal = true;
  if (!use_experimental_traversal) {
    const Eigen::Vector2f origin = range_data.origin.head<2>();
    const Eigen::Array2i origin_cell = tsdf->limits().GetCellIndex(origin);
    for (const Eigen::Vector3f& hit_3d : range_data.returns) {
      const Eigen::Vector2f hit = hit_3d.head<2>();
      const Eigen::Vector2f direction = (hit - origin).normalized();
      float ray_length = (hit - origin).norm();
      const Eigen::Vector2f end_position =
          hit + truncation_distance * direction;
      const Eigen::Array2i end_cell = tsdf->limits().GetCellIndex(end_position);
      const Eigen::Array2i delta = end_cell - origin_cell;
      const int num_samples = delta.cwiseAbs().maxCoeff();
      CHECK_LT(num_samples, 1 << 15);
      // 'num_samples' is the number of samples we equi-distantly place on the
      // line between 'origin' and 'hit'. (including a fractional part for sub-
      // voxels) It is chosen so that between two samples we change from one
      // voxel to the next on the fastest changing dimension.
      for (int position = 0; position < num_samples; ++position) {
        const Eigen::Array2i cell =
            origin_cell + delta * position / num_samples;
        // float resolution = tsdf->limits().resolution();
        const Eigen::Vector2f cell_position =
            tsdf->limits().GetCellCenter(cell);
        // float distance = (hit - cell_position).dot(direction);
        float distance_sign =
            (hit - cell_position).dot(direction) > 0 ? 1.0 : -1.0;
        float distance = distance_sign * (hit - cell_position).norm();
        if (distance > truncation_distance) {
          distance = truncation_distance;
        } else if (distance < -truncation_distance) {
          distance = -truncation_distance;
        }
        UpdateCell(tsdf, cell, distance, ray_length);
      }
    }
  } else {
    const Eigen::Vector2f origin = range_data.origin.head<2>();
//    const Eigen::Array2i origin_cell = tsdf->limits().GetCellIndex(origin);
    for (const Eigen::Vector3f& hit_3d : range_data.returns) {
      const Eigen::Vector2f hit = hit_3d.head<2>();

      const double range = (hit - origin).norm();
      // Start and end of truncation parameter, normalized by range.
      const double range_inv = 1.0 / range;
      const double t_truncation_distance = truncation_distance * range_inv;
      const float t_start = 0.0;
      const float t_end = 1.0 + t_truncation_distance;

      // Calculate helper variables for raytracing loop.
      const Eigen::Vector2f ray = hit - origin;
      float resolution = tsdf->limits().resolution();
      float voxel_size_inv = 1. / resolution;
      Eigen::Vector2f voxel_index, voxel_step;
      Eigen::Vector2f t_max, t_delta;
      GetRaytracingHelperVariables(origin, ray, t_start, t_end, voxel_size_inv,
                                   &voxel_index, &voxel_step, &t_max, &t_delta);

      float t = 0.0;

      while (t < 1.0f + t_truncation_distance) {
        int min_coeff_idx;
        const float t_next = t_max.minCoeff(&min_coeff_idx);
        //const Eigen::Array2i origin_cell = tsdf->limits().GetCellIndex(origin);
        //const Eigen::Vector2f cell_position = tsdf->limits().GetCellCenter(cell);

        Eigen::Vector2f sampling_point = origin + t * ray;
        Eigen::Array2i cell_index = tsdf->limits().GetCellIndex(sampling_point);
        Eigen::Vector2f cell_center = tsdf->limits().GetCellCenter(cell_index);

        // Distance in meters from sensor origin to the effective traversal point.
        // Can be estimated in two different ways.
        float distance;
        constexpr bool kUseCentroidTraversalDistance = true;
        if (kUseCentroidTraversalDistance) {
          // Estimate using the distance to the voxel centroid.
          distance = (cell_center - origin).norm();
        } else {
          // Estimate using the distance midpoint of the voxel entry and voxel exit
          // points.
          const float t_mid = (t + t_next) / 2.0;
          const float t_adj = t_start + t_mid * (t_end - t_start);
          distance = t_adj * range;
        }

        UpdateCell(tsdf, cell_index, range - distance, range);


        // Advance to next voxel.
        t = t_next;
        //voxel_index(min_coeff_idx) += voxel_step(min_coeff_idx);
        //local_voxel_index(min_coeff_idx) += voxel_step(min_coeff_idx);
        t_max(min_coeff_idx) += t_delta(min_coeff_idx);


      }

    }
  }
  tsdf->FinishUpdate();
}

void RangeDataInserter2DTSDF::UpdateCell(TSDF2D* const tsdf,
                                         const Eigen::Array2i& cell,
                                         float update_sdf,
                                         float ray_length) const {
  float update_weight = 0.f;
  if (proto::CONSTANT_WEIGHT == options_.tsdf().range_data_inserter_type()) {
    update_weight = ComputeWeightConstant(update_sdf);
  } else if (proto::LINEAR_WEIGHT ==
      options_.tsdf().range_data_inserter_type()) {
    update_weight = ComputeWeightLinear(update_sdf, ray_length);
  } else if (proto::QUADRATIC_WEIGHT ==
      options_.tsdf().range_data_inserter_type()) {
    update_weight = ComputeWeightQuadratic(update_sdf, ray_length);
  }

  // float scaling = std::abs(tsdf->GetTSDF(cell) / (0.3f + update_sdf));
  // update_weight *= scaling;
  float updated_weight = tsdf->GetWeight(cell) + update_weight;
  float updated_sdf = updated_weight > 0.f
                          ? (tsdf->GetTSDF(cell) * tsdf->GetWeight(cell) +
                             update_sdf * update_weight) /
                                updated_weight
                          : tsdf->GetTSDF(cell);
  tsdf->UpdateCell(cell, updated_sdf, updated_weight);
}

float RangeDataInserter2DTSDF::ComputeWeightConstant(float sdf) const {
  float behind_surface_factor = 1.0;
  if (-options_.tsdf().behind_surface_distance() > sdf) {
    behind_surface_factor = (sdf + options_.tsdf().truncation_distance()) /
                            (options_.tsdf().truncation_distance() -
                             options_.tsdf().behind_surface_distance());
  }
  return behind_surface_factor * options_.tsdf().update_weight();
}

float RangeDataInserter2DTSDF::ComputeWeightLinear(float sdf,
                                                      float ray_length) const {
  float weight = 0.f;
  if (-options_.tsdf().behind_surface_distance() < sdf) {
    weight = options_.tsdf().update_weight() / ray_length;
  } else if (-options_.tsdf().truncation_distance() < sdf) {
    float behind_surface_factor =
        (sdf + options_.tsdf().truncation_distance()) /
            (options_.tsdf().truncation_distance() -
                options_.tsdf().behind_surface_distance());
    weight = behind_surface_factor * options_.tsdf().update_weight() /
        ray_length;
  }
  return weight;
}

float RangeDataInserter2DTSDF::ComputeWeightQuadratic(float sdf,
                                                      float ray_length) const {
  float weight = 0.f;
  if (-options_.tsdf().behind_surface_distance() < sdf) {
    weight = options_.tsdf().update_weight() / (std::pow(ray_length, 2));
  } else if (-options_.tsdf().truncation_distance() < sdf) {
    float behind_surface_factor =
        (sdf + options_.tsdf().truncation_distance()) /
        (options_.tsdf().truncation_distance() -
         options_.tsdf().behind_surface_distance());
    weight = behind_surface_factor * options_.tsdf().update_weight() /
             (std::pow(ray_length, 2));
  }
  return weight;
}

}  // namespace mapping
}  // namespace cartographer
