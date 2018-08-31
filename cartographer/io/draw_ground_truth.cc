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

#include "cartographer/io/draw_ground_truth.h"

#include "cartographer/io/image.h"
#include "cartographer/transform/transform.h"
#include "cartographer/transform/transform_interpolation_buffer.h"

namespace cartographer {
namespace io {

namespace {
transform::Rigid3d LookupTransform(
    const transform::TransformInterpolationBuffer&
        transform_interpolation_buffer,
    const common::Time time) {
  const common::Time earliest_time =
      transform_interpolation_buffer.earliest_time();
  if (transform_interpolation_buffer.Has(time)) {
    return transform_interpolation_buffer.Lookup(time);
  } else if (time < earliest_time) {
    return transform_interpolation_buffer.Lookup(earliest_time);
  }
  return transform_interpolation_buffer.Lookup(
      transform_interpolation_buffer.latest_time());
}
}  // namespace

void DrawGroundTruth(const mapping::proto::Trajectory& trajectory,
                     const ground_truth::proto::GroundTruth& ground_truth,
                     const FloatColor& color,
                     const PoseToPixelFunction& pose_to_pixel,
                     cairo_surface_t* surface) {
  if (trajectory.node_size() == 0) {
    return;
  }
  constexpr double kTrajectoryWidth = 4.;
  constexpr double kAlpha = 0.7;

  auto cr = ::cartographer::io::MakeUniqueCairoPtr(cairo_create(surface));

  cairo_set_source_rgba(cr.get(), color[0], color[1], color[2], kAlpha);
  cairo_set_line_width(cr.get(), kTrajectoryWidth);

  const transform::TransformInterpolationBuffer transform_interpolation_buffer(
      trajectory);

  for (const auto& relation : ground_truth.relation()) {
    const auto pose1 =
        LookupTransform(transform_interpolation_buffer,
                        common::FromUniversal(relation.timestamp1()));
    const auto pose2 =
        LookupTransform(transform_interpolation_buffer,
                        common::FromUniversal(relation.timestamp2()));
    transform::Rigid3d expected = transform::ToRigid3(relation.expected());
    /*transform::Rigid3d expected_corrected =
        transform::Rigid3d(-expected.translation(), expected.rotation());*/
    transform::Rigid3d expected_corrected = expected;
    const Eigen::Array2i pixel_pose1 = pose_to_pixel(pose1);
    const Eigen::Array2i pixel_expected =
        pose_to_pixel(pose1 * expected_corrected);
    const Eigen::Array2i pixel_pose2 = pose_to_pixel(pose2);

    cairo_set_source_rgba(cr.get(), color[0], color[1], color[2], kAlpha);
    cairo_move_to(cr.get(), pixel_pose1.x(), pixel_pose1.y());
    cairo_line_to(cr.get(), pixel_expected.x(), pixel_expected.y());
    cairo_stroke(cr.get());
  }

  for (const auto& relation : ground_truth.relation()) {
    const auto pose1 =
        LookupTransform(transform_interpolation_buffer,
                        common::FromUniversal(relation.timestamp1()));
    const auto pose2 =
        LookupTransform(transform_interpolation_buffer,
                        common::FromUniversal(relation.timestamp2()));
    transform::Rigid3d expected = transform::ToRigid3(relation.expected());
    /*transform::Rigid3d expected_corrected =
        transform::Rigid3d(-expected.translation(), expected.rotation());*/
    transform::Rigid3d expected_corrected = expected;
    const Eigen::Array2i pixel_pose1 = pose_to_pixel(pose1);
    const Eigen::Array2i pixel_expected =
        pose_to_pixel(pose1 * expected_corrected);
    const Eigen::Array2i pixel_pose2 = pose_to_pixel(pose2);

    cairo_set_source_rgba(cr.get(), 1, 0, 0, 1);
    cairo_move_to(cr.get(), pixel_pose2.x(), pixel_pose2.y());
    cairo_line_to(cr.get(), pixel_expected.x(), pixel_expected.y());
    cairo_stroke(cr.get());

    const transform::Rigid3d error =
        (pose1.inverse() * pose2) * expected_corrected.inverse();
    if (error.translation().squaredNorm() > 1.0)
      LOG(WARNING) << "HUGE ERROR " << error.translation().squaredNorm();
  }

  cairo_surface_flush(surface);
}

}  // namespace io
}  // namespace cartographer
