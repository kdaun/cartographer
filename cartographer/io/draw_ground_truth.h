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

#ifndef CARTOGRAPHER_IO_DRAW_GROUND_TRUTH_H_
#define CARTOGRAPHER_IO_DRAW_GROUND_TRUTH_H_

#include "cairo/cairo.h"
#include "cartographer/ground_truth/proto/relations.pb.h"
#include "cartographer/io/color.h"
#include "cartographer/mapping/proto/trajectory.pb.h"
#include "cartographer/transform/rigid_transform.h"

namespace cartographer {
namespace io {

using PoseToPixelFunction =
    std::function<Eigen::Array2i(const transform::Rigid3d& pose)>;

// Draws the 'trajectory' with the given 'color' onto 'surface'. The
// 'pose_to_pixel' function must translate a trajectory node's position into the
// pixel on 'surface'.
void DrawGroundTruth(const mapping::proto::Trajectory& trajectory,
                     const ground_truth::proto::GroundTruth& ground_truth,
                     const FloatColor& color,
                     const PoseToPixelFunction& pose_to_pixel,
                     cairo_surface_t* surface);

}  // namespace io
}  // namespace cartographer

#endif  // CARTOGRAPHER_IO_DRAW_GROUND_TRUTH_H_
