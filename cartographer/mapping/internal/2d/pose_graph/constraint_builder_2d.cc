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

#include "cartographer/mapping/internal/2d/pose_graph/constraint_builder_2d.h"

#include <cartographer/mapping/map_builder.h>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <string>

#include "Eigen/Eigenvalues"
#include "cartographer/common/make_unique.h"
#include "cartographer/common/math.h"
#include "cartographer/common/thread_pool.h"
#include "cartographer/mapping/2d/scan_matching/proto/ceres_scan_matcher_options_2d.pb.h"
#include "cartographer/mapping/2d/scan_matching/proto/fast_correlative_scan_matcher_options_2d.pb.h"
#include "cartographer/metrics/counter.h"
#include "cartographer/metrics/gauge.h"
#include "cartographer/metrics/histogram.h"
#include "cartographer/transform/transform.h"
#include "glog/logging.h"

namespace cartographer {
namespace mapping {
namespace pose_graph {
namespace {
static std::random_device rd;
static std::default_random_engine e1(42);
}  // namespace

static auto* kConstraintsSearchedMetric = metrics::Counter::Null();
static auto* kConstraintsFoundMetric = metrics::Counter::Null();
static auto* kGlobalConstraintsSearchedMetric = metrics::Counter::Null();
static auto* kGlobalConstraintsFoundMetric = metrics::Counter::Null();
static auto* kQueueLengthMetric = metrics::Gauge::Null();
static auto* kConstraintScoresMetric = metrics::Histogram::Null();
static auto* kGlobalConstraintScoresMetric = metrics::Histogram::Null();

transform::Rigid2d ComputeSubmapPose(const Submap2D& submap) {
  return transform::Project2D(submap.local_pose());
}

ConstraintBuilder2D::ConstraintBuilder2D(
    const pose_graph::proto::ConstraintBuilderOptions& options,
    common::ThreadPool* const thread_pool)
    : options_(options),
      thread_pool_(thread_pool),
      sampler_(options.sampling_ratio()),
      ceres_scan_matcher_(options.ceres_scan_matcher_options()) {}

ConstraintBuilder2D::~ConstraintBuilder2D() {
  common::MutexLocker locker(&mutex_);
  CHECK_EQ(constraints_.size(), 0) << "WhenDone() was not called";
  CHECK_EQ(pending_computations_.size(), 0);
  CHECK_EQ(submap_queued_work_items_.size(), 0);
  CHECK(when_done_ == nullptr);
}

void ConstraintBuilder2D::MaybeAddConstraint(
    const SubmapId& submap_id, const Submap2D* const submap,
    const NodeId& node_id, const TrajectoryNode::Data* const constant_data,
    const transform::Rigid2d& initial_relative_pose) {
  if (initial_relative_pose.translation().norm() >
      options_.max_constraint_distance()) {
    return;
  }
  if (sampler_.Pulse()) {
    common::MutexLocker locker(&mutex_);
    constraints_.emplace_back();
    kQueueLengthMetric->Set(constraints_.size());
    auto* const constraint = &constraints_.back();
    ++pending_computations_[current_computation_];
    const int current_computation = current_computation_;
    ScheduleSubmapScanMatcherConstructionAndQueueWorkItem(
        submap_id, &submap->probability_grid(), [=]() EXCLUDES(mutex_) {
          ComputeConstraint(submap_id, submap, node_id,
                            false, /* match_full_submap */
                            constant_data, initial_relative_pose, constraint);
          FinishComputation(current_computation);
        });
  }
}

void ConstraintBuilder2D::MaybeAddGlobalConstraint(
    const SubmapId& submap_id, const Submap2D* const submap,
    const NodeId& node_id, const TrajectoryNode::Data* const constant_data) {
  common::MutexLocker locker(&mutex_);
  constraints_.emplace_back();
  kQueueLengthMetric->Set(constraints_.size());
  auto* const constraint = &constraints_.back();
  ++pending_computations_[current_computation_];
  const int current_computation = current_computation_;
  ScheduleSubmapScanMatcherConstructionAndQueueWorkItem(
      submap_id, &submap->probability_grid(), [=]() EXCLUDES(mutex_) {
        ComputeConstraint(
            submap_id, submap, node_id, true, /* match_full_submap */
            constant_data, transform::Rigid2d::Identity(), constraint);
        FinishComputation(current_computation);
      });
}

void ConstraintBuilder2D::NotifyEndOfNode() {
  common::MutexLocker locker(&mutex_);
  ++current_computation_;
}

void ConstraintBuilder2D::WhenDone(
    const std::function<void(const ConstraintBuilder2D::Result&)>& callback) {
  common::MutexLocker locker(&mutex_);
  CHECK(when_done_ == nullptr);
  when_done_ =
      common::make_unique<std::function<void(const Result&)>>(callback);
  ++pending_computations_[current_computation_];
  const int current_computation = current_computation_;
  thread_pool_->Schedule(
      [this, current_computation] { FinishComputation(current_computation); });
}

void ConstraintBuilder2D::ScheduleSubmapScanMatcherConstructionAndQueueWorkItem(
    const SubmapId& submap_id, const Grid2D* const submap,
    const std::function<void()>& work_item) {
  if (submap_scan_matchers_[submap_id].fast_correlative_scan_matcher !=
      nullptr) {
    thread_pool_->Schedule(work_item);
  } else {
    submap_queued_work_items_[submap_id].push_back(work_item);
    if (submap_queued_work_items_[submap_id].size() == 1) {
      thread_pool_->Schedule(
          [=]() { ConstructSubmapScanMatcher(submap_id, submap); });
    }
  }
}

void ConstraintBuilder2D::ConstructSubmapScanMatcher(
    const SubmapId& submap_id, const Grid2D* const submap) {
  auto submap_scan_matcher =
      common::make_unique<scan_matching::FastCorrelativeScanMatcher2D>(
          *submap, options_.fast_correlative_scan_matcher_options());
  common::MutexLocker locker(&mutex_);
  submap_scan_matchers_[submap_id] = {submap, std::move(submap_scan_matcher)};
  for (const std::function<void()>& work_item :
       submap_queued_work_items_[submap_id]) {
    thread_pool_->Schedule(work_item);
  }
  submap_queued_work_items_.erase(submap_id);
}

const ConstraintBuilder2D::SubmapScanMatcher*
ConstraintBuilder2D::GetSubmapScanMatcher(const SubmapId& submap_id) {
  common::MutexLocker locker(&mutex_);
  const SubmapScanMatcher* submap_scan_matcher =
      &submap_scan_matchers_[submap_id];
  CHECK(submap_scan_matcher->fast_correlative_scan_matcher != nullptr);
  return submap_scan_matcher;
}

void ConstraintBuilder2D::ComputeConstraint(
    const SubmapId& submap_id, const Submap2D* const submap,
    const NodeId& node_id, bool match_full_submap,
    const TrajectoryNode::Data* const constant_data,
    const transform::Rigid2d& initial_relative_pose,
    std::unique_ptr<ConstraintBuilder2D::Constraint>* constraint) {
  const transform::Rigid2d initial_pose =
      ComputeSubmapPose(*submap) * initial_relative_pose;
  const SubmapScanMatcher* const submap_scan_matcher =
      GetSubmapScanMatcher(submap_id);

  // The 'constraint_transform' (submap i <- node j) is computed from:
  // - a 'filtered_gravity_aligned_point_cloud' in node j,
  // - the initial guess 'initial_pose' for (map <- node j),
  // - the result 'pose_estimate' of Match() (map <- node j).
  // - the ComputeSubmapPose() (map <- submap i)
  float score = 0.;
  transform::Rigid2d pose_estimate = transform::Rigid2d::Identity();

  // Compute 'pose_estimate' in three stages:
  // 1. Fast estimate using the fast correlative scan matcher.
  // 2. Prune if the score is too low.
  // 3. Refine.
  if (match_full_submap) {
    kGlobalConstraintsSearchedMetric->Increment();
    if (submap_scan_matcher->fast_correlative_scan_matcher->MatchFullSubmap(
            constant_data->filtered_gravity_aligned_point_cloud,
            options_.global_localization_min_score(), &score, &pose_estimate)) {
      CHECK_GT(score, options_.global_localization_min_score());
      CHECK_GE(node_id.trajectory_id, 0);
      CHECK_GE(submap_id.trajectory_id, 0);
      kGlobalConstraintsFoundMetric->Increment();
      kGlobalConstraintScoresMetric->Observe(score);
    } else {
      return;
    }
  } else {
    kConstraintsSearchedMetric->Increment();
    if (submap_scan_matcher->fast_correlative_scan_matcher->Match(
            initial_pose, constant_data->filtered_gravity_aligned_point_cloud,
            options_.min_score(), &score, &pose_estimate)) {
      // We've reported a successful local match.
      CHECK_GT(score, options_.min_score());
      kConstraintsFoundMetric->Increment();
      kConstraintScoresMetric->Observe(score);
    } else {
      return;
    }
  }
  {
    common::MutexLocker locker(&mutex_);
    score_histogram_.Add(score);
  }

  // Use the CSM estimate as both the initial and previous pose. This has the
  // effect that, in the absence of better information, we prefer the original
  // CSM estimate.
  ceres::Solver::Summary unused_summary;
  ceres_scan_matcher_.Match(pose_estimate.translation(), pose_estimate,
                            constant_data->filtered_gravity_aligned_point_cloud,
                            *submap_scan_matcher->grid, &pose_estimate,
                            &unused_summary);
  for (const auto& evaluation_constraint :
       cartographer::mapping::evaluation_constraints) {
    // todo(kdaun) add evaluation related stuff
    if (evaluation_constraint.submap_id == submap_id &&
        evaluation_constraint.node_id == node_id) {
      LOG(INFO) << "evaluating constraint...";
      // map cost
      float min_trans = -0.02;
      float max_trans = 0.02;
      float resolution_trans = 0.0001;
      float min_rot = 0.f;
      float max_rot = 2.f * M_PI;
      float resolution_rot = 0.1;
      std::ofstream log_file;
      std::string log_file_path;
      time_t seconds;
      time(&seconds);
      log_file_path = "gradient_pg_" + std::to_string(submap_id.submap_index) +
                      std::to_string(node_id.node_index) +
                      std::to_string(seconds) + ".csv";
      log_file.open(log_file_path);

      LOG(INFO) << "writing " << log_file_path;
      for (float x = min_trans; x < max_trans; x += resolution_trans) {
        for (float y = min_trans; y < max_trans; y += resolution_trans) {
          const Eigen::Vector2d target_translation = {0., 0.};
          ceres::Solver::Summary summary;
          double cost;
          std::vector<double> jacobians;
          std::vector<Eigen::Vector3f> sample_cloud = {{x, y, 0.f}};

          transform::Rigid2d pose_estimate_sample =
          transform::Rigid2d({(double)x, double(y)}, 0.0) * pose_estimate ;
          ceres_scan_matcher_.Evaluate(
              pose_estimate_sample.translation(), pose_estimate_sample,
              constant_data->filtered_gravity_aligned_point_cloud,
              submap->probability_grid(), &cost, NULL, &jacobians);
          double theta = 0.;
          std::vector<double> result = {x, y, theta, cost};
          result.insert(result.end(), jacobians.begin(), jacobians.end());

          for (auto& element : result) {
            log_file << element << ",";
          }
          log_file << "\n";
        }
      }
      log_file.close();

      // WHOLE CLOUD TSDF
      log_file_path =
          "gradient_tsdf_" + std::to_string(submap_id.submap_index) +
          std::to_string(node_id.node_index) + std::to_string(seconds) + ".csv";
      log_file.open(log_file_path);
      LOG(INFO) << "writing " << log_file_path;

      for (float x = min_trans; x < max_trans; x += resolution_trans) {
        for (float y = min_trans; y < max_trans; y += resolution_trans) {
          const Eigen::Vector2d target_translation = {0., 0.};
          ceres::Solver::Summary summary;
          double cost;
          std::vector<double> jacobians;

          transform::Rigid2d pose_estimate_sample =
              transform::Rigid2d({(double)x, double(y)}, 0.0) * pose_estimate;
          ceres_scan_matcher_.Evaluate(
              pose_estimate_sample.translation(), pose_estimate_sample,
              constant_data->filtered_gravity_aligned_point_cloud,
              submap->tsdf(), &cost, NULL, &jacobians);
          double theta = 0.;
          std::vector<double> result = {x, y, theta, cost};
          result.insert(result.end(), jacobians.begin(), jacobians.end());

          for (auto& element : result) {
            log_file << element << ",";
          }
          log_file << "\n";
        }
      }
      log_file.close();

      min_trans = -20.0;
      max_trans = 20.0;
      resolution_trans = 0.05;

      // SINGLE POINT TSDF
      log_file_path = "single_point_gradient_tsdf_" +
                      std::to_string(submap_id.submap_index) +
                      std::to_string(node_id.node_index) +
                      std::to_string(seconds) + ".csv";
      log_file.open(log_file_path);
      LOG(INFO) << "writing " << log_file_path;
      for (float x = min_trans; x < max_trans; x += resolution_trans) {
        for (float y = min_trans; y < max_trans; y += resolution_trans) {
          const Eigen::Vector2d target_translation = {0., 0.};
          ceres::Solver::Summary summary;
          double cost;
          std::vector<double> jacobians;
          std::vector<Eigen::Vector3f> sample_cloud = {{x, y, 0.f}};

          transform::Rigid2d pose_estimate_sample =
              pose_estimate * transform::Rigid2d({(double)x, double(y)}, 0.0);
          ceres_scan_matcher_.Evaluate(pose_estimate_sample.translation(),
                                       pose_estimate_sample, {{0.0, 0.0, 0.0}},
                                       submap->tsdf(), &cost, NULL, &jacobians);
          double theta = 0.;
          std::vector<double> result = {x, y, theta, cost};
          result.insert(result.end(), jacobians.begin(), jacobians.end());

          for (auto& element : result) {
            log_file << element << ",";
          }
          log_file << "\n";
        }
      }
      log_file.close();

      // SINGLE POINT PG
      log_file_path =
          "single_point_gradient_pg_" + std::to_string(submap_id.submap_index) +
          std::to_string(node_id.node_index) + std::to_string(seconds) + ".csv";
      log_file.open(log_file_path);
      LOG(INFO) << "writing " << log_file_path;
      for (float x = min_trans; x < max_trans; x += resolution_trans) {
        for (float y = min_trans; y < max_trans; y += resolution_trans) {
          const Eigen::Vector2d target_translation = {0., 0.};
          ceres::Solver::Summary summary;
          double cost;
          std::vector<double> jacobians;
          std::vector<Eigen::Vector3f> sample_cloud = {{x, y, 0.f}};

          transform::Rigid2d pose_estimate_sample =
              pose_estimate * transform::Rigid2d({(double)x, double(y)}, 0.0);
          ceres_scan_matcher_.Evaluate(pose_estimate_sample.translation(),
                                       pose_estimate_sample, {{0.0, 0.0, 0.0}},
                                       submap->probability_grid(), &cost, NULL,
                                       &jacobians);
          double theta = 0.;
          std::vector<double> result = {x, y, theta, cost};
          result.insert(result.end(), jacobians.begin(), jacobians.end());

          for (auto& element : result) {
            log_file << element << ",";
          }
          log_file << "\n";
        }
      }
      log_file.close();

      // DISPLACEMENT
      log_file_path = "displacement_" + std::to_string(submap_id.submap_index) +
                      std::to_string(node_id.node_index) +
                      std::to_string(seconds) + ".csv";
      log_file.open(log_file_path);
      LOG(INFO) << "writing " << log_file_path;

      log_file << "grid_type,grid_resolution,initial_error_trans,initial_error_"
                  "angle,matched_error_trans,matched_error_x,matched_error_y,"
                  "matched_error_angle,solver_iterations,"
                  "matching_time\n";

      // std::vector<double> trans_errors = {0.05, 0.1, 0.25, 0.5};
      std::vector<double> trans_errors = {0.05, 0.1};
      // std::vector<double> rot_errors = {0.2 * M_PI_4, 0.4 * M_PI_4,
      //                                  0.6 * M_PI_4, 0.8 * M_PI_4};
      std::vector<double> rot_errors = {0.0};
      int n_samples = 100;

      ceres_scan_matcher_.options_.set_translation_weight(0.0);
      ceres_scan_matcher_.options_.set_rotation_weight(0.0);

      for (double error_trans : trans_errors) {
        for (double error_rot : rot_errors) {
          for (int i_sample = 0; i_sample < n_samples; ++i_sample) {
            std::uniform_real_distribution<double> error_translation_direction(
                -M_PI, M_PI);
            double e_scale = error_trans == 0.0 ? 0.0 : error_trans;
            double e_orientation = error_translation_direction(e1);
            double e_x = std::cos(e_orientation) * e_scale;
            double e_y = std::sin(e_orientation) * e_scale;
            double e_rotation_direction = error_translation_direction(e1);
            double e_theta = e_rotation_direction > 0 ? error_rot : -error_rot;
            const cartographer::transform::Rigid2d initial_displacement =
                cartographer::transform::Rigid2d({e_x, e_y}, e_theta);
            ceres::Solver::Summary summary;

            transform::Rigid2d displaced_pose_estimate =
                initial_displacement * pose_estimate;
            transform::Rigid2d matched_displaced_pose_estimate;
            ceres_scan_matcher_.Match(
                displaced_pose_estimate.translation(), displaced_pose_estimate,
                constant_data->filtered_gravity_aligned_point_cloud,
                submap->probability_grid(), &matched_displaced_pose_estimate,
                &summary);
            double initial_trans_error =
                initial_displacement.translation().norm();
            double matching_trans_error =
                (pose_estimate.inverse() * matched_displaced_pose_estimate)
                    .translation()
                    .norm();
            double initial_error_angle =
                initial_displacement.rotation().smallestAngle();
            double matched_error_angle =
                (pose_estimate.inverse() * matched_displaced_pose_estimate)
                    .rotation()
                    .smallestAngle();
            double solver_iterations =
                summary.num_successful_steps + summary.num_unsuccessful_steps;
            double matching_time = summary.minimizer_time_in_seconds;
            log_file
                << "PROBABILITY_GRID,0.05," << initial_trans_error << ","
                << initial_error_angle << "," << matching_trans_error << ","
                << (pose_estimate.inverse() * matched_displaced_pose_estimate)
                       .translation()[0]
                << ","
                << (pose_estimate.inverse() * matched_displaced_pose_estimate)
                       .translation()[1]
                << "," << matched_error_angle << "," << solver_iterations << ","
                << matching_time << "\n";

            displaced_pose_estimate = initial_displacement * pose_estimate;

            double default_occupied_space_weight =
                ceres_scan_matcher_.options_.occupied_space_weight();
            ceres_scan_matcher_.ceres_solver_options_.max_num_iterations = 200;
            ceres_scan_matcher_.ceres_solver_options_.function_tolerance = 1e-9;
            // LOG(INFO)<<"default_occupied_space_weight
            // "<<default_occupied_space_weight;
            ceres_scan_matcher_.options_.set_occupied_space_weight(60.);
            ceres_scan_matcher_.Match(
                displaced_pose_estimate.translation(), displaced_pose_estimate,
                constant_data->filtered_gravity_aligned_point_cloud,
                submap->tsdf(), &matched_displaced_pose_estimate, &summary);
            ceres_scan_matcher_.ceres_solver_options_.function_tolerance = 1e-6;
            ceres_scan_matcher_.ceres_solver_options_.max_num_iterations = 10;
            ceres_scan_matcher_.options_.set_occupied_space_weight(
                default_occupied_space_weight);
            initial_trans_error = initial_displacement.translation().norm();
            matching_trans_error =
                (pose_estimate.inverse() * matched_displaced_pose_estimate)
                    .translation()
                    .norm();
            initial_error_angle =
                initial_displacement.rotation().smallestAngle();
            matched_error_angle =
                (pose_estimate.inverse() * matched_displaced_pose_estimate)
                    .rotation()
                    .smallestAngle();
            solver_iterations =
                summary.num_successful_steps + summary.num_unsuccessful_steps;
            matching_time = summary.minimizer_time_in_seconds;
            if (summary.termination_type != 1) {
              LOG(INFO) << "TERMINATION TYPE: " << summary.termination_type
                        << " " << summary.message;
            }
            // LOG(INFO)<<summary.FullReport();
            log_file << "TSDF,0.05," << initial_trans_error << ","
                     << initial_error_angle << "," << matching_trans_error << ","
                     << (pose_estimate.inverse() *
                         matched_displaced_pose_estimate)
                            .translation()[0]
                     << ","
                     << (pose_estimate.inverse() *
                         matched_displaced_pose_estimate)
                            .translation()[1]
                     << "," << matched_error_angle << "," << solver_iterations
                     << "," << matching_time << "\n";
          }
        }
      }
      log_file.close();
    }
  }

  const transform::Rigid2d constraint_transform =
      ComputeSubmapPose(*submap).inverse() * pose_estimate;
  constraint->reset(new Constraint{submap_id,
                                   node_id,
                                   {transform::Embed3D(constraint_transform),
                                    options_.loop_closure_translation_weight(),
                                    options_.loop_closure_rotation_weight()},
                                   Constraint::INTER_SUBMAP});

  if (options_.log_matches()) {
    std::ostringstream info;
    info << "Node " << node_id << " with "
         << constant_data->filtered_gravity_aligned_point_cloud.size()
         << " points on submap " << submap_id << std::fixed;
    if (match_full_submap) {
      info << " matches";
    } else {
      const transform::Rigid2d difference =
          initial_pose.inverse() * pose_estimate;
      info << " differs by translation " << std::setprecision(2)
           << difference.translation().norm() << " rotation "
           << std::setprecision(3) << std::abs(difference.normalized_angle());
    }
    info << " with score " << std::setprecision(1) << 100. * score << "%.";
    LOG(INFO) << info.str();
  }
}

void ConstraintBuilder2D::FinishComputation(const int computation_index) {
  Result result;
  std::unique_ptr<std::function<void(const Result&)>> callback;
  {
    common::MutexLocker locker(&mutex_);
    if (--pending_computations_[computation_index] == 0) {
      pending_computations_.erase(computation_index);
    }
    if (pending_computations_.empty()) {
      CHECK_EQ(submap_queued_work_items_.size(), 0);
      if (when_done_ != nullptr) {
        for (const std::unique_ptr<Constraint>& constraint : constraints_) {
          if (constraint != nullptr) {
            result.push_back(*constraint);
          }
        }
        if (options_.log_matches()) {
          LOG(INFO) << constraints_.size() << " computations resulted in "
                    << result.size() << " additional constraints.";
          LOG(INFO) << "Score histogram:\n" << score_histogram_.ToString(10);
        }
        constraints_.clear();
        callback = std::move(when_done_);
        when_done_.reset();
      }
    }
    kQueueLengthMetric->Set(constraints_.size());
  }
  if (callback != nullptr) {
    (*callback)(result);
  }
}

int ConstraintBuilder2D::GetNumFinishedNodes() {
  common::MutexLocker locker(&mutex_);
  if (pending_computations_.empty()) {
    return current_computation_;
  }
  return pending_computations_.begin()->first;
}

void ConstraintBuilder2D::DeleteScanMatcher(const SubmapId& submap_id) {
  common::MutexLocker locker(&mutex_);
  CHECK(pending_computations_.empty());
  submap_scan_matchers_.erase(submap_id);
}

void ConstraintBuilder2D::RegisterMetrics(metrics::FamilyFactory* factory) {
  auto* counts = factory->NewCounterFamily(
      "/mapping/2d/pose_graph/constraint_builder/constraints",
      "Constraints computed");
  kConstraintsSearchedMetric =
      counts->Add({{"search_region", "local"}, {"matcher", "searched"}});
  kConstraintsFoundMetric =
      counts->Add({{"search_region", "local"}, {"matcher", "found"}});
  kGlobalConstraintsSearchedMetric =
      counts->Add({{"search_region", "global"}, {"matcher", "searched"}});
  kGlobalConstraintsFoundMetric =
      counts->Add({{"search_region", "global"}, {"matcher", "found"}});
  auto* queue_length = factory->NewGaugeFamily(
      "/mapping/2d/pose_graph/constraint_builder/queue_length", "Queue length");
  kQueueLengthMetric = queue_length->Add({});
  auto boundaries = metrics::Histogram::FixedWidth(0.05, 20);
  auto* scores = factory->NewHistogramFamily(
      "/mapping/2d/pose_graph/constraint_builder/scores",
      "Constraint scores built", boundaries);
  kConstraintScoresMetric = scores->Add({{"search_region", "local"}});
  kGlobalConstraintScoresMetric = scores->Add({{"search_region", "global"}});
}

}  // namespace pose_graph
}  // namespace mapping
}  // namespace cartographer
