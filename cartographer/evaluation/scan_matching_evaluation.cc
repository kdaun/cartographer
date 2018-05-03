#include <cartographer/mapping/2d/scan_matching/proto/ceres_scan_matcher_options_2d.pb.h>

#include <fstream>
#include <random>
#include <string>

#include "cairo.h"

#include "cartographer/common/lua_parameter_dictionary.h"
#include "cartographer/common/lua_parameter_dictionary_test_helpers.h"
#include "cartographer/common/make_unique.h"
#include "cartographer/common/time.h"
#include "cartographer/evaluation/scan_cloud_generator.h"
#include "cartographer/mapping/2d/probability_grid.h"
#include "cartographer/mapping/2d/range_data_inserter_2d.h"
#include "cartographer/mapping/2d/range_data_inserter_2d_probability_grid.h"
#include "cartographer/mapping/2d/range_data_inserter_2d_tsdf.h"
#include "cartographer/mapping/internal/2d/scan_matching/ceres_scan_matcher_2d.h"
#include "cartographer/mapping/internal/2d/scan_matching/occupied_space_cost_function_2d.h"

namespace cartographer {
namespace evaluation {
namespace {
static int rendered_grid_id = 0;
static std::random_device rd;
static std::default_random_engine e1(42);
}  // namespace

struct Sample {
  cartographer::sensor::RangeData range_data;
  cartographer::transform::Rigid2d ground_truth_pose;
};

struct SampleResult {
  double initial_trans_error;
  double initial_rot_error;
  double matching_trans_error;
  double matching_rot_error;
  double matching_time;
  int matching_iterations;
};

void GenerateRangeData(const ScanCloudGenerator::ModelType model_type,
                       const Eigen::Vector2d size, const double resolution,
                       cartographer::sensor::RangeData& range_data) {
  cartographer::sensor::PointCloud scan_cloud;
  ScanCloudGenerator test_set_generator(resolution);
  switch (model_type) {
    case ScanCloudGenerator::ModelType::CIRCLE:
      test_set_generator.generateCircle(scan_cloud, size[0]);
    case ScanCloudGenerator::ModelType::SQUARE:
      test_set_generator.generateSquare(scan_cloud, size[0]);
    case ScanCloudGenerator::ModelType::RECTANGLE:
      test_set_generator.generateRectangle(scan_cloud, size[0], size[1]);
      range_data.returns = scan_cloud;
      range_data.origin = Eigen::Vector3f{0, 0, 0};
  }
}

Sample GenerateSample(double error_trans, double error_rot,
                      const ScanCloudGenerator::ModelType model_type,
                      const Eigen::Vector2d size, const double resolution) {
  cartographer::sensor::RangeData range_data;
  GenerateRangeData(model_type, size, resolution, range_data);

  std::uniform_real_distribution<double> error_translation_direction(-M_PI,
                                                                     M_PI);
  double orientation = error_translation_direction(e1);
  double scale = error_trans == 0.0 ? 0.0 : error_trans;
  double x = std::cos(orientation) * scale;
  double y = std::sin(orientation) * scale;
  Sample sample;
  sample.ground_truth_pose =
      cartographer::transform::Rigid2d({x, y}, error_rot);
  sensor::RangeData initial_pose_estimate_range_data =
      cartographer::sensor::TransformRangeData(
          range_data,
          transform::Embed3D(sample.ground_truth_pose.cast<float>()));
  sample.range_data = initial_pose_estimate_range_data;
  return sample;
}

void GenerateSampleSet(int n_training, int n_test, double error_trans,
                       double error_rot,
                       const ScanCloudGenerator::ModelType model_type,
                       const Eigen::Vector2d size, const double resolution,
                       std::vector<Sample>& training_set,
                       std::vector<Sample>& test_set) {
  for (int i = 0; i < n_training; ++i) {
    training_set.push_back(
        GenerateSample(0.0, 0.0, model_type, size, resolution));
  }
  for (int i = 0; i < n_test; ++i) {
    test_set.push_back(
        GenerateSample(error_trans, error_rot, model_type, size, resolution));
  }
}
template <typename GridType>
std::unique_ptr<GridType> generateGrid() {
  std::unique_ptr<GridType> grid =
      common::make_unique<GridType>(cartographer::mapping::MapLimits(
          0.05, Eigen::Vector2d(1., 1.),
          cartographer::mapping::CellLimits(40, 40)));
  return std::move(grid);
}

template <>
std::unique_ptr<cartographer::mapping::TSDF2D>
generateGrid<cartographer::mapping::TSDF2D>() {
  std::unique_ptr<cartographer::mapping::TSDF2D> grid =
      common::make_unique<cartographer::mapping::TSDF2D>(
          cartographer::mapping::MapLimits(
              0.05, Eigen::Vector2d(1., 1.),
              cartographer::mapping::CellLimits(40, 40)),
          0.3, 1.0);
  return std::move(grid);
}

void renderGridwithScan(
    const cartographer::mapping::ProbabilityGrid& grid, const Sample& sample,
    const cartographer::transform::Rigid2d& initial_transform,
    const cartographer::transform::Rigid2d& matched_transform) {
  sensor::RangeData initial_pose_estimate_range_data =
      cartographer::sensor::TransformRangeData(
          sample.range_data,
          transform::Embed3D(initial_transform.cast<float>()));
  sensor::RangeData matched_range_data =
      cartographer::sensor::TransformRangeData(
          sample.range_data,
          transform::Embed3D(matched_transform.cast<float>()));

  const cartographer::mapping::MapLimits& limits = grid.limits();
  double scale = 1. / limits.resolution();
  cairo_surface_t* grid_surface;
  cairo_t* grid_surface_context;

  int scaled_num_x_cells = limits.cell_limits().num_x_cells * scale;
  int scaled_num_y_cells = limits.cell_limits().num_y_cells * scale;
  grid_surface = cairo_image_surface_create(
      CAIRO_FORMAT_ARGB32, scaled_num_x_cells, scaled_num_y_cells);
  grid_surface_context = cairo_create(grid_surface);
  cairo_device_to_user_distance(grid_surface_context, &scale, &scale);
  for (int ix = 0; ix < scaled_num_x_cells; ++ix) {
    for (int iy = 0; iy < scaled_num_y_cells; ++iy) {
      float p = 1. - grid.GetProbability({iy, ix});
      cairo_set_source_rgb(grid_surface_context, p, p, p);
      cairo_rectangle(grid_surface_context, scale * (float(ix)),
                      scale * ((float)iy), scale, scale);
      cairo_fill(grid_surface_context);
    }
  }

  cairo_set_source_rgb(grid_surface_context, 0.8, 0.0, 0);
  for (auto& scan : initial_pose_estimate_range_data.returns) {
    float x = scale * (limits.max().x() - scan[0]);
    float y = scale * (limits.max().y() - scan[1]);
    cairo_rectangle(grid_surface_context, (x - 0.15) * scale,
                    (y - 0.15) * scale, 0.3 * scale, 0.3 * scale);
  }
  cairo_fill(grid_surface_context);

  cairo_set_source_rgb(grid_surface_context, 0.0, 0.8, 0);
  for (auto& scan : matched_range_data.returns) {
    float x = scale * (limits.max().x() - scan[0]);
    float y = scale * (limits.max().y() - scan[1]);
    cairo_rectangle(grid_surface_context, (x - 0.15) * scale,
                    (y - 0.15) * scale, 0.3 * scale, 0.3 * scale);
  }

  cairo_fill(grid_surface_context);

  time_t seconds;
  time(&seconds);
  std::string filename = "grid_with_inserted_cloud" + std::to_string(seconds) +
                         std::to_string(rendered_grid_id) + ".png";
  rendered_grid_id++;
  cairo_surface_write_to_png(grid_surface, filename.c_str());
}

void renderGridwithScan(
    const cartographer::mapping::TSDF2D& grid, const Sample& sample,
    const cartographer::transform::Rigid2d& initial_transform,
    const cartographer::transform::Rigid2d& matched_transform) {
  sensor::RangeData initial_pose_estimate_range_data =
      cartographer::sensor::TransformRangeData(
          sample.range_data,
          transform::Embed3D(initial_transform.cast<float>()));
  sensor::RangeData
      matched_range_data = cartographer::sensor::TransformRangeData( sample.range_data,
                                                                     transform::Embed3D(matched_transform.cast<float>()));

  const cartographer::mapping::MapLimits& limits = grid.limits();
  double scale = 1. / limits.resolution();
  cairo_surface_t* grid_surface;
  cairo_t* grid_surface_context;

  int scaled_num_x_cells = limits.cell_limits().num_x_cells * scale;
  int scaled_num_y_cells = limits.cell_limits().num_y_cells * scale;
  grid_surface = cairo_image_surface_create(
      CAIRO_FORMAT_ARGB32, scaled_num_x_cells, scaled_num_y_cells);
  grid_surface_context = cairo_create(grid_surface);
  cairo_device_to_user_distance(grid_surface_context, &scale, &scale);
  for (int ix = 0; ix < scaled_num_x_cells; ++ix) {
    for (int iy = 0; iy < scaled_num_y_cells; ++iy) {
      float r = 1.f;
      float g = 1.f;
      float b = 1.f;
      float normalized_tsdf = grid.GetTSDF({iy, ix}) / grid.GetMaxTSDF();
      if (normalized_tsdf > 0.f) {
        g = 1. - std::pow(std::abs(normalized_tsdf), 0.5);
        b = g;
      } else {
        r = 1. - std::pow(std::abs(normalized_tsdf), 0.5);
        g = r;
      }
      cairo_set_source_rgb(grid_surface_context, r, g, b);
      cairo_rectangle(grid_surface_context, scale * (float(ix)),
                      scale * ((float)iy), scale, scale);
      cairo_fill(grid_surface_context);
    }
  }

  cairo_set_source_rgb(grid_surface_context, 0.8, 0.0, 0);
  for (auto& scan : initial_pose_estimate_range_data.returns) {
    float x = scale * (limits.max().x() - scan[0]);
    float y = scale * (limits.max().y() - scan[1]);
    cairo_rectangle(grid_surface_context, (x - 0.15) * scale,
                    (y - 0.15) * scale, 0.3 * scale, 0.3 * scale);
  }
  cairo_fill(grid_surface_context);

  cairo_set_source_rgb(grid_surface_context, 0.0, 0.8, 0);
  for (auto& scan : matched_range_data.returns) {
    float x = scale * (limits.max().x() - scan[0]);
    float y = scale * (limits.max().y() - scan[1]);
    cairo_rectangle(grid_surface_context, (x - 0.15) * scale,
                    (y - 0.15) * scale, 0.3 * scale, 0.3 * scale);
  }

  cairo_fill(grid_surface_context);

  time_t seconds;
  time(&seconds);
  std::string filename = "grid_with_inserted_cloud" + std::to_string(seconds) +
                         std::to_string(rendered_grid_id) + ".png";
  rendered_grid_id++;
  cairo_surface_write_to_png(grid_surface, filename.c_str());
}

template <typename GridType, typename RangeDataInserter>
void EvaluateScanMatcher(
    const std::vector<Sample>& training_set,
    const std::vector<Sample>& test_set,
    const cartographer::mapping::proto::RangeDataInserterOptions2D&
        range_data_inserter_options,
    const cartographer::mapping::scan_matching::proto::
        CeresScanMatcherOptions2D& ceres_scan_matcher_options,
    std::vector<SampleResult>* results) {
  RangeDataInserter range_data_inserter(range_data_inserter_options);

  std::unique_ptr<GridType> grid = generateGrid<GridType>();

  for (auto sample : training_set) {
    range_data_inserter.Insert(sample.range_data, grid.get());
  }

  for (auto sample : test_set) {
    SampleResult sample_result;
    MatchScan(sample, ceres_scan_matcher_options, *grid.get(), &sample_result);
    results->push_back(sample_result);
  }
}

template <typename GridType>
void MatchScan(const Sample& sample,
               const cartographer::mapping::scan_matching::proto::
                   CeresScanMatcherOptions2D& ceres_scan_matcher_options,
               const GridType& grid, SampleResult* sample_result) {
  cartographer::mapping::scan_matching::CeresScanMatcher2D scan_matcher(
      ceres_scan_matcher_options);

  const Eigen::Vector2d target_translation = {0., 0.};
  const cartographer::transform::Rigid2d initial_pose_estimate =
      cartographer::transform::Rigid2d::Translation({0.0, 0.0});
  cartographer::transform::Rigid2d matched_pose_estimate;
  ceres::Solver::Summary summary;

  scan_matcher.Match(target_translation, initial_pose_estimate,
                     sample.range_data.returns, grid, &matched_pose_estimate,
                     &summary);
  const auto initial_error =
      initial_pose_estimate * sample.ground_truth_pose;
  const auto matching_error =
      matched_pose_estimate * sample.ground_truth_pose;
  sample_result->initial_trans_error = initial_error.translation().norm();
  sample_result->matching_trans_error = matching_error.translation().norm();
  sample_result->initial_rot_error = initial_error.rotation().smallestAngle();
  sample_result->matching_rot_error = matching_error.rotation().smallestAngle();
  sample_result->matching_iterations =
      summary.num_successful_steps + summary.num_unsuccessful_steps;
  sample_result->matching_time = summary.minimizer_time_in_seconds;
  LOG(INFO) << "Matching error " << sample_result->matching_trans_error << " \t"
            << sample_result->initial_trans_error << " \t"
            << sample_result->matching_iterations << " \t"
            << sample_result->matching_time;

  renderGridwithScan(grid, sample, initial_pose_estimate,
                     matched_pose_estimate);
}

template <typename GridType, typename RangeDataInserter>
void EvaluateScanMatcherGradient(
    const std::vector<Sample>& training_set,
    const std::vector<Sample>& test_set,
    const cartographer::mapping::proto::RangeDataInserterOptions2D&
        range_data_inserter_options,
    const cartographer::mapping::scan_matching::proto::
        CeresScanMatcherOptions2D& ceres_scan_matcher_options,
    std::vector<std::vector<double>>* gradients) {
  RangeDataInserter range_data_inserter(range_data_inserter_options);

  std::unique_ptr<GridType> grid = generateGrid<GridType>();

  for (auto sample : training_set) {
    range_data_inserter.Insert(sample.range_data, grid.get());
  }
  float min_trans = -1.8;
  float max_trans = 1.8;
  float resolution_trans = 0.025;
  float min_rot = 0.f;
  float max_rot = 2.f * M_PI;
  float resolution_rot = resolution_trans * (max_rot-min_rot)/(max_trans-min_trans);
  for (auto sample : test_set) {
    for (float x = min_trans; x < max_trans; x += resolution_trans) {
      for (float y = min_trans; y < max_trans; y += resolution_trans) {
        for (float theta = min_rot; theta < max_rot; theta += resolution_rot) {
          cartographer::mapping::scan_matching::CeresScanMatcher2D scan_matcher(
              ceres_scan_matcher_options);
          const Eigen::Vector2d target_translation = {0., 0.};
          const cartographer::transform::Rigid2d initial_pose_estimate =
              cartographer::transform::Rigid2d({x, y}, theta);
          cartographer::transform::Rigid2d matched_pose_estimate;
          ceres::Solver::Summary summary;
          double cost;
          std::vector<double> jacobians;
          scan_matcher.Evaluate(target_translation, initial_pose_estimate,
                                sample.range_data.returns, *grid.get(), &cost,
                                NULL, &jacobians);
          std::vector<double> result = {x, y, theta, cost};
          result.insert(result.end(), jacobians.begin(), jacobians.end());
          gradients->push_back(result);
        }
      }
      LOG(INFO)<<100.0*(1.0 -(max_trans - x)/(max_trans - min_trans))<<"%";
    }
  }
}

void RunScanMatchingEvaluation() {
  cartographer::mapping::proto::RangeDataInserterOptions2D
      range_data_inserter_options;
  auto parameter_dictionary_range_data_inserter = common::MakeDictionary(
      "return { "
      "probability_grid = {"
      "insert_free_space = true, "
      "hit_probability = 0.7, "
      "miss_probability = 0.4, "
      "},"
      "tsdf = {"
      "range_data_inserter_type = \"CONSTANT_WEIGHT\","
      "truncation_distance = 0.3,"
      "behind_surface_distance = 0.3,"
      "update_weight = 1.0,"
      "maximum_weight = 50.,"
      "},"
      "}");
  range_data_inserter_options =
      cartographer::mapping::CreateRangeDataInserterOptions2D(
          parameter_dictionary_range_data_inserter.get());
  auto parameter_dictionary = common::MakeDictionary(R"text(
        return {
          occupied_space_weight = 1.,
          translation_weight = 0.0,
          rotation_weight = 0.0,
          ceres_solver_options = {
            use_nonmonotonic_steps = true,
            max_num_iterations = 50,
            num_threads = 1,
          },
        })text");
  const cartographer::mapping::scan_matching::proto::CeresScanMatcherOptions2D
      ceres_scan_matcher_options =
          cartographer::mapping::scan_matching::CreateCeresScanMatcherOptions2D(
              parameter_dictionary.get());
  int n_training = 25;
  int n_test = 10000;

  std::ofstream log_file;
  std::string log_file_path;
  time_t seconds;
  time(&seconds);
  log_file_path = "scan_matching_benchmark" + std::to_string(seconds) + ".csv";
  log_file.open(log_file_path);
  log_file << "grid_type,grid_resolution,initial_error_trans,initial_error_"
              "angle,matched_error_trans,matched_error_angle,solver_iterations,"
              "matching_time\n";

  // std::vector<double> trans_errors = {0.05, 0.1, 0.25, 0.5};
  std::vector<double> trans_errors = {0.0};
  std::vector<double> rot_errors = {0.2 * M_PI_4, 0.4 * M_PI_4, 0.6 * M_PI_4,
                                    0.8 * M_PI_4};

  for (double error_trans : trans_errors) {
    for (double error_rot : rot_errors) {
      const ScanCloudGenerator::ModelType model_type =
          ScanCloudGenerator::ModelType::RECTANGLE;
      const Eigen::Vector2d size = {0.5, 0.5};
      const double resolution = 0.05;
      std::vector<Sample> training_set;
      std::vector<Sample> test_set;
      GenerateSampleSet(n_training, n_test, error_trans, error_rot, model_type,
                        size, resolution, training_set, test_set);
      std::vector<SampleResult> probability_grid_results;
      LOG(INFO) << "Evaluating probability grid:";
      EvaluateScanMatcher<
          cartographer::mapping::ProbabilityGrid,
          cartographer::mapping::RangeDataInserter2DProbabilityGrid>(
          training_set, test_set, range_data_inserter_options,
          ceres_scan_matcher_options, &probability_grid_results);
      for (const auto& res : probability_grid_results) {
        log_file << "PROBABILITY_GRID"
                 << "," << resolution << "," << res.initial_trans_error << ","
                 << res.initial_rot_error << "," << res.matching_trans_error
                 << "," << res.matching_rot_error << ","
                 << res.matching_iterations << "," << res.matching_time << "\n";
      }
      std::vector<SampleResult> tsdf_results;
      LOG(INFO) << "Evaluating TSDF:";
      EvaluateScanMatcher<cartographer::mapping::TSDF2D,
                          cartographer::mapping::RangeDataInserter2DTSDF>(
          training_set, test_set, range_data_inserter_options,
          ceres_scan_matcher_options, &tsdf_results);
      for (const auto& res : tsdf_results) {
        log_file << "TSDF"
                 << "," << resolution << "," << res.initial_trans_error << ","
                 << res.initial_rot_error << "," << res.matching_trans_error
                 << "," << res.matching_rot_error << ","
                 << res.matching_iterations << "," << res.matching_time << "\n";
      }
    }
  }
  log_file.close();
}

void EvaluateGradients() {
  cartographer::mapping::proto::RangeDataInserterOptions2D
      range_data_inserter_options;
  auto parameter_dictionary_range_data_inserter = common::MakeDictionary(
      "return { "
      "probability_grid = {"
      "insert_free_space = true, "
      "hit_probability = 0.7, "
      "miss_probability = 0.4, "
      "},"
      "tsdf = {"
      "range_data_inserter_type = \"CONSTANT_WEIGHT\","
      "truncation_distance = 0.3,"
      "behind_surface_distance = 0.3,"
      "update_weight = 1.0,"
      "maximum_weight = 50.,"
      "},"
      "}");
  range_data_inserter_options =
      cartographer::mapping::CreateRangeDataInserterOptions2D(
          parameter_dictionary_range_data_inserter.get());
  auto parameter_dictionary = common::MakeDictionary(R"text(
        return {
          occupied_space_weight = 1.,
          translation_weight = 0.0,
          rotation_weight = 0.0,
          ceres_solver_options = {
            use_nonmonotonic_steps = true,
            max_num_iterations = 50,
            num_threads = 1,
          },
        })text");
  const cartographer::mapping::scan_matching::proto::CeresScanMatcherOptions2D
      ceres_scan_matcher_options =
          cartographer::mapping::scan_matching::CreateCeresScanMatcherOptions2D(
              parameter_dictionary.get());
  int n_training = 25;
  int n_test = 1;

  std::vector<double> trans_errors = {0.0};
  std::vector<double> rot_errors = {0.0};
  for (double error_trans : trans_errors) {
    for (double error_rot : rot_errors) {
      const ScanCloudGenerator::ModelType model_type =
          ScanCloudGenerator::ModelType::RECTANGLE;
      const Eigen::Vector2d size = {1.0, 1.0};
      const double resolution = 0.05;
      std::vector<Sample> training_set;
      std::vector<Sample> test_set;
      GenerateSampleSet(n_training, n_test, error_trans, error_rot, model_type,
                        size, resolution, training_set, test_set);
      std::vector<std::vector<double>> gradients;
      LOG(INFO) << "Evaluating probability grid:";
      EvaluateScanMatcherGradient<
          cartographer::mapping::ProbabilityGrid,
          cartographer::mapping::RangeDataInserter2DProbabilityGrid>(
          training_set, test_set, range_data_inserter_options,
          ceres_scan_matcher_options, &gradients);

      std::ofstream log_file;
      std::string log_file_path;
      time_t seconds;
      time(&seconds);
      log_file_path = "gradient_pg_" + std::to_string(seconds) + ".csv";
      log_file.open(log_file_path);
      for (auto& row : gradients) {
        for (auto& element : row) {
          log_file << element << ",";
        }
        log_file << "\n";
      }
      log_file.close();

      LOG(INFO) << "Evaluating TSDF:";
      gradients.clear();
      EvaluateScanMatcherGradient<
          cartographer::mapping::TSDF2D,
          cartographer::mapping::RangeDataInserter2DTSDF>(
          training_set, test_set, range_data_inserter_options,
          ceres_scan_matcher_options, &gradients);
      log_file_path = "gradient_tsdf_" + std::to_string(seconds) + ".csv";
      log_file.open(log_file_path);
      for (auto& row : gradients) {
        for (auto& element : row) {
          log_file << element << ",";
        }
        log_file << "\n";
      }
      log_file.close();

      std::vector<Sample> test_set_single_point;
      Sample single_point_sample;
      single_point_sample.range_data.returns.push_back({0.f, 0.f, 0.f});
      test_set_single_point.push_back(single_point_sample);

      LOG(INFO) << "Evaluating probability grid Single Point:";
      gradients.clear();
      EvaluateScanMatcherGradient<
          cartographer::mapping::ProbabilityGrid,
          cartographer::mapping::RangeDataInserter2DProbabilityGrid>(
          training_set, test_set_single_point, range_data_inserter_options,
          ceres_scan_matcher_options, &gradients);
      log_file_path =
          "gradient_pg_single_point_" + std::to_string(seconds) + ".csv";
      log_file.open(log_file_path);
      for (auto& row : gradients) {
        for (auto& element : row) {
          log_file << element << ",";
        }
        log_file << "\n";
      }
      log_file.close();

      LOG(INFO) << "Evaluating TSDF Single Point:";
      gradients.clear();
      EvaluateScanMatcherGradient<
          cartographer::mapping::TSDF2D,
          cartographer::mapping::RangeDataInserter2DTSDF>(
          training_set, test_set_single_point, range_data_inserter_options,
          ceres_scan_matcher_options, &gradients);
      log_file_path =
          "gradient_tsdf_single_point_" + std::to_string(seconds) + ".csv";
      log_file.open(log_file_path);
      for (auto& row : gradients) {
        for (auto& element : row) {
          log_file << element << ",";
        }
        log_file << "\n";
      }
      log_file.close();
    }
  }
}

}  // namespace evaluation
}  // namespace cartographer

int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = true;
  google::ParseCommandLineFlags(&argc, &argv, true);
  //cartographer::evaluation::RunScanMatchingEvaluation();
  cartographer::evaluation::EvaluateGradients();
}
