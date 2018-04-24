#include <cartographer/mapping/2d/scan_matching/proto/ceres_scan_matcher_options_2d.pb.h>

#include <random>

#include "cairo.h"

#include "cartographer/common/lua_parameter_dictionary.h"
#include "cartographer/common/lua_parameter_dictionary_test_helpers.h"
#include "cartographer/common/make_unique.h"
#include "cartographer/evaluation/scan_cloud_generator.h"
#include "cartographer/mapping/2d/probability_grid.h"
#include "cartographer/mapping/2d/range_data_inserter_2d.h"
#include "cartographer/mapping/2d/range_data_inserter_2d_probability_grid.h"
#include "cartographer/mapping/2d/range_data_inserter_2d_tsdf.h"
#include "cartographer/mapping/internal/2d/scan_matching/ceres_scan_matcher_2d.h"
#include "cartographer/mapping/internal/2d/scan_matching/occupied_space_cost_function_2d.h"

namespace cartographer {
namespace evaluation {

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
static std::random_device rd;
static std::default_random_engine e1(rd());

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

Sample GenerateSample(double error_trans,
                      const ScanCloudGenerator::ModelType model_type,
                      const Eigen::Vector2d size, const double resolution) {
  cartographer::sensor::RangeData range_data;
  GenerateRangeData(model_type, size, resolution, range_data);

  std::uniform_real_distribution<double> error_distribution_translation(
      error_trans - resolution / 2., error_trans - resolution / 2.);
  std::uniform_real_distribution<double> error_translation_direction(-M_PI,
                                                                     M_PI);
  double orientation = error_translation_direction(e1);
  double scale = error_trans == 0.0 ? 0.0 : error_distribution_translation(e1);
  double x = std::cos(orientation) * scale;
  double y = std::sin(orientation) * scale;
  Sample sample;
  sample.ground_truth_pose =
      cartographer::transform::Rigid2d::Translation({x, y});
  sensor::RangeData initial_pose_estimate_range_data =
      cartographer::sensor::TransformRangeData(
          range_data,
          transform::Embed3D(sample.ground_truth_pose.cast<float>()));
  sample.range_data = range_data;
  return sample;
}

void GenerateSampleSet(int n_training, int n_test, double error_trans,
                       const ScanCloudGenerator::ModelType model_type,
                       const Eigen::Vector2d size, const double resolution,
                       std::vector<Sample>& training_set,
                       std::vector<Sample>& test_set) {
  for (int i = 0; i < n_training; ++i) {
    training_set.push_back(GenerateSample(0.0, model_type, size, resolution));
  }
  for (int i = 0; i < n_test; ++i) {
    test_set.push_back(
        GenerateSample(error_trans, model_type, size, resolution));
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

template <typename GridType, typename RangeDataInserter>
void EvaluateScanMatcher(
    const std::vector<Sample>& training_set,
    const std::vector<Sample>& test_set,
    const cartographer::mapping::proto::RangeDataInserterOptions2D&
        range_data_inserter_options,
    const cartographer::mapping::scan_matching::proto::
        CeresScanMatcherOptions2D& ceres_scan_matcher_options,
    std::vector<SampleResult> results) {
  RangeDataInserter range_data_inserter(range_data_inserter_options);

  std::unique_ptr<GridType> grid = generateGrid<GridType>();

  for (auto sample : training_set) {
    range_data_inserter.Insert(sample.range_data, grid.get());
  }

  for (auto sample : test_set) {
    SampleResult sample_result;
    MatchScan(sample, ceres_scan_matcher_options, *grid.get(), &sample_result);
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
      initial_pose_estimate * sample.ground_truth_pose.inverse();
  const auto matching_error =
      matched_pose_estimate * sample.ground_truth_pose.inverse();
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
  int n_training = 5;
  int n_test = 5;
  double error_trans = 0.05;
  const ScanCloudGenerator::ModelType model_type =
      ScanCloudGenerator::ModelType::RECTANGLE;
  const Eigen::Vector2d size = {0.6, 0.6};
  const double resolution = 0.05;
  std::vector<Sample> training_set;
  std::vector<Sample> test_set;
  GenerateSampleSet(n_training, n_test, error_trans, model_type, size,
                    resolution, training_set, test_set);
  std::vector<SampleResult> probability_grid_results;
  LOG(INFO) << "Evaluating probability grid:";
  EvaluateScanMatcher<
      cartographer::mapping::ProbabilityGrid,
      cartographer::mapping::RangeDataInserter2DProbabilityGrid>(
      training_set, test_set, range_data_inserter_options,
      ceres_scan_matcher_options, probability_grid_results);
  LOG(INFO) << "Evaluating TSDF:";
  EvaluateScanMatcher<cartographer::mapping::TSDF2D,
                      cartographer::mapping::RangeDataInserter2DTSDF>(
      training_set, test_set, range_data_inserter_options,
      ceres_scan_matcher_options, probability_grid_results);

  /*
  sensor::RangeData initial_pose_estimate_range_data =
      cartographer::sensor::TransformRangeData(
          range_data, transform::Embed3D(initial_pose_estimate.cast<float>()));
  sensor::RangeData matched_range_data =
      cartographer::sensor::TransformRangeData(
          range_data, transform::Embed3D(matched_pose_estimate.cast<float>()));

  const cartographer::mapping::MapLimits& limits = probability_grid.limits();
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
      float p = 1. - probability_grid.GetProbability({iy, ix});
      cairo_set_source_rgb(grid_surface_context, p, p, p);
      cairo_rectangle(grid_surface_context, scale * (float(ix)),
                      scale * ((float)iy), scale, scale);
      cairo_fill(grid_surface_context);
    }
  }

  cairo_set_source_rgb(grid_surface_context, 0.0, 0.0, 0.8);
  for (auto& scan : scan_cloud) {
    float x = scale * (limits.max().x() - scan[0]);
    float y = scale * (limits.max().y() - scan[1]);
    cairo_rectangle(grid_surface_context, (x - 0.15) * scale,
                    (y - 0.15) * scale, 0.3 * scale, 0.3 * scale);
  }
  cairo_fill(grid_surface_context);

  LOG(INFO) << "grid: "
            << cairo_surface_write_to_png(grid_surface,
                                          "grid_with_inserted_cloud.png");

  const cartographer::mapping::scan_matching::GridArrayAdapter adapter(
      probability_grid);
  ceres::BiCubicInterpolator<
      cartographer::mapping::scan_matching::GridArrayAdapter>
      interpolator(adapter);

  cairo_surface_t* interpolated_grid_surface;
  cairo_t* interpolated_grid_surface_context;
  interpolated_grid_surface = cairo_image_surface_create(
      CAIRO_FORMAT_ARGB32, scaled_num_x_cells, scaled_num_y_cells);
  interpolated_grid_surface_context = cairo_create(interpolated_grid_surface);

  for (double ix = 0; ix < scaled_num_x_cells; ix += 0.2) {
    for (double iy = 0; iy < scaled_num_y_cells; iy += 0.2) {
      double interpolated_value;
      interpolator.Evaluate(ix + double(INT_MAX / 4) - 0.5,
                            iy + double(INT_MAX / 4) - 0.5,
                            &interpolated_value);
      float p = interpolated_value;
      cairo_set_source_rgb(interpolated_grid_surface_context, p, p, p);
      cairo_rectangle(interpolated_grid_surface_context, scale * (float(ix)),
                      scale * ((float)iy), scale, scale);
      cairo_fill(interpolated_grid_surface_context);
    }
  }

  cairo_set_source_rgb(interpolated_grid_surface_context, 0.8, 0.0, 0);
  for (auto& scan : initial_pose_estimate_range_data.returns) {
    float x = scale * (limits.max().x() - scan[0]);
    float y = scale * (limits.max().y() - scan[1]);
    cairo_rectangle(interpolated_grid_surface_context, (x - 0.15) * scale,
                    (y - 0.15) * scale, 0.3 * scale, 0.3 * scale);
  }
  cairo_fill(interpolated_grid_surface_context);

  cairo_set_source_rgb(interpolated_grid_surface_context, 0.0, 0.8, 0);
  for (auto& scan : matched_range_data.returns) {
    float x = scale * (limits.max().x() - scan[0]);
    float y = scale * (limits.max().y() - scan[1]);
    cairo_rectangle(interpolated_grid_surface_context, (x - 0.15) * scale,
                    (y - 0.15) * scale, 0.3 * scale, 0.3 * scale);
  }
  cairo_fill(interpolated_grid_surface_context);

  LOG(INFO) << "interpolated_grid: "
            << cairo_surface_write_to_png(interpolated_grid_surface,
                                          "interpolated_grid_with.png");
                                          */
}

}  // namespace evaluation
}  // namespace cartographer

int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = true;
  google::ParseCommandLineFlags(&argc, &argv, true);
  cartographer::evaluation::RunScanMatchingEvaluation();
}
