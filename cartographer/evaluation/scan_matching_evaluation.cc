#include <cartographer/mapping/2d/scan_matching/proto/ceres_scan_matcher_options_2d.pb.h>

#include "cartographer/common/lua_parameter_dictionary.h"
#include "cartographer/common/lua_parameter_dictionary_test_helpers.h"
#include "cartographer/common/make_unique.h"
#include "cartographer/evaluation/scan_cloud_generator.h"
#include "cartographer/mapping/2d/range_data_inserter_2d.h"
#include "cartographer/mapping/2d/range_data_inserter_2d_probability_grid.h"
#include "cartographer/mapping/2d/range_data_inserter_2d_tsdf.h"
#include "cartographer/mapping/internal/2d/scan_matching/ceres_scan_matcher_2d.h"
#include "cartographer/mapping/internal/2d/scan_matching/occupied_space_cost_function_2d.h"

#include "cartographer/mapping/2d/probability_grid.h"

#include "cairo.h"

namespace cartographer {
namespace evaluation {

void RunScanMatchingEvaluation() {
  cartographer::sensor::PointCloud scan_cloud;
  ScanCloudGenerator test_set_generator(0.01);
  test_set_generator.generateRectangle(scan_cloud, 1.0, 0.5);

  const cartographer::transform::Rigid2d initial_pose_estimate =
      cartographer::transform::Rigid2d::Translation({0.05, 0.05});
  const Eigen::Vector2d target_translation = {0., 0.};
  cartographer::transform::Rigid2d matched_pose_estimate;
  ceres::Solver::Summary summary;
  cartographer::mapping::ProbabilityGrid probability_grid(
      cartographer::mapping::MapLimits(
          0.05, Eigen::Vector2d(1., 1.),
          cartographer::mapping::CellLimits(40, 40)));
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
  cartographer::mapping::RangeDataInserter2DProbabilityGrid range_data_inserter(
      range_data_inserter_options);

  cartographer::sensor::RangeData range_data;
  range_data.returns = scan_cloud;
  range_data.origin = Eigen::Vector3f{0, 0, 0};

  for (int i = 0; i < 10; ++i) {
    range_data_inserter.Insert(range_data, &probability_grid);
  }
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

  cartographer::mapping::scan_matching::CeresScanMatcher2D scan_matcher(
      ceres_scan_matcher_options);
  scan_matcher.Match(target_translation, initial_pose_estimate, scan_cloud,
                     probability_grid, &matched_pose_estimate, &summary);
  LOG(INFO) << summary.FullReport();
  LOG(INFO) << "matched_pose_estimate " << matched_pose_estimate;

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
}

}  // namespace evaluation
}  // namespace cartographer

int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = true;
  google::ParseCommandLineFlags(&argc, &argv, true);
  cartographer::evaluation::RunScanMatchingEvaluation();
}
