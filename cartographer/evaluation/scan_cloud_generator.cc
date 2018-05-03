#include "cartographer/evaluation/scan_cloud_generator.h"

#include <random>

namespace cartographer {
namespace evaluation {
namespace {
static std::default_random_engine e1(42);
}  // namespace

ScanCloudGenerator::ScanCloudGenerator() : resolution_(0.02) {}

ScanCloudGenerator::ScanCloudGenerator(float resolution)
    : resolution_(resolution) {}

void ScanCloudGenerator::generateSquare(cartographer::sensor::PointCloud& cloud,
                                        float size, float noise_std_dev) {
  generateRectangle(cloud, size, size, noise_std_dev);
}

void ScanCloudGenerator::generateRectangle(
    cartographer::sensor::PointCloud& cloud, float size_x, float size_y,
    float noise_std_dev) {
  // std::random_device r;
  //
  std::normal_distribution<float> normal_distribution(0,
                                                      noise_std_dev);  // 0.01

  cloud.clear();
  float x_min = -size_x / 2.0;
  float x_max = size_x / 2.0;
  float y_min = -size_y / 2.0;
  float y_max = size_y / 2.0;

  std::uniform_real_distribution<double> error_translation_direction(-M_PI,
                                                                     M_PI);
  double e_orientation;
  double e_scale = noise_std_dev == 0.0 ? 0.0 : noise_std_dev;
  double e_x;
  double e_y;

  for (float x = x_min; x <= x_max + 1e-5; x += resolution_) {
    float y = y_min;
    e_orientation = error_translation_direction(e1);
    e_x = std::cos(e_orientation) * e_scale;
    e_y = std::sin(e_orientation) * e_scale;
    cloud.emplace_back(Eigen::Vector3f(x + e_x, y + e_y, 0.));
    y = y_max;
    e_orientation = error_translation_direction(e1);
    e_x = std::cos(e_orientation) * e_scale;
    e_y = std::sin(e_orientation) * e_scale;
    cloud.emplace_back(Eigen::Vector3f(x + e_x, y + e_y, 0.));
  }
  for (float y = y_min; y <= y_max + 1e-5; y += resolution_) {
    float x = x_min;
    e_orientation = error_translation_direction(e1);
    e_x = std::cos(e_orientation) * e_scale;
    e_y = std::sin(e_orientation) * e_scale;
    cloud.emplace_back(Eigen::Vector3f(x + e_x, y + e_y, 0.));
    x = x_max;
    e_orientation = error_translation_direction(e1);
    e_x = std::cos(e_orientation) * e_scale;
    e_y = std::sin(e_orientation) * e_scale;
    cloud.emplace_back(Eigen::Vector3f(x + e_x, y + e_y, 0.));
  }
}

void ScanCloudGenerator::generateCircle(cartographer::sensor::PointCloud& cloud,
                                        float radius, float noise_std_dev) {
  // std::random_device r;
  // std::default_random_engine e1(42);
  std::normal_distribution<float> normal_distribution(0, noise_std_dev);

  cloud.clear();
  float angular_resolution = 2.0 * std::asin(resolution_ / (2 * radius));

  for (float angle = 0; angle < 2.0 * M_PI; angle += angular_resolution) {
    float x = std::cos(angle) * radius;
    float y = std::sin(angle) * radius;
    cloud.emplace_back(Eigen::Vector3f(x + normal_distribution(e1),
                                       y + normal_distribution(e1), 0.));
  }
}

}  // namespace evaluation
}  // namespace cartographer