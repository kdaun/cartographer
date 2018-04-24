#include "cartographer/evaluation/scan_cloud_generator.h"

#include <random>

namespace cartographer {
namespace evaluation {

static std::default_random_engine e1(42);

ScanCloudGenerator::ScanCloudGenerator() : resolution_(0.025) {}

ScanCloudGenerator::ScanCloudGenerator(float resolution)
    : resolution_(resolution) {}

void ScanCloudGenerator::generateSquare(cartographer::sensor::PointCloud& cloud,
                                        float size) {
  generateRectangle(cloud, size, size);
}

void ScanCloudGenerator::generateRectangle(
    cartographer::sensor::PointCloud& cloud, float size_x, float size_y) {
  // std::random_device r;
  //
  std::normal_distribution<float> normal_distribution(0, 0.01);

  cloud.clear();
  float x_min = -size_x / 2.0;
  float x_max = size_x / 2.0;
  float y_min = -size_y / 2.0;
  float y_max = size_y / 2.0;

  for (float x = x_min; x <= x_max; x += resolution_) {
    float y = y_min;
    cloud.emplace_back(Eigen::Vector3f(x + normal_distribution(e1),
                                       y + normal_distribution(e1), 0.));
    y = y_max;
    cloud.emplace_back(Eigen::Vector3f(x + normal_distribution(e1),
                                       y + normal_distribution(e1), 0.));
  }
  for (float y = y_min; y <= y_max; y += resolution_) {
    float x = x_min;
    cloud.emplace_back(Eigen::Vector3f(x + normal_distribution(e1),
                                       y + normal_distribution(e1), 0.));
    x = x_max;
    cloud.emplace_back(Eigen::Vector3f(x + normal_distribution(e1),
                                       y + normal_distribution(e1), 0.));
  }
}

void ScanCloudGenerator::generateCircle(cartographer::sensor::PointCloud& cloud,
                                        float radius) {
  // std::random_device r;
  // std::default_random_engine e1(42);
  std::normal_distribution<float> normal_distribution(0, 0.01);

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