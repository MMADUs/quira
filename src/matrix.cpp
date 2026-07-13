#include <quira/matrix.hpp>

namespace quira {

Eigen::Matrix2d add_matrices(const Eigen::Matrix2d& a, const Eigen::Matrix2d& b) {
  return a + b;
}

}  // namespace quira
