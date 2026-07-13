#include <quira/matrix.hpp>

#include <Eigen/Dense>
#include <cassert>

int main() {
    Eigen::Matrix2d a;
    a << 1.0, 2.0,
         3.0, 4.0;

    Eigen::Matrix2d b;
    b << 5.0, 6.0,
         7.0, 8.0;

    Eigen::Matrix2d expected;
    expected << 6.0, 8.0,
                10.0, 12.0;

    Eigen::Matrix2d result = quira::add_matrices(a, b);

    assert(result.isApprox(expected));

    return 0;
}
