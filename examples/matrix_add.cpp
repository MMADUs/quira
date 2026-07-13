#include <quira/matrix.hpp>

#include <Eigen/Dense>
#include <iostream>

int main() {
    Eigen::Matrix2d a;
    a << 1.0, 2.0,
         3.0, 4.0;

    Eigen::Matrix2d b;
    b << 5.0, 6.0,
         7.0, 8.0;

    Eigen::Matrix2d c = quira::add_matrices(a, b);

    std::cout << "Matrix A:\n" << a << "\n\n";
    std::cout << "Matrix B:\n" << b << "\n\n";
    std::cout << "A + B:\n" << c << '\n';

    return 0;
}
