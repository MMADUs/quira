// Copyright 2026 Muhammad Nizwa
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "quira/linalg.hpp"

#include "quira/exception.hpp"

#include <cmath>
#include <limits>
#include <sstream>

namespace quira::linalg {

namespace {

void validate_tolerance(types::real_n tolerance, const char* caller) {
  if (!(tolerance >= types::real_n{0.0})) {
    std::ostringstream oss;
    oss << " requires a non-negative tolerance, got " << tolerance << " instead";

    throw exception::InvalidArgument(caller, oss.str());
  }
}

void validate_non_empty_matrix(const types::c_mat& matrix, const char* caller) {
  if (matrix.size() == 0) {
    std::ostringstream oss;
    oss << "requires a non-empty matrix, got " << matrix.size() << " instead";

    throw exception::InvalidArgument(caller, oss.str());
  }
}

void validate_non_empty_vector(const types::ket& vector, const char* caller) {
  if (vector.size() == 0) {
    std::ostringstream oss;
    oss << "requires a non-empty vector, got " << vector.size() << " instead";

    throw exception::InvalidArgument(caller, oss.str());
  }
}

void validate_square_matrix(const types::c_mat& matrix, const char* caller) {
  validate_non_empty_matrix(matrix, caller);

  if (matrix.rows() != matrix.cols()) {
    std::ostringstream oss;
    oss << "requires a square matrix, got " << "(rows= " << matrix.rows()
        << ", cols=" << matrix.cols() << ") instead";

    throw exception::InvalidArgument(caller, oss.str());
  }
}

const char* eigen_info_message(Eigen::ComputationInfo info) {
  switch (info) {
    case Eigen::Success:
      return "success";
    case Eigen::NumericalIssue:
      return "numerical issue";
    case Eigen::NoConvergence:
      return "no convergence";
    case Eigen::InvalidInput:
      return "invalid input";
  }

  return "unknown Eigen computation error";
}

}  // namespace

std::size_t basis_size(std::size_t num_qubits) {
  if (num_qubits >= std::numeric_limits<std::size_t>::digits) {
    std::ostringstream oss;
    oss << "cannot represent 2^" << num_qubits
        << " as std::size_t; maximum supported qubit count is "
        << (std::numeric_limits<std::size_t>::digits - 1);

    throw exception::OutOfRange("linalg::basis_size()", oss.str());
  }

  return std::size_t{1} << num_qubits;
}

types::c_mat adjoint(const types::c_mat& matrix) {
  return matrix.adjoint();
}

types::cplx_n trace(const types::c_mat& matrix) {
  validate_square_matrix(matrix, "linalg::trace()");
  return matrix.trace();
}

types::real_n norm(const types::ket& vector) {
  validate_non_empty_vector(vector, "linalg::norm()");
  return vector.norm();
}

types::real_n norm(const types::c_mat& matrix) {
  validate_non_empty_matrix(matrix, "linalg::norm()");
  return matrix.norm();
}

types::c_mat kron(const types::c_mat& lhs, const types::c_mat& rhs) {
  validate_non_empty_matrix(lhs, "linalg::kron()");
  validate_non_empty_matrix(rhs, "linalg::kron()");

  types::c_mat result(lhs.rows() * rhs.rows(), lhs.cols() * rhs.cols());

  for (Eigen::Index row = 0; row < lhs.rows(); ++row) {
    for (Eigen::Index col = 0; col < lhs.cols(); ++col) {
      result.block(row * rhs.rows(), col * rhs.cols(), rhs.rows(), rhs.cols()) =
          lhs(row, col) * rhs;
    }
  }

  return result;
}

types::ket kron(const types::ket& lhs, const types::ket& rhs) {
  validate_non_empty_vector(lhs, "linalg::kron()");
  validate_non_empty_vector(rhs, "linalg::kron()");

  types::ket result(lhs.size() * rhs.size());

  for (Eigen::Index row = 0; row < lhs.size(); ++row) {
    result.segment(row * rhs.size(), rhs.size()) = lhs(row) * rhs;
  }

  return result;
}

types::c_mat kron_pow(const types::c_mat& matrix, std::size_t power) {
  validate_non_empty_matrix(matrix, "linalg::kron_pow()");

  if (power == 0) {
    std::ostringstream oss;
    oss << "requires power greater than 0, got " << power << " instead";

    throw exception::InvalidArgument("linalg::kron_pow()", oss.str());
  }

  types::c_mat result = matrix;

  for (std::size_t i = 1; i < power; ++i) {
    result = kron(result, matrix);
  }

  return result;
}

types::ket kron_pow(const types::ket& vector, std::size_t power) {
  validate_non_empty_vector(vector, "linalg::kron_pow()");

  if (power == 0) {
    std::ostringstream oss;
    oss << "requires power greater than 0, got " << power << " instead";

    throw exception::InvalidArgument("linalg::kron_pow()", oss.str());
  }

  types::ket result = vector;

  for (std::size_t i = 1; i < power; ++i) {
    result = kron(result, vector);
  }

  return result;
}

types::c_mat projector(const types::ket& state) {
  validate_non_empty_vector(state, "linalg::projector()");

  const types::real_n norm_squared = state.squaredNorm();

  if (norm_squared == types::real_n{0.0}) {
    std::ostringstream oss;
    oss << "cannot build projector from a zero state vector";

    throw exception::InvalidState("linalg::projector()", oss.str());
  }

  return (state * state.adjoint()) / norm_squared;
}

bool is_unitary(const types::c_mat& matrix, types::real_n tolerance) {
  validate_tolerance(tolerance, "linalg::is_unitary()");

  if (matrix.rows() == 0 || matrix.cols() == 0) {
    return false;
  }

  if (matrix.rows() != matrix.cols()) {
    return false;
  }

  const types::c_mat identity = types::c_mat::Identity(matrix.rows(), matrix.cols());
  const types::c_mat product = matrix.adjoint() * matrix;

  return product.isApprox(identity, tolerance);
}

bool is_hermitian(const types::c_mat& matrix, types::real_n tolerance) {
  validate_tolerance(tolerance, "linalg::is_hermitian()");

  if (matrix.rows() == 0 || matrix.cols() == 0) {
    return false;
  }

  if (matrix.rows() != matrix.cols()) {
    return false;
  }

  return matrix.isApprox(matrix.adjoint(), tolerance);
}

bool is_normalized(const types::ket& state, types::real_n tolerance) {
  validate_tolerance(tolerance, "linalg::is_normalized()");

  if (state.size() == 0) {
    return false;
  }

  return std::abs(state.squaredNorm() - 1.0) <= tolerance;
}

types::ket normalize(const types::ket& state) {
  validate_non_empty_vector(state, "linalg::normalize()");

  const types::real_n current_norm = state.norm();

  if (current_norm == types::real_n{0.0}) {
    std::ostringstream oss;
    oss << "cannot normalize a zero state vector";

    throw exception::InvalidState("linalg::normalize()", oss.str());
  }

  return state / current_norm;
}

types::c_mat commutator(const types::c_mat& lhs, const types::c_mat& rhs) {
  validate_square_matrix(lhs, "linalg::commutator()");
  validate_square_matrix(rhs, "linalg::commutator()");

  // checking only rows is enough because validate square matrix will check rows == cols
  if (lhs.rows() != rhs.rows()) {
    std::ostringstream oss;
    oss << "requires matrices with matching dimension, got "
        << "(lhs rows=" << lhs.rows() << ", lhs cols=" << lhs.cols()
        << ", rhs rows=" << rhs.rows() << ", rhs cols=" << rhs.cols() << ") instead";

    throw exception::DimensionMismatch("linalg::commutator()", oss.str());
  }

  return lhs * rhs - rhs * lhs;
}

types::ket eigenvalues(const types::c_mat& matrix) {
  validate_square_matrix(matrix, "linalg::eigenvalues()");

  Eigen::ComplexEigenSolver<types::c_mat> solver(matrix, false);

  if (solver.info() != Eigen::Success) {
    std::ostringstream oss;
    oss << "failed with " << eigen_info_message(solver.info())
        << " for matrix with shape " << "(rows=" << matrix.rows()
        << ", cols=" << matrix.cols() << ")";

    throw exception::NumericError("linalg::eigenvalues()", oss.str());
  }

  return solver.eigenvalues();
}

std::vector<types::real_n> singular_values(const types::c_mat& matrix) {
  validate_non_empty_matrix(matrix, "linalg::singular_values()");

  Eigen::JacobiSVD<types::c_mat> svd(matrix);
  const auto& values = svd.singularValues();

  std::vector<types::real_n> result;
  result.reserve(static_cast<std::size_t>(values.size()));

  for (Eigen::Index index = 0; index < values.size(); ++index) {
    result.push_back(values(index));
  }

  return result;
}

}  // namespace quira::linalg
