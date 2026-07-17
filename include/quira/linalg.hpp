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

#pragma once

#include "quira/constants.hpp"
#include "quira/types.hpp"

namespace quira::linalg {

/**
 * @brief Returns the Hilbert-space dimension for an n-qubit state.
 *
 * @param num_qubits Number of qubits.
 * @return Dimension equal to 2^num_qubits.
 *
 * @throws std::out_of_range If the result cannot fit in Index.
 */
[[nodiscard]] std::size_t basis_size(std::size_t num_qubits);

/**
 * @brief Returns the conjugate transpose of a matrix.
 *
 * @param matrix Input matrix.
 * @return Adjoint matrix.
 */
[[nodiscard]] types::c_mat adjoint(const types::c_mat& matrix);

/**
 * @brief Returns the trace of a square matrix.
 *
 * @param matrix Input square matrix.
 * @return Sum of the diagonal entries.
 *
 * @throws std::invalid_argument If matrix is empty or not square.
 */
[[nodiscard]] types::cplx_n trace(const types::c_mat& matrix);

/**
 * @brief Returns the Euclidean norm of a vector.
 *
 * @param vector Input vector.
 * @return Vector norm.
 *
 * @throws std::invalid_argument If vector is empty.
 */
[[nodiscard]] types::real_n norm(const types::ket& vector);

/**
 * @brief Returns the Frobenius norm of a matrix.
 *
 * @param matrix Input matrix.
 * @return Matrix Frobenius norm.
 *
 * @throws std::invalid_argument If matrix is empty.
 */
[[nodiscard]] double norm(const types::c_mat& matrix);

/**
 * @brief Returns the kronecker product (tensor product) of two matrices.
 *
 * @param lhs Left matrix.
 * @param rhs Right matrix.
 * @return Kronecker product lhs tensor rhs.
 *
 * @throws std::invalid_argument If either matrix has zero size.
 */
[[nodiscard]] types::c_mat kron(const types::c_mat& lhs, const types::c_mat& rhs);

/**
 * @brief Returns the kronecker product (tensor product) of two vectors.
 *
 * @param lhs Left vector.
 * @param rhs Right vector.
 * @return Kronecker product lhs tensor rhs.
 *
 * @throws std::invalid_argument If either vector has zero size.
 */
[[nodiscard]] types::ket kron(const types::ket& lhs, const types::ket& rhs);

/**
 * @brief Returns repeated kronecker product (tensor product) of a matrix.
 *
 * @param matrix Matrix to repeat.
 * @param power Number of tensor factors.
 * @return matrix tensor matrix tensor ... tensor matrix.
 *
 * @throws std::invalid_argument If matrix has zero size or count is zero.
 */
[[nodiscard]] types::c_mat kron_pow(const types::c_mat& matrix, std::size_t power);

/**
 * @brief Returns repeated kronecker product (tensor product) of a vector.
 *
 * @param vector Vector to repeat.
 * @param power Number of tensor factors.
 * @return vector tensor vector tensor ... tensor vector.
 *
 * @throws std::invalid_argument If vector has zero size or count is zero.
 */
[[nodiscard]] types::ket kron_pow(const types::ket& vector, std::size_t power);

/**
 * @brief Returns the projector |psi><psi|.
 *
 * @param state State vector.
 * @return Normalized rank-one projector.
 *
 * @throws std::invalid_argument If state has zero size.
 * @throws std::runtime_error If state has zero norm.
 */
[[nodiscard]] types::c_mat projector(const types::ket& state);

/**
 * @brief Checks whether a matrix is unitary.
 *
 * @param matrix Matrix to check.
 * @param tolerance Numerical tolerance.
 * @return True if matrix^dagger * matrix is approximately identity.
 *
 * @throws std::invalid_argument If tolerance is negative or NaN.
 */
[[nodiscard]] bool is_unitary(const types::c_mat& matrix,
                              double tolerance = constants::EPS);

/**
 * @brief Checks whether a matrix is Hermitian.
 *
 * @param matrix Matrix to check.
 * @param tolerance Numerical tolerance.
 * @return True if matrix is approximately equal to its adjoint.
 *
 * @throws std::invalid_argument If tolerance is negative or NaN.
 */
[[nodiscard]] bool is_hermitian(const types::c_mat& matrix,
                                double tolerance = constants::EPS);

/**
 * @brief Checks whether a state vector is normalized.
 *
 * @param state State vector.
 * @param tolerance Numerical tolerance.
 * @return True if the squared norm is approximately one.
 */
[[nodiscard]] bool is_normalized(const types::ket& state,
                                 double tolerance = constants::EPS);

/**
 * @brief Returns a normalized copy of a state vector.
 *
 * @param state State vector.
 * @return Normalized state vector.
 *
 * @throws std::invalid_argument If state has zero size.
 * @throws std::runtime_error If state has zero norm.
 */
[[nodiscard]] types::ket normalize(const types::ket& state);

/**
 * @brief Returns the matrix commutator [A, B] = AB - BA.
 *
 * @param lhs Left square matrix.
 * @param rhs Right square matrix.
 * @return Matrix commutator.
 *
 * @throws std::invalid_argument If matrices are empty, not square, or have
 * incompatible dimensions.
 */
[[nodiscard]] types::c_mat commutator(const types::c_mat& lhs, const types::c_mat& rhs);

/**
 * @brief Returns the eigenvalues of a square matrix.
 *
 * @param matrix Input square matrix.
 * @return Complex eigenvalues.
 *
 * @throws std::invalid_argument If matrix is empty or not square.
 * @throws std::runtime_error If Eigen cannot compute the decomposition.
 */
[[nodiscard]] types::ket eigenvalues(const types::c_mat& matrix);

/**
 * @brief Returns the singular values of a matrix.
 *
 * @param matrix Input matrix.
 * @return Singular values in descending order.
 *
 * @throws std::invalid_argument If matrix is empty.
 */
[[nodiscard]] std::vector<double> singular_values(const types::c_mat& matrix);

}  // namespace quira::linalg
