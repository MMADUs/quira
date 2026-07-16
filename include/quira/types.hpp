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

#include <Eigen/Dense>

#include <complex>
#include <cstdint>

namespace quira::types {

/**
 * @brief Bit index type
 */
#if defined(SMALL_BIT)
using qubit = std::uint8_t;
using clbit = std::uint8_t;
#elif defined(LARGE_BIT)
using qubit = std::size_t;
using clbit = std::size_t;
#else
using qubit = std::uint8_t;
using clbit = std::uint8_t;
#endif

/**
 * @brief Real number precision type
 */
#if defined(REAL_FLOAT)
using real_n = float;
#elif defined(REAL_DOUBLE)
using real_n = double;
#elif defined(REAL_LONG_DOUBLE)
using real_n = long double;
#else
using real_n = double;
#endif

/**
 * @brief Complex number type
 */
using cplx_n = std::complex<real_n>;

/**
 * @brief Dynamic Eigen matrix over the field specified by \a Scalar
 */
template<typename Scalar>
using dyn_mat = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

/**
 * @brief Dynamic Eigen column vector over the field specified by Scalar
 */
template<typename Scalar>
using dyn_col_vect = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

/**
 * @brief Dynamic Eigen row vector over the field specified by Scalar
 */
template<typename Scalar>
using dyn_row_vect = Eigen::Matrix<Scalar, 1, Eigen::Dynamic>;

/**
 * @brief Complex number with real number precision dynamic Eigen column vector
 */
using ket = dyn_col_vect<cplx_n>;

/**
 * @brief Complex number with real number precision dynamic Eigen row vector
 */
using bra = dyn_row_vect<cplx_n>;

/**
 * @brief Complex number with real number precision dynamic Eigen matrix
 */
using c_mat = dyn_mat<cplx_n>;

/**
 * @brief Real number precision dynamic Eigen matrix
 */
using r_mat = dyn_mat<real_n>;

/**
 * @brief Textual representation (Dirac notation) of a quantum state
 */
template<typename Scalar>
// NOLINTNEXTLINE(readability-identifier-naming)
struct dirac_t {
  std::vector<std::size_t> dims_rows{};  ///< row dimensions
  std::vector<std::size_t> dims_cols{};  ///< column dimensions
  std::vector<std::pair<Scalar, std::vector<std::size_t>>>
      states{};  ///< vector of (amplitude, dits)

  /**
   * @brief Equality operator
   *
   * @param rhs dirac_t object against which the equality is being tested
   * @return True if the dirac_t objects are equal (component-wise), false
   * otherwise
   */
  bool operator==(const dirac_t& rhs) const {
    return std::tie(dims_rows, dims_cols, states) ==
           std::tie(rhs.dims_rows, rhs.dims_cols, rhs.states);
  }

  /**
   * @brief Inequality operator
   *
   * @param rhs dirac_t object against which the inequality is being tested
   * @return True if the dirac_t objects are not equal (component-wise),
   * false otherwise
   */

  bool operator!=(const dirac_t& rhs) const { return !(*this == rhs); }
};

/**
 * @brief Variant type-matching utility for std::visit
 * @tparam Ts Type list
 */
template<class... Ts>
struct Overloaded : Ts... {
  using Ts::operator()...;
};

/**
 * @brief Template deduction rule
 */
template<class... Ts>
Overloaded(Ts...) -> Overloaded<Ts...>;

}  // namespace quira::types
