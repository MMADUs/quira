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

#include "quira/types.hpp"

#include <cstddef>
#include <vector>

namespace quira::statistics {

/**
 * @brief Creates a uniform probability distribution.
 *
 * @param size Number of outcomes.
 * @return Probability vector with every entry equal to 1 / size.
 *
 * @throws std::invalid_argument If size is zero.
 */
[[nodiscard]] std::vector<types::real_n> uniform_distribution(std::size_t size);

/**
 * @brief Computes the row marginal of a joint probability table.
 *
 * For a table P(X, Y), this returns P(X).
 *
 * @param probability_table Joint probability table with rows for X and columns
 * for Y.
 * @return Marginal probability vector over X.
 *
 * @throws std::invalid_argument If the table is empty, contains invalid
 * probabilities, or is not normalized.
 */
[[nodiscard]] std::vector<types::real_n> marginal_x(
    const types::r_mat& probability_table);

/**
 * @brief Computes the column marginal of a joint probability table.
 *
 * For a table P(X, Y), this returns P(Y).
 *
 * @param probability_table Joint probability table with rows for X and columns
 * for Y.
 * @return Marginal probability vector over Y.
 *
 * @throws std::invalid_argument If the table is empty, contains invalid
 * probabilities, or is not normalized.
 */
[[nodiscard]] std::vector<types::real_n> marginal_y(
    const types::r_mat& probability_table);

/**
 * @brief Computes the mean of values under a probability distribution.
 *
 * @param probabilities Probability vector.
 * @param values Random-variable values matching probabilities.
 * @return Expected value.
 *
 * @throws std::invalid_argument If inputs are empty, sizes differ, or
 * probabilities are invalid.
 */
[[nodiscard]] types::real_n mean(const std::vector<types::real_n>& probabilities,
                                 const std::vector<types::real_n>& values);

/**
 * @brief Computes the classical variance of values under a distribution.
 *
 * @param probabilities Probability vector.
 * @param values Random-variable values matching probabilities.
 * @return Variance.
 *
 * @throws std::invalid_argument If inputs are empty, sizes differ, or
 * probabilities are invalid.
 */
[[nodiscard]] types::real_n variance(const std::vector<types::real_n>& probabilities,
                                     const std::vector<types::real_n>& values);

/**
 * @brief Computes covariance from a joint probability distribution.
 *
 * @param joint_probabilities Joint probability table P(X, Y).
 * @param x_values Values associated with rows.
 * @param y_values Values associated with columns.
 * @return Covariance of X and Y.
 *
 * @throws std::invalid_argument If inputs are empty, dimensions mismatch, or
 * probabilities are invalid.
 */
[[nodiscard]] types::real_n covariance(const types::r_mat& joint_probabilities,
                                       const std::vector<types::real_n>& x_values,
                                       const std::vector<types::real_n>& y_values);

/**
 * @brief Computes normalized covariance from a joint distribution.
 *
 * @param joint_probabilities Joint probability table P(X, Y).
 * @param x_values Values associated with rows.
 * @param y_values Values associated with columns.
 * @return Correlation coefficient.
 *
 * @throws std::invalid_argument If inputs are invalid or either variable has
 * zero variance.
 */
[[nodiscard]] types::real_n correlation(const types::r_mat& joint_probabilities,
                                        const std::vector<types::real_n>& x_values,
                                        const std::vector<types::real_n>& y_values);

}  // namespace quira::statistics
