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

#include "quira/statistics.hpp"

#include "quira/constants.hpp"
#include "quira/exception.hpp"

#include <cmath>
#include <sstream>

namespace quira::statistics {
namespace {

constexpr types::real_n PROBABILITY_TOLERANCE = constants::EPS * 1000.0;

bool approximately_one(types::real_n value) {
  return std::abs(value - types::real_n{1.0}) <= PROBABILITY_TOLERANCE;
}

void validate_probability_value(types::real_n value, const char* caller,
                                const char* label) {
  if (!std::isfinite(value) || value < types::real_n{0.0}) {
    std::ostringstream oss;
    oss << label << " must be finite and non-negative";

    throw exception::InvalidArgument(caller, oss.str());
  }
}

void validate_probability_vector(const std::vector<types::real_n>& probabilities,
                                 const char* caller) {
  if (probabilities.empty()) {
    throw exception::InvalidArgument(caller, "probability vector cannot be empty");
  }

  types::real_n total = types::real_n{0.0};

  for (std::size_t index = 0; index < probabilities.size(); ++index) {
    std::ostringstream label;
    label << "probability at index " << index;
    validate_probability_value(probabilities[index], caller, label.str().c_str());
    total += probabilities[index];
  }

  if (!approximately_one(total)) {
    std::ostringstream oss;
    oss << "probabilities must sum to 1, got " << total;

    throw exception::InvalidArgument(caller, oss.str());
  }
}

void validate_values(const std::vector<types::real_n>& values, const char* caller,
                     const char* label) {
  if (values.empty()) {
    std::ostringstream oss;
    oss << label << " cannot be empty";

    throw exception::InvalidArgument(caller, oss.str());
  }

  for (std::size_t index = 0; index < values.size(); ++index) {
    if (!std::isfinite(values[index])) {
      std::ostringstream oss;
      oss << label << " at index " << index << " must be finite";

      throw exception::InvalidArgument(caller, oss.str());
    }
  }
}

void validate_probability_values_pair(const std::vector<types::real_n>& probabilities,
                                      const std::vector<types::real_n>& values,
                                      const char* caller) {
  validate_probability_vector(probabilities, caller);
  validate_values(values, caller, "values");

  if (probabilities.size() != values.size()) {
    std::ostringstream oss;
    oss << "probabilities and values must have matching sizes, got "
        << probabilities.size() << " and " << values.size();

    throw exception::DimensionMismatch(caller, oss.str());
  }
}

void validate_joint_probability_table(const types::r_mat& table, const char* caller) {
  if (table.rows() == 0 || table.cols() == 0) {
    throw exception::InvalidArgument(caller, "joint probability table cannot be empty");
  }

  types::real_n total = types::real_n{0.0};

  for (Eigen::Index row = 0; row < table.rows(); ++row) {
    for (Eigen::Index col = 0; col < table.cols(); ++col) {
      const types::real_n value = table(row, col);

      if (!std::isfinite(value) || value < types::real_n{0.0}) {
        std::ostringstream oss;
        oss << "joint probability at row " << row << ", column " << col
            << " must be finite and non-negative";

        throw exception::InvalidArgument(caller, oss.str());
      }

      total += value;
    }
  }

  if (!approximately_one(total)) {
    std::ostringstream oss;
    oss << "joint probabilities must sum to 1, got " << total;

    throw exception::InvalidArgument(caller, oss.str());
  }
}

void validate_joint_values(const types::r_mat& joint_probabilities,
                           const std::vector<types::real_n>& x_values,
                           const std::vector<types::real_n>& y_values,
                           const char* caller) {
  validate_joint_probability_table(joint_probabilities, caller);
  validate_values(x_values, caller, "x_values");
  validate_values(y_values, caller, "y_values");

  if (static_cast<Eigen::Index>(x_values.size()) != joint_probabilities.rows()) {
    std::ostringstream oss;
    oss << "x_values size must match joint probability rows, got " << x_values.size()
        << " and " << joint_probabilities.rows();

    throw exception::DimensionMismatch(caller, oss.str());
  }

  if (static_cast<Eigen::Index>(y_values.size()) != joint_probabilities.cols()) {
    std::ostringstream oss;
    oss << "y_values size must match joint probability columns, got " << y_values.size()
        << " and " << joint_probabilities.cols();

    throw exception::DimensionMismatch(caller, oss.str());
  }
}

}  // namespace

std::vector<types::real_n> uniform_distribution(std::size_t size) {
  if (size == 0) {
    throw exception::InvalidArgument("statistics::uniform_distribution()",
                                     "size must be greater than zero");
  }

  return std::vector<types::real_n>(
      size, types::real_n{1.0} / static_cast<types::real_n>(size));
}

std::vector<types::real_n> marginal_x(const types::r_mat& probability_table) {
  constexpr const char* CALLER = "statistics::marginal_x()";
  validate_joint_probability_table(probability_table, CALLER);

  std::vector<types::real_n> result(static_cast<std::size_t>(probability_table.rows()),
                                    types::real_n{0.0});

  for (Eigen::Index row = 0; row < probability_table.rows(); ++row) {
    for (Eigen::Index col = 0; col < probability_table.cols(); ++col) {
      result[static_cast<std::size_t>(row)] += probability_table(row, col);
    }
  }

  return result;
}

std::vector<types::real_n> marginal_y(const types::r_mat& probability_table) {
  constexpr const char* CALLER = "statistics::marginal_y()";
  validate_joint_probability_table(probability_table, CALLER);

  std::vector<types::real_n> result(static_cast<std::size_t>(probability_table.cols()),
                                    types::real_n{0.0});

  for (Eigen::Index row = 0; row < probability_table.rows(); ++row) {
    for (Eigen::Index col = 0; col < probability_table.cols(); ++col) {
      result[static_cast<std::size_t>(col)] += probability_table(row, col);
    }
  }

  return result;
}

types::real_n mean(const std::vector<types::real_n>& probabilities,
                   const std::vector<types::real_n>& values) {
  constexpr const char* CALLER = "statistics::mean()";
  validate_probability_values_pair(probabilities, values, CALLER);

  types::real_n result = types::real_n{0.0};

  for (std::size_t index = 0; index < probabilities.size(); ++index) {
    result += probabilities[index] * values[index];
  }

  return result;
}

types::real_n variance(const std::vector<types::real_n>& probabilities,
                       const std::vector<types::real_n>& values) {
  constexpr const char* CALLER = "statistics::variance()";
  validate_probability_values_pair(probabilities, values, CALLER);

  const types::real_n expected = mean(probabilities, values);

  types::real_n result = types::real_n{0.0};

  for (std::size_t index = 0; index < probabilities.size(); ++index) {
    const types::real_n delta = values[index] - expected;
    result += probabilities[index] * delta * delta;
  }

  return result;
}

types::real_n covariance(const types::r_mat& joint_probabilities,
                         const std::vector<types::real_n>& x_values,
                         const std::vector<types::real_n>& y_values) {
  constexpr const char* CALLER = "statistics::covariance()";
  validate_joint_values(joint_probabilities, x_values, y_values, CALLER);

  const std::vector<types::real_n> x_marginal = marginal_x(joint_probabilities);
  const std::vector<types::real_n> y_marginal = marginal_y(joint_probabilities);

  const types::real_n x_mean = mean(x_marginal, x_values);
  const types::real_n y_mean = mean(y_marginal, y_values);

  types::real_n expected_xy = types::real_n{0.0};

  for (Eigen::Index row = 0; row < joint_probabilities.rows(); ++row) {
    for (Eigen::Index col = 0; col < joint_probabilities.cols(); ++col) {
      expected_xy += joint_probabilities(row, col) *
                     x_values[static_cast<std::size_t>(row)] *
                     y_values[static_cast<std::size_t>(col)];
    }
  }

  return expected_xy - x_mean * y_mean;
}

types::real_n correlation(const types::r_mat& joint_probabilities,
                          const std::vector<types::real_n>& x_values,
                          const std::vector<types::real_n>& y_values) {
  constexpr const char* CALLER = "statistics::correlation()";
  validate_joint_values(joint_probabilities, x_values, y_values, CALLER);

  const std::vector<types::real_n> x_marginal = marginal_x(joint_probabilities);
  const std::vector<types::real_n> y_marginal = marginal_y(joint_probabilities);

  const types::real_n x_variance = variance(x_marginal, x_values);
  const types::real_n y_variance = variance(y_marginal, y_values);

  if (x_variance <= PROBABILITY_TOLERANCE || y_variance <= PROBABILITY_TOLERANCE) {
    throw exception::InvalidArgument(CALLER, "correlation requires non-zero variance");
  }

  return covariance(joint_probabilities, x_values, y_values) /
         std::sqrt(x_variance * y_variance);
}

}  // namespace quira::statistics
