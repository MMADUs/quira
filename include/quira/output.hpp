// Copyright 2026 Quira contributors
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

#include "quira/register.hpp"
#include "quira/types.hpp"

#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace quira {

/**
 * @brief Aggregated result from running a circuit for one or more shots.
 *
 * RunResult stores the classical register produced by each shot and provides
 * Qiskit-like helpers for counts, probabilities, and the most frequent
 * measurement outcome.
 */
class SimulationOutput {
public:
  /**
   * @brief Creates a run result from per-shot classical registers.
   *
   * @param shots Classical registers produced by individual shots.
   * @param shot_count Number of requested shots.
   */
  SimulationOutput(std::vector<ClassicalRegister> shots, Index shot_count);

  /**
   * @brief Returns the number of requested shots.
   */
  Index shots() const noexcept;

  /**
   * @brief Returns all per-shot classical registers.
   *
   * @return Read-only shot register storage.
   */
  const std::vector<ClassicalRegister>& shot_registers() const noexcept;

  /**
   * @brief Returns one shot's classical register.
   *
   * @param index Shot index.
   * @return Read-only classical register for the selected shot.
   *
   * @throws std::out_of_range If index is outside the stored shots.
   */
  const ClassicalRegister& shot(Index index) const;

  /**
   * @brief Counts measurement outcomes across all shots.
   *
   * @return Map from measured bit string to occurrence count.
   *
   * @note Shots with no measured classical bits are ignored.
   */
  std::unordered_map<std::string, Index> get_counts() const;

  /**
   * @brief Computes measurement outcome probabilities.
   *
   * @return Map from measured bit string to empirical probability.
   *
   * @note Probabilities are computed as count divided by the requested shot
   * count.
   */
  std::unordered_map<std::string, double> get_probabilities() const;

  /**
   * @brief Returns the most frequent measurement outcome.
   *
   * @return Empty optional when there are no measured outcomes.
   */
  std::optional<std::pair<std::string, Index>> get_most_frequent() const;

private:
  std::vector<ClassicalRegister> shot_registers_;
  Index shots_{};
};

}  // namespace quira
