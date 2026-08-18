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

#include "quira/register.hpp"

#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace quira {

/**
 * @brief Aggregated output from running a circuit.
 *
 * SimulationOutput stores one ClassicalRegister for each accepted shot. If a
 * circuit contains post-selection, some requested shots may be rejected and are
 * tracked separately from accepted shots.
 *
 * Counts and probabilities are computed from accepted shots only.
 */
class SimulationOutput {
public:
  /**
   * @brief Creates simulation output from accepted shot registers.
   *
   * @param shots Classical registers produced by accepted shots.
   * @param requested_shots Number of shots requested by the caller.
   * @param rejected_shots Number of shots rejected by post-selection.
   */
  SimulationOutput(std::vector<ClassicalRegister> shots, std::size_t requested_shots,
                   std::size_t rejected_shots = 0);

  /**
   * @brief Returns the number of requested shots.
   *
   * @return Requested shot count, including accepted and rejected shots.
   */
  [[nodiscard]] std::size_t shots() const noexcept;

  /**
   * @brief Returns the number of accepted shots.
   *
   * @return Number of stored shot registers.
   */
  [[nodiscard]] std::size_t accepted_shots() const noexcept;

  /**
   * @brief Returns the number of rejected shots.
   *
   * @return Number of shots rejected by post-selection.
   */
  [[nodiscard]] std::size_t rejected_shots() const noexcept;

  /**
   * @brief Returns all accepted per-shot classical registers.
   *
   * @return Read-only accepted shot register storage.
   */
  [[nodiscard]] const std::vector<ClassicalRegister>& shot_registers() const noexcept;

  /**
   * @brief Returns one accepted shot's classical register.
   *
   * @param index Accepted shot index.
   * @return Read-only classical register for the selected accepted shot.
   *
   * @throws std::out_of_range If index is outside accepted shots.
   */
  [[nodiscard]] const ClassicalRegister& shot(std::size_t index) const;

  /**
   * @brief Counts accepted measurement outcomes.
   *
   * @return Map from measured bit string to occurrence count.
   *
   * @note Shots with no measured classical bits are ignored.
   */
  [[nodiscard]] std::unordered_map<std::string, std::size_t> get_counts() const;

  /**
   * @brief Computes empirical probabilities from accepted shots.
   *
   * @return Map from measured bit string to probability.
   *
   * @note Probabilities use accepted_shots() as the denominator, not shots().
   */
  [[nodiscard]] std::unordered_map<std::string, double> get_probabilities() const;

  /**
   * @brief Returns the most frequent accepted measurement outcome.
   *
   * @return Empty optional when there are no measured outcomes.
   */
  [[nodiscard]] std::optional<std::pair<std::string, std::size_t>> get_most_frequent()
      const;

private:
  std::vector<ClassicalRegister> shot_registers_;
  std::size_t requested_shots_{};
  std::size_t rejected_shots_{};
};

}  // namespace quira
