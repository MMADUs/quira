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

#include "quira/gate.hpp"

#include <string>
#include <vector>

namespace quira {

/**
 * @brief Generic explicit unitary gate.
 *
 * Unitary stores a unitary matrix and an ordered target-qubit list directly.
 * It is useful for user-provided gates, circuit composition, and circuit
 * adjoints where preserving the original concrete gate type is not required.
 *
 * @note The target order follows the same convention as QuantumGate::targets():
 * targets[0] maps to local bit 0 in the gate matrix.
 */
class Unitary final : public CloneableGate<Unitary> {
public:
  /**
   * @brief Creates an explicit unitary gate.
   *
   * @param name Human-readable gate name.
   * @param unitary Unitary matrix with dimension 2^k by 2^k, where k is the
   * number of targets.
   * @param targets Qubits acted on by the unitary.
   *
   * @throws std::invalid_argument If targets are empty, contain duplicates, the
   * unitary dimension does not match target count, or the matrix is not unitary.
   */
  Unitary(std::string name, types::c_mat unitary, std::vector<types::qubit> targets);

  /**
   * @brief Returns the human-readable gate name.
   */
  [[nodiscard]] std::string name() const override;

  /**
   * @brief Returns the explicit unitary matrix.
   */
  [[nodiscard]] types::c_mat unitary() const override;

  /**
   * @brief Returns the qubits acted on by the gate.
   *
   * @return Target qubits in the order expected by the unitary matrix.
   */
  [[nodiscard]] std::vector<types::qubit> targets() const override;

private:
  std::string name_;
  types::c_mat unitary_;
  std::vector<types::qubit> targets_;
};

}  // namespace quira
