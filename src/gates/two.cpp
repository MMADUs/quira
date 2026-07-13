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

#include "quira/gates/two.hpp"

#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace quira {
namespace {

void validate_distinct_qubits(const char* gate_name, const std::vector<Qubit>& qubits) {
  std::unordered_set<Qubit> seen;
  for (Qubit qubit : qubits) {
    if (!seen.insert(qubit).second) {
      std::ostringstream oss;
      oss << gate_name << " received duplicate qubit index " << qubit;
      throw std::invalid_argument(oss.str());
    }
  }
}

}  // namespace

/**
 * Controlled-NOT.
 *
 * Flips the target qubit when the control qubit is in the |1> state.
 */
ControlledNot::ControlledNot(Qubit control, Qubit target)
    : TwoQubit<ControlledNot>(control, target) {
  validate_distinct_qubits("CNOT", {first_, second_});
}

std::string ControlledNot::name() const {
  std::ostringstream oss;
  oss << "CNOT(control=" << control() << ", target=" << target() << ")";
  return oss.str();
}

Matrix ControlledNot::unitary() const {
  Matrix matrix = Matrix::Zero(4, 4);

  matrix(0, 0) = 1.0;
  matrix(3, 1) = 1.0;
  matrix(2, 2) = 1.0;
  matrix(1, 3) = 1.0;

  return matrix;
}

Qubit ControlledNot::control() const noexcept {
  return first_;
}

Qubit ControlledNot::target() const noexcept {
  return second_;
}

/**
 * Controlled-Z.
 *
 * Applies a phase flip when both qubits are in the |1> state.
 */
ControlledZ::ControlledZ(Qubit control, Qubit target)
    : TwoQubit<ControlledZ>(control, target) {
  validate_distinct_qubits("CZ", {first_, second_});
}

std::string ControlledZ::name() const {
  std::ostringstream oss;
  oss << "CZ(control=" << control() << ", target=" << target() << ")";
  return oss.str();
}

Matrix ControlledZ::unitary() const {
  Matrix matrix = Matrix::Identity(4, 4);
  matrix(3, 3) = -1.0;
  return matrix;
}

Qubit ControlledZ::control() const noexcept {
  return first_;
}

Qubit ControlledZ::target() const noexcept {
  return second_;
}

/**
 * Swap.
 *
 * Exchanges the quantum states of two qubits.
 */
Swap::Swap(Qubit first, Qubit second) : TwoQubit<Swap>(first, second) {
  validate_distinct_qubits("SWAP", {first_, second_});
}

std::string Swap::name() const {
  std::ostringstream oss;
  oss << "SWAP(first=" << first() << ", second=" << second() << ")";
  return oss.str();
}

Matrix Swap::unitary() const {
  Matrix matrix = Matrix::Identity(4, 4);
  matrix(1, 1) = 0.0;
  matrix(2, 2) = 0.0;
  matrix(1, 2) = 1.0;
  matrix(2, 1) = 1.0;
  return matrix;
}

}  // namespace quira
