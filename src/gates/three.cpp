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

#include "quira/gates/three.hpp"

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
 * Toffoli.
 *
 * Flips the target qubit when both control qubits are in the |1> state.
 */
Toffoli::Toffoli(Qubit first_control, Qubit second_control, Qubit target)
    : ThreeQubit<Toffoli>(first_control, second_control, target) {
  validate_distinct_qubits("TOFFOLI", {first_, second_, third_});
}

std::string Toffoli::name() const {
  std::ostringstream oss;
  oss << "TOFFOLI(first_control=" << first_control()
      << ", second_control=" << second_control() << ", target=" << target() << ")";
  return oss.str();
}

Matrix Toffoli::unitary() const {
  Matrix matrix = Matrix::Identity(8, 8);

  matrix(3, 3) = 0.0;
  matrix(7, 7) = 0.0;
  matrix(7, 3) = 1.0;
  matrix(3, 7) = 1.0;

  return matrix;
}

Qubit Toffoli::first_control() const noexcept {
  return first_;
}

Qubit Toffoli::second_control() const noexcept {
  return second_;
}

Qubit Toffoli::target() const noexcept {
  return third_;
}

/**
 * Fredkin.
 *
 * Swaps the two target qubits when the control qubit is in the |1> state.
 */
Fredkin::Fredkin(Qubit control, Qubit first_target, Qubit second_target)
    : ThreeQubit<Fredkin>(control, first_target, second_target) {
  validate_distinct_qubits("FREDKIN", {first_, second_, third_});
}

std::string Fredkin::name() const {
  std::ostringstream oss;
  oss << "FREDKIN(control=" << control() << ", first_target=" << first_target()
      << ", second_target=" << second_target() << ")";
  return oss.str();
}

Matrix Fredkin::unitary() const {
  Matrix matrix = Matrix::Identity(8, 8);

  // Qiskit-style little endian:
  // local index = control + 2 * first_target + 4 * second_target.
  // When control is 1, swap target states 01 and 10: 3 <-> 5.
  matrix(3, 3) = 0.0;
  matrix(5, 5) = 0.0;
  matrix(5, 3) = 1.0;
  matrix(3, 5) = 1.0;

  return matrix;
}

Qubit Fredkin::control() const noexcept {
  return first_;
}

Qubit Fredkin::first_target() const noexcept {
  return second_;
}

Qubit Fredkin::second_target() const noexcept {
  return third_;
}

}  // namespace quira
