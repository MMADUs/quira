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

#include "quira/gates/two.hpp"

#include "quira/constants.hpp"
#include "quira/exception.hpp"

#include <sstream>
#include <unordered_set>
#include <vector>

namespace quira {

namespace {

void validate_distinct_qubits(const char* gate_name,
                              const std::vector<types::qubit>& qubits) {
  std::unordered_set<types::qubit> seen;

  for (types::qubit qubit : qubits) {
    if (!seen.insert(qubit).second) {
      std::ostringstream oss;
      oss << gate_name << " received duplicate qubit index " << qubit;
      throw exception::InvalidArgument(gate_name, oss.str());
    }
  }
}

}  // namespace

CNOT::CNOT(types::qubit control, types::qubit target)
    : TwoQubit<CNOT>(control, target) {
  validate_distinct_qubits("CNOT", {control, target});
}

std::string CNOT::name() const {
  std::ostringstream oss;
  oss << "CNOT(control=" << control() << ", target=" << target() << ")";
  return oss.str();
}

types::c_mat CNOT::unitary() const {
  types::c_mat matrix = types::c_mat::Zero(4, 4);

  matrix(0, 0) = constants::ONE;
  matrix(3, 1) = constants::ONE;
  matrix(2, 2) = constants::ONE;
  matrix(1, 3) = constants::ONE;

  return matrix;
}

types::qubit CNOT::control() const noexcept {
  return first();
}

types::qubit CNOT::target() const noexcept {
  return second();
}

CZ::CZ(types::qubit control, types::qubit target) : TwoQubit<CZ>(control, target) {
  validate_distinct_qubits("CZ", {control, target});
}

std::string CZ::name() const {
  std::ostringstream oss;
  oss << "CZ(control=" << control() << ", target=" << target() << ")";
  return oss.str();
}

types::c_mat CZ::unitary() const {
  types::c_mat matrix = types::c_mat::Identity(4, 4);
  matrix(3, 3) = -constants::ONE;
  return matrix;
}

types::qubit CZ::control() const noexcept {
  return first();
}

types::qubit CZ::target() const noexcept {
  return second();
}

Swap::Swap(types::qubit first, types::qubit second) : TwoQubit<Swap>(first, second) {
  validate_distinct_qubits("SWAP", {first, second});
}

std::string Swap::name() const {
  std::ostringstream oss;
  oss << "SWAP(first=" << first() << ", second=" << second() << ")";
  return oss.str();
}

types::c_mat Swap::unitary() const {
  types::c_mat matrix = types::c_mat::Identity(4, 4);

  matrix(1, 1) = constants::ZERO;
  matrix(2, 2) = constants::ZERO;
  matrix(1, 2) = constants::ONE;
  matrix(2, 1) = constants::ONE;

  return matrix;
}

}  // namespace quira
