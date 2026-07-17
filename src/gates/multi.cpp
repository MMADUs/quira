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

#include "quira/gates/multi.hpp"

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

Toff::Toff(types::qubit first_control, types::qubit second_control, types::qubit target)
    : MultiQubit<Toff>({first_control, second_control, target}) {
  validate_distinct_qubits("Toff", target_list());
}

std::string Toff::name() const {
  std::ostringstream oss;
  oss << "Toff(first_control=" << first_control()
      << ", second_control=" << second_control() << ", target=" << target() << ")";
  return oss.str();
}

types::c_mat Toff::unitary() const {
  types::c_mat matrix = types::c_mat::Identity(8, 8);

  matrix(6, 6) = constants::ZERO;
  matrix(7, 7) = constants::ZERO;
  matrix(6, 7) = constants::ONE;
  matrix(7, 6) = constants::ONE;

  return matrix;
}

types::qubit Toff::first_control() const noexcept {
  return target_list()[0];
}

types::qubit Toff::second_control() const noexcept {
  return target_list()[1];
}

types::qubit Toff::target() const noexcept {
  return target_list()[2];
}

Fred::Fred(types::qubit control, types::qubit first_target, types::qubit second_target)
    : MultiQubit<Fred>({control, first_target, second_target}) {
  validate_distinct_qubits("Fred", target_list());
}

std::string Fred::name() const {
  std::ostringstream oss;
  oss << "Fred(control=" << control() << ", first_target=" << first_target()
      << ", second_target=" << second_target() << ")";
  return oss.str();
}

types::c_mat Fred::unitary() const {
  types::c_mat matrix = types::c_mat::Identity(8, 8);

  matrix(5, 5) = constants::ZERO;
  matrix(6, 6) = constants::ZERO;
  matrix(5, 6) = constants::ONE;
  matrix(6, 5) = constants::ONE;

  return matrix;
}

types::qubit Fred::control() const noexcept {
  return target_list()[0];
}

types::qubit Fred::first_target() const noexcept {
  return target_list()[1];
}

types::qubit Fred::second_target() const noexcept {
  return target_list()[2];
}

}  // namespace quira
