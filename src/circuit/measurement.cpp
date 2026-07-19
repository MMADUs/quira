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

#include "quira/circuit.hpp"
#include "quira/exception.hpp"
#include "quira/instruction.hpp"
#include "quira/linalg.hpp"

#include <unordered_set>

namespace quira {

namespace {

void validate_distinct_qubits(const std::vector<types::qubit>& qubits,
                              const char* caller) {
  std::unordered_set<types::qubit> seen;

  for (types::qubit qubit : qubits) {
    if (!seen.insert(qubit).second) {
      throw exception::InvalidArgument(caller, "contains duplicate qubits");
    }
  }
}

}  // namespace

QuantumCircuit& QuantumCircuit::measure(types::qubit qubit, types::clbit clbit) {
  validate_qubit(qubit, "QuantumCircuit::measure()");
  validate_clbit(clbit, "QuantumCircuit::measure()");

  instructions_.emplace_back(MeasureInstruction{.qubit = qubit, .clbit = clbit});

  return *this;
}

QuantumCircuit& QuantumCircuit::measure(const std::vector<types::qubit>& qubits,
                                        const std::vector<types::qubit>& clbits) {
  validate_measurement_lists(qubits, clbits, "QuantumCircuit::measure()");

  for (std::size_t index = 0; index < qubits.size(); ++index) {
    instructions_.emplace_back(
        MeasureInstruction{.qubit = qubits[index], .clbit = clbits[index]});
  }

  return *this;
}

QuantumCircuit& QuantumCircuit::measure_all() {
  if (num_clbits_ < num_qubits_) {
    throw exception::OutOfRange("QuantumCircuit::measure_all()",
                                "requires at least as many classical bits as qubits");
  }

  for (std::size_t i = 0; i < num_qubits_; ++i) {
    instructions_.emplace_back(MeasureInstruction{
        .qubit = static_cast<types::qubit>(i), .clbit = static_cast<types::clbit>(i)});
  }

  return *this;
}

QuantumCircuit& QuantumCircuit::measure_basis(const types::c_mat& basis,
                                              const std::vector<types::qubit>& qubits,
                                              const std::vector<types::clbit>& clbits) {
  validate_measurement_lists(qubits, clbits, "QuantumCircuit::measure_basis()");

  const std::size_t dimension = linalg::basis_size(qubits.size());

  if (basis.rows() != static_cast<Eigen::Index>(dimension) ||
      basis.cols() != static_cast<Eigen::Index>(dimension)) {
    std::ostringstream oss;
    oss << "measurement basis dimension does not match qubits, got shape ("
        << basis.rows() << ", " << basis.cols() << ") for expected basis size "
        << dimension;

    throw exception::DimensionMismatch("QuantumCircuit::measure_basis()", oss.str());
  }

  if (!linalg::is_unitary(basis)) {
    throw exception::InvalidArgument("QuantumCircuit::measure_basis()",
                                     "requires a unitary basis");
  }

  instructions_.emplace_back(
      MeasureBasisInstruction{.basis = basis, .qubits = qubits, .clbits = clbits});

  return *this;
}

QuantumCircuit& QuantumCircuit::post_select(types::qubit qubit, bool expected) {
  validate_qubit(qubit, "QuantumCircuit::post_select()");
  instructions_.emplace_back(
      PostSelectInstruction{.qubits = {qubit}, .expected = {expected}});

  return *this;
}

QuantumCircuit& QuantumCircuit::post_select(const std::vector<types::qubit>& qubits,
                                            const std::vector<bool>& expected) {
  if (qubits.empty()) {
    throw exception::InvalidArgument("QuantumCircuit::post_select()",
                                     "requires qubits, found empty");
  }

  if (qubits.size() != expected.size()) {
    throw exception::InvalidArgument("QuantumCircuit::post_select()",
                                     "requires matching vector size");
  }

  validate_qubits(qubits, "QuantumCircuit::post_select()");
  validate_distinct_qubits(qubits, "QuantumCircuit::post_select()");

  instructions_.emplace_back(
      PostSelectInstruction{.qubits = qubits, .expected = expected});

  return *this;
}

}  // namespace quira
