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
#include "quira/linalg.hpp"

#include <sstream>
#include <unordered_set>

namespace quira {

namespace {

template<typename T>
void validate_distinct(const std::vector<T>& values, const exception::ErrorCode code,
                       const char* message, const char* caller) {
  std::unordered_set<T> seen;

  for (T value : values) {
    if (!seen.insert(value).second) {
      if (code == exception::ErrorCode::InvalidQubit) {
        throw exception::InvalidQubit(caller, message);
      } else if (code == exception::ErrorCode::InvalidClbit) {
        throw exception::InvalidClbit(caller, message);
      } else {
      }
    }
  }
}

}  // namespace

void QuantumCircuit::validate_qubit(types::qubit qubit, const char* caller) const {
  if (qubit >= num_qubits_) {
    std::ostringstream oss;
    oss << "qubit index " << qubit << " is out of range for circuit with "
        << num_qubits_ << " qubits";

    throw exception::OutOfRange(caller, oss.str());
  }
}

void QuantumCircuit::validate_clbit(types::clbit clbit, const char* caller) const {
  if (clbit >= num_clbits_) {
    std::ostringstream oss;
    oss << "classical bit index " << clbit << " is out of range for circuit with "
        << num_clbits_ << " classical bits";

    throw exception::OutOfRange(caller, oss.str());
  }
}

void QuantumCircuit::validate_qubits(const std::vector<types::qubit>& qubits,
                                     const char* caller) const {
  for (types::qubit qubit : qubits) {
    validate_qubit(qubit, caller);
  }
}

void QuantumCircuit::validate_clbits(const std::vector<types::clbit>& clbits,
                                     const char* caller) const {
  for (types::clbit clbit : clbits) {
    validate_clbit(clbit, caller);
  }
}

void QuantumCircuit::validate_gate(const QuantumGate& gate, const char* caller) const {
  const auto targets = gate.targets();

  if (targets.empty()) {
    throw exception::InvalidArgument(caller,
                                     "quantum gate must target at least one qubit");
  }

  validate_qubits(targets, caller);
  validate_distinct(targets, exception::ErrorCode::InvalidQubit, caller,
                    "quantum gate contains duplicate target qubits");
}

void QuantumCircuit::validate_measurement_lists(const std::vector<types::qubit>& qubits,
                                                const std::vector<types::clbit>& clbits,
                                                const char* caller) const {
  if (qubits.empty()) {
    throw exception::InvalidArgument(caller, "measurement requires at least one qubit");
  }

  if (qubits.size() != clbits.size()) {
    throw exception::InvalidArgument(
        caller, "measurement requires matching qubit and classical bit counts");
  }

  validate_qubits(qubits, caller);
  validate_clbits(clbits, caller);

  validate_distinct(qubits, exception::ErrorCode::InvalidClbit, caller,
                    "Measurement contains duplicate qubits");
  validate_distinct(clbits, exception::ErrorCode::InvalidQubit, caller,
                    "Measurement contains duplicate classical bits");
}

bool QuantumCircuit::is_valid() const {
  try {
    for (const Instruction& instruction : instructions_) {
      std::visit(
          types::Overloaded{
              [&](const GateInstruction& item) {
                validate_gate(item.gate(), "QuantumCircuit::is_valid()");

                if (item.condition()) {
                  validate_clbit(item.condition()->clbit, "QuantumCircuit::is_valid()");
                }
              },
              [&](const MeasureInstruction& item) {
                validate_qubit(item.qubit, "QuantumCircuit::is_valid()");
                validate_clbit(item.clbit, "QuantumCircuit::is_valid()");
              },
              [&](const MeasureBasisInstruction& item) {
                validate_measurement_lists(item.qubits, item.clbits,
                                           "QuantumCircuit::is_valid()");

                const std::size_t dimension = linalg::basis_size(item.qubits.size());

                if (item.basis.rows() != static_cast<Eigen::Index>(dimension) ||
                    item.basis.cols() != static_cast<Eigen::Index>(dimension)) {
                  std::ostringstream oss;
                  oss << "measurement basis dimension does not match qubits, got shape "
                         "("
                      << item.basis.rows() << ", " << item.basis.cols()
                      << ") for expected basis size " << dimension;

                  throw exception::DimensionMismatch("QuantumCircuit::is_valid()",
                                                     oss.str());
                }

                if (!linalg::is_unitary(item.basis)) {
                  throw exception::InvalidArgument(
                      "QuantumCircuit::is_valid()",
                      "MeasureBasisInstruction requires a unitary basis");
                }
              },
              [&](const PostSelectInstruction& item) {
                if (item.qubits.empty()) {
                  throw exception::InvalidArgument(
                      "QuantumCircuit::is_valid()",
                      "PostSelectInstruction requires at least one qubit");
                }

                if (item.qubits.size() != item.expected.size()) {
                  throw exception::InvalidArgument(
                      "QuantumCircuit::is_valid()",
                      "PostSelectInstruction qubits and expected sizes must match");
                }

                validate_qubits(item.qubits, "QuantumCircuit::is_valid()");
                validate_distinct(item.qubits, exception::ErrorCode::InvalidQubit,
                                  "QuantumCircuit::is_valid()",
                                  "PostSelectInstruction contains duplicate qubits");
              },
              [&](const ResetInstruction& item) {
                validate_qubit(item.qubit, "QuantumCircuit::is_valid()");
              },
              [](const BarrierInstruction&) {
                // Barriers carry no indices or payload to validate.
              },
          },
          instruction);
    }

    return true;
  } catch (...) {
    return false;
  }
}

}  // namespace quira
