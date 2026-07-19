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

namespace quira {

QuantumCircuit::QuantumCircuit(std::size_t num_qubits, std::size_t num_clbits)
    : num_qubits_(num_qubits), num_clbits_(num_clbits) {
}

QuantumCircuit::QuantumCircuit(const QuantumRegister& qreg,
                               const ClassicalRegister& creg)
    : QuantumCircuit(qreg.size(), creg.size()) {
}

std::size_t QuantumCircuit::num_qubits() const noexcept {
  return num_qubits_;
}

std::size_t QuantumCircuit::num_clbits() const noexcept {
  return num_clbits_;
}

const std::vector<Instruction>& QuantumCircuit::instructions() const noexcept {
  return instructions_;
}

QuantumCircuit& QuantumCircuit::add(const QuantumGate& gate) {
  validate_gate(gate, "QuantumCircuit::add()");
  instructions_.emplace_back(GateInstruction{gate.clone()});

  return *this;
}

QuantumCircuit& QuantumCircuit::add(std::unique_ptr<QuantumGate> gate) {
  if (!gate) {
    throw exception::InvalidArgument("QuantumCircuit::add()", "received a null gate");
  }

  validate_gate(*gate, "QuantumCircuit::add()");
  instructions_.emplace_back(GateInstruction{std::move(gate)});

  return *this;
}

QuantumCircuit& QuantumCircuit::conditional_add(types::clbit clbit, bool value,
                                                const QuantumGate& gate) {
  validate_clbit(clbit, "QuantumCircuit::conditional_add()");
  validate_gate(gate, "QuantumCircuit::conditional_add()");
  instructions_.emplace_back(GateInstruction{
      gate.clone(), ClassicalCondition{.clbit = clbit, .value = value}});

  return *this;
}

QuantumCircuit& QuantumCircuit::conditional_add(types::clbit clbit, bool value,
                                                std::unique_ptr<QuantumGate> gate) {
  validate_clbit(clbit, "QuantumCircuit::conditional_add()");

  if (!gate) {
    throw exception::InvalidArgument("QuantumCircuit::conditional_add()",
                                     "received a null gate");
  }

  validate_gate(*gate, "QuantumCircuit::conditional_add()");
  instructions_.emplace_back(GateInstruction{
      std::move(gate), ClassicalCondition{.clbit = clbit, .value = value}});

  return *this;
}

QuantumCircuit& QuantumCircuit::reset(types::qubit qubit) {
  validate_qubit(qubit, "QuantumCircuit::reset()");
  instructions_.emplace_back(ResetInstruction{qubit});

  return *this;
}

QuantumCircuit& QuantumCircuit::barrier() {
  instructions_.emplace_back(BarrierInstruction{});

  return *this;
}

}  // namespace quira
