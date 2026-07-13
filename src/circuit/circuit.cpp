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

#include "quira/circuit.hpp"

#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace quira {

GateInstruction::GateInstruction(std::unique_ptr<QuantumGate> gate)
    : gate_(std::move(gate)) {
  if (!gate_) {
    throw std::invalid_argument("GateInstruction requires a non-null gate");
  }
}

GateInstruction::GateInstruction(std::unique_ptr<QuantumGate> gate,
                                 ClassicalCondition condition)
    : gate_(std::move(gate)), condition_(condition) {
  if (!gate_) {
    throw std::invalid_argument("GateInstruction requires a non-null gate");
  }
}

GateInstruction::GateInstruction(const GateInstruction& other)
    : gate_(other.gate_ ? other.gate_->clone() : nullptr),
      condition_(other.condition_) {
}

GateInstruction& GateInstruction::operator=(const GateInstruction& other) {
  if (this == &other) {
    return *this;
  }

  gate_ = other.gate_ ? other.gate_->clone() : nullptr;
  condition_ = other.condition_;
  return *this;
}

const QuantumGate& GateInstruction::gate() const noexcept {
  return *gate_;
}

const std::optional<ClassicalCondition>& GateInstruction::condition() const noexcept {
  return condition_;
}

QuantumCircuit::QuantumCircuit(Index num_qubits, Index num_clbits)
    : num_qubits_(num_qubits), num_clbits_(num_clbits) {
}

Index QuantumCircuit::num_qubits() const noexcept {
  return num_qubits_;
}

Index QuantumCircuit::num_clbits() const noexcept {
  return num_clbits_;
}

const std::vector<Instruction>& QuantumCircuit::instructions() const noexcept {
  return instructions_;
}

QuantumCircuit& QuantumCircuit::add(const QuantumGate& gate) {
  validate_gate(gate);
  instructions_.emplace_back(GateInstruction{gate.clone()});
  return *this;
}

QuantumCircuit& QuantumCircuit::add(std::unique_ptr<QuantumGate> gate) {
  if (!gate) {
    throw std::invalid_argument("QuantumCircuit::add received a null gate");
  }

  validate_gate(*gate);
  instructions_.emplace_back(GateInstruction{std::move(gate)});
  return *this;
}

QuantumCircuit& QuantumCircuit::conditional_add(Clbit clbit, bool value,
                                                const QuantumGate& gate) {
  validate_clbit(clbit);
  validate_gate(gate);

  instructions_.emplace_back(
      GateInstruction{gate.clone(), ClassicalCondition{clbit, value}});
  return *this;
}

QuantumCircuit& QuantumCircuit::conditional_add(Clbit clbit, bool value,
                                                std::unique_ptr<QuantumGate> gate) {
  validate_clbit(clbit);

  if (!gate) {
    throw std::invalid_argument("QuantumCircuit::conditional_add received a null gate");
  }

  validate_gate(*gate);

  instructions_.emplace_back(
      GateInstruction{std::move(gate), ClassicalCondition{clbit, value}});
  return *this;
}

QuantumCircuit& QuantumCircuit::measure(Qubit qubit, Clbit clbit) {
  validate_qubit(qubit);
  validate_clbit(clbit);

  instructions_.emplace_back(MeasureInstruction{qubit, clbit});
  return *this;
}

QuantumCircuit& QuantumCircuit::measure_all() {
  if (num_clbits_ < num_qubits_) {
    throw std::out_of_range(
        "QuantumCircuit::measure_all requires at least as many classical "
        "bits as qubits");
  }

  for (Qubit qubit = 0; qubit < num_qubits_; ++qubit) {
    instructions_.emplace_back(MeasureInstruction{qubit, qubit});
  }

  return *this;
}

QuantumCircuit& QuantumCircuit::reset(Qubit qubit) {
  validate_qubit(qubit);
  instructions_.emplace_back(ResetInstruction{qubit});
  return *this;
}

QuantumCircuit& QuantumCircuit::barrier() {
  instructions_.emplace_back(BarrierInstruction{});
  return *this;
}

void QuantumCircuit::validate_qubit(Qubit qubit) const {
  if (qubit >= num_qubits_) {
    std::ostringstream oss;
    oss << "Qubit index " << qubit << " is out of range for circuit with "
        << num_qubits_ << " qubits";
    throw std::out_of_range(oss.str());
  }
}

void QuantumCircuit::validate_clbit(Clbit clbit) const {
  if (clbit >= num_clbits_) {
    std::ostringstream oss;
    oss << "Classical bit index " << clbit << " is out of range for circuit with "
        << num_clbits_ << " classical bits";
    throw std::out_of_range(oss.str());
  }
}

void QuantumCircuit::validate_gate(const QuantumGate& gate) const {
  const auto targets = gate.targets();

  if (targets.empty()) {
    throw std::invalid_argument("Quantum gate must target at least one qubit");
  }

  std::unordered_set<Qubit> seen_targets;
  for (Qubit target : targets) {
    validate_qubit(target);

    const auto inserted = seen_targets.insert(target).second;
    if (!inserted) {
      std::ostringstream oss;
      oss << "Quantum gate '" << gate.name() << "' contains duplicate target qubit "
          << target;
      throw std::invalid_argument(oss.str());
    }
  }
}

}  // namespace quira
