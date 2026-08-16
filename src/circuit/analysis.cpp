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
#include "quira/unitary.hpp"

#include <algorithm>
#include <ranges>
#include <unordered_set>
#include <variant>

namespace quira {

namespace {

std::vector<types::qubit> touched_qubits(const Instruction& instruction) {
  return std::visit(
      types::Overloaded{
          [](const GateInstruction& item) -> std::vector<types::qubit> {
            return item.gate().targets();
          },
          [](const MeasureInstruction& item) -> std::vector<types::qubit> {
            return {item.qubit};
          },
          [](const MeasureBasisInstruction& item) -> std::vector<types::qubit> {
            return item.qubits;
          },
          [](const PostSelectInstruction& item) -> std::vector<types::qubit> {
            return item.qubits;
          },
          [](const ResetInstruction& item) -> std::vector<types::qubit> {
            return {item.qubit};
          },
          [](const BarrierInstruction&) -> std::vector<types::qubit> { return {}; },
      },
      instruction);
}

}  // namespace

QuantumCircuit QuantumCircuit::adjoint() const {
  QuantumCircuit result{num_qubits_, num_clbits_};

  for (const Instruction& instruction : std::views::reverse(instructions_)) {
    std::visit(types::Overloaded{
                   [&](const GateInstruction& item) {
                     if (item.condition()) {
                       throw exception::InvalidArgument(
                           "QuantumCircuit::adjoint()",
                           "cannot adjoint a circuit with conditional gates");
                     }

                     const QuantumGate& gate = item.gate();

                     auto adjoint_gate = std::make_unique<Unitary>(
                         gate.name() + "^dagger", linalg::adjoint(gate.unitary()),
                         gate.targets());

                     result.add(std::move(adjoint_gate));
                   },
                   [&](const BarrierInstruction&) { result.barrier(); },
                   [](const MeasureInstruction&) {
                     throw exception::InvalidArgument(
                         "QuantumCircuit::adjoint()",
                         "cannot adjoint a circuit with non-unitary instructions");
                   },
                   [](const MeasureBasisInstruction&) {
                     throw exception::InvalidArgument(
                         "QuantumCircuit::adjoint()",
                         "cannot adjoint a circuit with non-unitary instructions");
                   },
                   [](const PostSelectInstruction&) {
                     throw exception::InvalidArgument(
                         "QuantumCircuit::adjoint()",
                         "cannot adjoint a circuit with non-unitary instructions");
                   },
                   [](const ResetInstruction&) {
                     throw exception::InvalidArgument(
                         "QuantumCircuit::adjoint()",
                         "cannot adjoint a circuit with non-unitary instructions");
                   },
               },
               instruction);
  }

  return result;
}

QuantumCircuit QuantumCircuit::inverse() const {
  // inverse is just an alias of adjoint
  return adjoint();
}

std::size_t QuantumCircuit::depth() const {
  std::vector<std::size_t> layers(num_qubits_, 0);

  std::size_t result = 0;

  for (const Instruction& instruction : instructions_) {
    if (std::holds_alternative<BarrierInstruction>(instruction)) {
      ++result;
      std::ranges::fill(layers, result);
      continue;
    }

    const auto qubits = touched_qubits(instruction);

    if (qubits.empty()) {
      continue;
    }

    std::size_t layer = 0;

    for (types::qubit qubit : qubits) {
      layer = std::max(layer, layers[qubit]);
    }

    ++layer;

    for (types::qubit qubit : qubits) {
      layers[qubit] = layer;
    }

    result = std::max(result, layer);
  }

  return result;
}

std::size_t QuantumCircuit::gate_count() const {
  std::size_t count = 0;

  for (const Instruction& instruction : instructions_) {
    if (std::holds_alternative<GateInstruction>(instruction)) {
      ++count;
    }
  }

  return count;
}

std::unordered_map<std::string, std::size_t> QuantumCircuit::gate_counts() const {
  std::unordered_map<std::string, std::size_t> counts;

  for (const Instruction& instruction : instructions_) {
    if (const auto* gate = std::get_if<GateInstruction>(&instruction)) {
      ++counts[gate->gate().name()];
    }
  }

  return counts;
}

std::vector<types::qubit> QuantumCircuit::measured_qubits() const {
  std::unordered_set<types::qubit> measured;

  for (const Instruction& instruction : instructions_) {
    if (const auto* measure = std::get_if<MeasureInstruction>(&instruction)) {
      measured.insert(measure->qubit);
    } else if (const auto* basis = std::get_if<MeasureBasisInstruction>(&instruction)) {
      measured.insert(basis->qubits.begin(), basis->qubits.end());
    } else if (const auto* post = std::get_if<PostSelectInstruction>(&instruction)) {
      measured.insert(post->qubits.begin(), post->qubits.end());
    }
  }

  std::vector<types::qubit> result(measured.begin(), measured.end());
  std::ranges::sort(result);

  return result;
}

std::vector<types::qubit> QuantumCircuit::unmeasured_qubits() const {
  const auto measured = measured_qubits();
  const std::unordered_set<types::qubit> measured_set(measured.begin(), measured.end());

  std::vector<types::qubit> result;

  for (std::size_t i = 0; i < num_qubits_; ++i) {
    if (!measured_set.contains(i)) {
      result.push_back(i);
    }
  }

  return result;
}

}  // namespace quira
