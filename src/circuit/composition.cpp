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
#include "quira/unitary.hpp"

#include <stdexcept>
#include <unordered_set>

namespace quira {

namespace {

template<typename T>
void validate_distinct(const std::vector<T>& values, const char* message) {
  std::unordered_set<T> seen;

  for (T value : values) {
    if (!seen.insert(value).second) {
      throw std::invalid_argument(message);
    }
  }
}

std::vector<types::qubit> map_qubits(const std::vector<types::qubit>& source,
                                     const std::vector<types::qubit>& qubit_map,
                                     const char* caller) {
  std::vector<types::qubit> mapped;
  mapped.reserve(source.size());

  for (types::qubit qubit : source) {
    if (qubit >= qubit_map.size()) {
      throw exception::OutOfRange(caller, "qubit map is missing a source qubit");
    }

    mapped.push_back(qubit_map[qubit]);
  }

  return mapped;
}

std::vector<types::clbit> map_clbits(const std::vector<types::clbit>& source,
                                     const std::vector<types::clbit>& clbit_map,
                                     const char* caller) {
  std::vector<types::clbit> mapped;
  mapped.reserve(source.size());

  for (types::clbit clbit : source) {
    if (clbit >= clbit_map.size()) {
      throw exception::OutOfRange(caller, "classical bit map is missing a source bit");
    }

    mapped.push_back(clbit_map[clbit]);
  }

  return mapped;
}

}  // namespace

QuantumCircuit& QuantumCircuit::append(const QuantumCircuit& other) {
  return compose(other);
}

QuantumCircuit& QuantumCircuit::compose(const QuantumCircuit& other) {
  std::vector<types::qubit> qubit_map(other.num_qubits());
  std::vector<types::qubit> clbit_map(other.num_clbits());

  for (std::size_t i = 0; i < other.num_qubits(); ++i) {
    qubit_map[i] = i;
  }

  for (std::size_t i = 0; i < other.num_clbits(); ++i) {
    clbit_map[i] = i;
  }

  return compose(other, qubit_map, clbit_map);
}

QuantumCircuit& QuantumCircuit::compose(const QuantumCircuit& other,
                                        const std::vector<types::qubit>& qubit_map,
                                        const std::vector<types::clbit>& clbit_map) {
  if (qubit_map.size() != other.num_qubits()) {
    throw std::invalid_argument("Qubit map size must match source circuit");
  }

  if (clbit_map.size() != other.num_clbits()) {
    throw std::invalid_argument("Classical bit map size must match source circuit");
  }

  validate_qubits(qubit_map);
  validate_clbits(clbit_map);

  validate_distinct(qubit_map, "Qubit map contains duplicate targets");
  validate_distinct(clbit_map, "Classical bit map contains duplicate targets");

  for (const Instruction& instruction : other.instructions()) {
    std::visit(
        types::Overloaded{
            [&](const GateInstruction& item) {
              const QuantumGate& gate = item.gate();

              auto mapped_gate = std::make_unique<Unitary>(
                  gate.name(), gate.unitary(),
                  map_qubits(gate.targets(), qubit_map, "QuantumCircuit::compose()"));

              if (item.condition()) {
                const ClassicalCondition& condition = *item.condition();

                instructions_.emplace_back(GateInstruction{
                    std::move(mapped_gate),
                    ClassicalCondition{.clbit = clbit_map[condition.clbit],
                                       .value = condition.value}});
              } else {
                instructions_.emplace_back(GateInstruction{std::move(mapped_gate)});
              }
            },
            [&](const MeasureInstruction& item) {
              instructions_.emplace_back(MeasureInstruction{
                  .qubit = qubit_map[item.qubit], .clbit = clbit_map[item.clbit]});
            },
            [&](const MeasureBasisInstruction& item) {
              instructions_.emplace_back(MeasureBasisInstruction{
                  .basis = item.basis,
                  .qubits =
                      map_qubits(item.qubits, qubit_map, "QuantumCircuit::compose()"),
                  .clbits =
                      map_clbits(item.clbits, clbit_map, "QuantumCircuit::compose()")});
            },
            [&](const PostSelectInstruction& item) {
              instructions_.emplace_back(PostSelectInstruction{
                  .qubits =
                      map_qubits(item.qubits, qubit_map, "QuantumCircuit::compose()"),
                  .expected = item.expected});
            },
            [&](const ResetInstruction& item) {
              instructions_.emplace_back(
                  ResetInstruction{.qubit = qubit_map[item.qubit]});
            },
            [&](const BarrierInstruction&) {
              instructions_.emplace_back(BarrierInstruction{});
            },
        },
        instruction);
  }
  return *this;
}

QuantumCircuit QuantumCircuit::repeat(std::size_t count) const {
  QuantumCircuit result{num_qubits_, num_clbits_};

  for (std::size_t i = 0; i < count; ++i) {
    result.compose(*this);
  }

  return result;
}

}  // namespace quira
