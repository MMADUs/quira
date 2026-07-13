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

#include "quira/backend/statevec.hpp"
#include "quira/output.hpp"

#include <optional>
#include <type_traits>
#include <variant>
#include <vector>

namespace quira {

StateVectorSimulator::StateVectorSimulator() : rng_(std::random_device{}()) {
}

StateVectorSimulator::StateVectorSimulator(std::uint64_t seed) : rng_(seed) {
}

/**
 * Runs a circuit for the requested number of shots.
 *
 * Each shot starts from |0...0>, applies the recorded circuit instructions,
 * and stores classical measurement outcomes.
 */
SimulationOutput StateVectorSimulator::run(const QuantumCircuit& circuit, Index shots) {
  std::vector<ClassicalRegister> shot_registers;
  shot_registers.reserve(shots);

  for (Index shot = 0; shot < shots; ++shot) {
    shot_registers.push_back(run_single_shot(circuit));
  }

  return SimulationOutput{std::move(shot_registers), shots};
}

ClassicalRegister StateVectorSimulator::run_single_shot(const QuantumCircuit& circuit) {
  StateVector state{circuit.num_qubits()};
  ClassicalRegister classical{circuit.num_clbits()};

  for (const Instruction& instruction : circuit.instructions()) {
    std::visit(
        [&](const auto& item) {
          using T = std::decay_t<decltype(item)>;

          if constexpr (std::is_same_v<T, GateInstruction>) {
            const auto& condition = item.condition();
            if (condition) {
              const std::optional<bool> measured_value =
                  classical.get_raw(condition->clbit);

              if (!measured_value || *measured_value != condition->value) {
                return;
              }
            }

            const QuantumGate& gate = item.gate();
            state.apply(gate.unitary(), gate.targets());
          }

          if constexpr (std::is_same_v<T, MeasureInstruction>) {
            classical.set(item.clbit, state.measure(item.qubit, rng_));
          }

          if constexpr (std::is_same_v<T, ResetInstruction>) {
            state.reset_qubit(item.qubit, rng_);
          }

          if constexpr (std::is_same_v<T, BarrierInstruction>) {
            return;
          }
        },
        instruction);
  }

  return classical;
}

}  // namespace quira
