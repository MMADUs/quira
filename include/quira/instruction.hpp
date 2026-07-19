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

#pragma once

#include "quira/gate.hpp"
#include "quira/types.hpp"

#include <memory>
#include <optional>
#include <variant>
#include <vector>

namespace quira {
/**
 * @brief Classical condition attached to a gate instruction.
 *
 * A conditional gate is applied only when the selected classical bit has already
 * been measured and matches the requested value.
 */
struct ClassicalCondition {
  types::clbit clbit{};
  bool value{};
};

/**
 * @brief Gate operation recorded by a quantum circuit.
 *
 * GateInstruction owns a gate through exclusive polymorphic ownership. Copying
 * a GateInstruction clones the underlying concrete gate.
 */
class GateInstruction {
public:
  /**
   * @brief Creates an unconditional gate instruction.
   *
   * @param gate Gate object owned by the instruction.
   *
   * @throws std::invalid_argument If gate is null.
   */
  explicit GateInstruction(std::unique_ptr<QuantumGate> gate);

  /**
   * @brief Creates a classically conditional gate instruction.
   *
   * @param gate Gate object owned by the instruction.
   * @param condition Classical condition required before applying the gate.
   *
   * @throws std::invalid_argument If gate is null.
   */
  GateInstruction(std::unique_ptr<QuantumGate> gate, ClassicalCondition condition);

  ~GateInstruction() = default;

  GateInstruction(const GateInstruction& other);
  GateInstruction& operator=(const GateInstruction& other);

  GateInstruction(GateInstruction&&) noexcept = default;
  GateInstruction& operator=(GateInstruction&&) noexcept = default;

  /**
   * @brief Returns the recorded gate.
   *
   * @return Read-only reference to the owned gate.
   */
  [[nodiscard]] const QuantumGate& gate() const noexcept;

  /**
   * @brief Returns the optional classical condition.
   *
   * @return Empty optional for unconditional gates.
   */
  [[nodiscard]] const std::optional<ClassicalCondition>& condition() const noexcept;

private:
  std::unique_ptr<QuantumGate> gate_;
  std::optional<ClassicalCondition> condition_;
};

/**
 * @brief Standard Z-basis measurement instruction.
 *
 * Records measurement of one qubit into one classical bit. The simulator
 * collapses the qubit in the computational basis and writes the Boolean outcome
 * into clbit.
 */
struct MeasureInstruction {
  types::qubit qubit{};
  types::clbit clbit{};
};

/**
 * @brief Custom-basis measurement instruction.
 *
 * The basis matrix columns represent the measurement basis vectors. For k
 * measured qubits, basis must have dimension 2^k by 2^k and must be unitary.
 *
 * Simulation semantics are:
 * - apply basis dagger to the target qubits,
 * - measure in the standard Z basis,
 * - rotate the collapsed state back by applying basis.
 *
 * Outcomes are written to clbits in the same order as qubits.
 */
struct MeasureBasisInstruction {
  types::c_mat basis;
  std::vector<types::qubit> qubits;
  std::vector<types::clbit> clbits;
};

/**
 * @brief Post-selection measurement instruction.
 *
 * Post-selection measures qubits in the standard Z basis and accepts the shot
 * only when the measured outcomes match expected. Rejected shots are not added
 * to SimulationOutput counts.
 *
 * @note This instruction is a filter only. It does not write outcomes into
 * classical bits.
 */
struct PostSelectInstruction {
  std::vector<types::qubit> qubits;
  std::vector<bool> expected;
};

/**
 * @brief Reset instruction.
 *
 * Records a reset operation that returns one qubit to the |0> state.
 */
struct ResetInstruction {
  types::qubit qubit{};
};

/**
 * @brief Barrier instruction.
 *
 * Barriers preserve circuit structure and may be used by later compiler,
 * visualization, or scheduling passes. They do not affect simulation state.
 */
struct BarrierInstruction {};

/**
 * @brief Circuit instruction stored in insertion order.
 */
using Instruction =
    std::variant<GateInstruction, MeasureInstruction, MeasureBasisInstruction,
                 PostSelectInstruction, ResetInstruction, BarrierInstruction>;

}  // namespace quira
