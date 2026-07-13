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
 * A conditional gate is applied only when the selected classical bit matches
 * the requested Boolean value.
 */
struct ClassicalCondition {
  Clbit clbit{};
  bool value{};
};

/**
 * @brief Gate operation recorded by a quantum circuit.
 *
 * GateInstruction owns a cloned gate through exclusive polymorphic ownership.
 * This allows QuantumCircuit to store user-defined gate types without knowing
 * their concrete class.
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
   * @param condition Classical condition required for applying the gate.
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
 * @brief Measurement instruction.
 *
 * Records a measurement from one qubit into one classical bit.
 */
struct MeasureInstruction {
  Qubit qubit{};
  Clbit clbit{};
};

/**
 * @brief Reset instruction.
 *
 * Records a reset operation that returns a qubit to the |0> state.
 */
struct ResetInstruction {
  Qubit qubit{};
};

/**
 * @brief Barrier instruction.
 *
 * Barriers preserve circuit structure and may be used by later compiler,
 * visualization, or scheduling passes.
 */
struct BarrierInstruction {};

/**
 * @brief Instruction type stored by QuantumCircuit.
 *
 * Instructions are recorded in insertion order and executed by simulators in
 * the same order unless an optimizer/transpiler rewrites the circuit.
 */
using Instruction = std::variant<GateInstruction, MeasureInstruction, ResetInstruction,
                                 BarrierInstruction>;

/**
 * @brief Ordered description of a quantum circuit.
 *
 * QuantumCircuit records gate, measurement, reset, and barrier instructions.
 * It validates important API boundaries but does not execute or simulate the
 * circuit.
 */
class QuantumCircuit {
public:
  /**
   * @brief Creates a quantum circuit with quantum and classical bit capacity.
   *
   * @param num_qubits Number of qubits in the circuit.
   * @param num_clbits Number of classical bits in the circuit.
   */
  QuantumCircuit(Index num_qubits, Index num_clbits = 0);

  /**
   * @brief Returns the number of qubits.
   */
  Index num_qubits() const noexcept;

  /**
   * @brief Returns the number of classical bits.
   */
  Index num_clbits() const noexcept;

  /**
   * @brief Returns the recorded circuit instructions.
   *
   * @return Read-only instruction sequence.
   */
  const std::vector<Instruction>& instructions() const noexcept;

  /**
   * @brief Appends a gate by cloning it into the circuit.
   *
   * @param gate Gate to append.
   * @return Reference to this circuit.
   *
   * @throws std::invalid_argument If the gate has invalid targets.
   * @throws std::out_of_range If any target qubit is outside the circuit.
   */
  QuantumCircuit& add(const QuantumGate& gate);

  /**
   * @brief Appends a gate by transferring ownership into the circuit.
   *
   * @param gate Gate to append.
   * @return Reference to this circuit.
   *
   * @throws std::invalid_argument If gate is null or has invalid targets.
   * @throws std::out_of_range If any target qubit is outside the circuit.
   */
  QuantumCircuit& add(std::unique_ptr<QuantumGate> gate);

  /**
   * @brief Appends a classically conditional gate by cloning it.
   *
   * @param clbit Classical bit used by the condition.
   * @param value Required classical bit value.
   * @param gate Gate to append.
   * @return Reference to this circuit.
   *
   * @throws std::out_of_range If clbit or a target qubit is outside the
   * circuit.
   * @throws std::invalid_argument If the gate has invalid targets.
   */
  QuantumCircuit& conditional_add(Clbit clbit, bool value, const QuantumGate& gate);

  /**
   * @brief Appends a classically conditional gate by transferring ownership.
   *
   * @param clbit Classical bit used by the condition.
   * @param value Required classical bit value.
   * @param gate Gate to append.
   * @return Reference to this circuit.
   *
   * @throws std::invalid_argument If gate is null or has invalid targets.
   * @throws std::out_of_range If clbit or a target qubit is outside the
   * circuit.
   */
  QuantumCircuit& conditional_add(Clbit clbit, bool value,
                                  std::unique_ptr<QuantumGate> gate);

  /**
   * @brief Records a measurement from a qubit into a classical bit.
   *
   * @param qubit Qubit to measure.
   * @param clbit Classical bit that receives the measurement outcome.
   * @return Reference to this circuit.
   *
   * @throws std::out_of_range If qubit or clbit is outside the circuit.
   */
  QuantumCircuit& measure(Qubit qubit, Clbit clbit);

  /**
   * @brief Records measurement of every qubit into same-index classical bits.
   *
   * @return Reference to this circuit.
   *
   * @throws std::out_of_range If the circuit has fewer classical bits than
   * qubits.
   */
  QuantumCircuit& measure_all();

  /**
   * @brief Records a reset operation for one qubit.
   *
   * @param qubit Qubit to reset.
   * @return Reference to this circuit.
   *
   * @throws std::out_of_range If qubit is outside the circuit.
   */
  QuantumCircuit& reset(Qubit qubit);

  /**
   * @brief Records a barrier instruction.
   *
   * @return Reference to this circuit.
   */
  QuantumCircuit& barrier();

private:
  Index num_qubits_{};
  Index num_clbits_{};
  std::vector<Instruction> instructions_;

  void validate_qubit(Qubit qubit) const;
  void validate_clbit(Clbit clbit) const;
  void validate_gate(const QuantumGate& gate) const;
};

}  // namespace quira
