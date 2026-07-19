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
#include "quira/instruction.hpp"
#include "quira/register.hpp"
#include "quira/types.hpp"

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace quira {

/**
 * @brief Ordered description of a quantum circuit.
 *
 * QuantumCircuit records gates, measurements, post-selection, resets, and
 * barriers. It validates important API boundaries but does not execute the
 * circuit.
 */
class QuantumCircuit {
public:
  /**
   * @brief Creates a quantum circuit.
   *
   * @param num_qubits Number of qubits in the circuit.
   * @param num_clbits Number of classical bits in the circuit.
   */
  QuantumCircuit(std::size_t num_qubits, std::size_t num_clbits);

  /**
   * @brief Creates a quantum circuit from quantum and classical registers.
   *
   * The circuit uses the register sizes as its qubit and classical-bit capacity.
   * Register labels and per-bit values are not copied into the circuit.
   *
   * @param qreg Quantum register defining the number of qubits.
   * @param creg Classical register defining the number of classical bits.
   */
  QuantumCircuit(const QuantumRegister& qreg, const ClassicalRegister& creg);

  /**
   * @brief Returns the number of qubits.
   */
  [[nodiscard]] std::size_t num_qubits() const noexcept;

  /**
   * @brief Returns the number of classical bits.
   */
  [[nodiscard]] std::size_t num_clbits() const noexcept;

  /**
   * @brief Returns the recorded circuit instructions.
   *
   * @return Read-only instruction sequence.
   */
  [[nodiscard]] const std::vector<Instruction>& instructions() const noexcept;

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
   * @throws std::out_of_range If clbit or a target qubit is outside the circuit.
   * @throws std::invalid_argument If the gate has invalid targets.
   */
  QuantumCircuit& conditional_add(types::clbit clbit, bool value,
                                  const QuantumGate& gate);

  /**
   * @brief Appends a classically conditional gate by transferring ownership.
   *
   * @param clbit Classical bit used by the condition.
   * @param value Required classical bit value.
   * @param gate Gate to append.
   * @return Reference to this circuit.
   *
   * @throws std::invalid_argument If gate is null or has invalid targets.
   * @throws std::out_of_range If clbit or a target qubit is outside the circuit.
   */
  QuantumCircuit& conditional_add(types::clbit clbit, bool value,
                                  std::unique_ptr<QuantumGate> gate);

  /**
   * @brief Records one standard Z-basis measurement.
   *
   * @param qubit Qubit to measure.
   * @param clbit Classical bit that receives the measurement outcome.
   * @return Reference to this circuit.
   *
   * @throws std::out_of_range If qubit or clbit is outside the circuit.
   */
  QuantumCircuit& measure(types::qubit qubit, types::qubit clbit);

  /**
   * @brief Records multiple standard Z-basis measurements.
   *
   * @param qubits Qubits to measure.
   * @param clbits Classical bits receiving the outcomes.
   * @return Reference to this circuit.
   *
   * @throws std::invalid_argument If vectors are empty, have different sizes,
   * or contain duplicate qubits or classical bits.
   * @throws std::out_of_range If any qubit or classical bit is outside the
   * circuit.
   */
  QuantumCircuit& measure(const std::vector<types::qubit>& qubits,
                          const std::vector<types::clbit>& clbits);

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
   * @brief Records custom-basis measurement.
   *
   * @param basis Unitary basis matrix. Columns are measurement basis vectors.
   * @param qubits Qubits measured by the basis.
   * @param clbits Classical bits receiving the outcomes.
   * @return Reference to this circuit.
   *
   * @throws std::invalid_argument If vectors are empty, have different sizes,
   * contain duplicates, or if basis is not unitary with dimension 2^k by 2^k.
   * @throws std::out_of_range If any qubit or classical bit is outside the
   * circuit.
   */
  QuantumCircuit& measure_basis(const types::c_mat& basis,
                                const std::vector<types::qubit>& qubits,
                                const std::vector<types::clbit>& clbits);

  /**
   * @brief Records one standard Z-basis post-selection measurement.
   *
   * @param qubit Qubit to measure.
   * @param expected Required outcome for the shot to be accepted.
   * @return Reference to this circuit.
   *
   * @throws std::out_of_range If qubit is outside the circuit.
   */
  QuantumCircuit& post_select(types::qubit qubit, bool expected);

  /**
   * @brief Records multiple standard Z-basis post-selection measurements.
   *
   * @param qubits Qubits to measure.
   * @param expected Required outcomes in qubit order.
   * @return Reference to this circuit.
   *
   * @throws std::invalid_argument If vectors are empty, have different sizes,
   * or qubits contain duplicates.
   * @throws std::out_of_range If any qubit is outside the circuit.
   */
  QuantumCircuit& post_select(const std::vector<types::qubit>& qubits,
                              const std::vector<bool>& expected);

  /**
   * @brief Records a reset operation for one qubit.
   *
   * @param qubit Qubit to reset.
   * @return Reference to this circuit.
   *
   * @throws std::out_of_range If qubit is outside the circuit.
   */
  QuantumCircuit& reset(types::qubit qubit);

  /**
   * @brief Records a barrier instruction.
   *
   * @return Reference to this circuit.
   */
  QuantumCircuit& barrier();

  /**
   * @brief Appends another circuit to this circuit.
   *
   * This is equivalent to compose(other) with identity qubit and classical-bit
   * mapping.
   *
   * @param other Circuit to append.
   * @return Reference to this circuit.
   */
  QuantumCircuit& append(const QuantumCircuit& other);

  /**
   * @brief Composes another circuit onto this circuit using identity mapping.
   *
   * @param other Circuit to compose.
   * @return Reference to this circuit.
   *
   * @throws std::out_of_range If other uses qubits or classical bits that do
   * not exist in this circuit.
   */
  QuantumCircuit& compose(const QuantumCircuit& other);

  /**
   * @brief Composes another circuit onto this circuit with explicit remapping.
   *
   * @param other Circuit to compose.
   * @param qubit_map Maps each source qubit in other to a qubit in this circuit.
   * @param clbit_map Maps each source classical bit in other to a classical bit
   * in this circuit.
   * @return Reference to this circuit.
   *
   * @throws std::invalid_argument If mapping sizes do not match other, or if
   * mappings contain duplicate target bits.
   * @throws std::out_of_range If any mapped qubit or classical bit is outside
   * this circuit.
   */
  QuantumCircuit& compose(const QuantumCircuit& other,
                          const std::vector<types::qubit>& qubit_map,
                          const std::vector<types::clbit>& clbit_map);

  /**
   * @brief Returns a circuit containing this circuit repeated count times.
   *
   * @param count Number of repetitions.
   * @return New circuit with the same qubit and classical-bit capacity.
   *
   * @note count == 0 returns an empty circuit with the same capacity.
   */
  [[nodiscard]] QuantumCircuit repeat(std::size_t count) const;

  /**
   * @brief Returns the adjoint of this circuit.
   *
   * The adjoint reverses instruction order and replaces each unitary gate by
   * its conjugate-transpose MatrixGate.
   *
   * @return New circuit representing the unitary adjoint.
   *
   * @throws std::invalid_argument If the circuit contains non-unitary
   * instructions such as measurement, reset, post-selection, or conditional
   * gates.
   */
  [[nodiscard]] QuantumCircuit adjoint() const;

  /**
   * @brief Returns the inverse of this circuit.
   *
   * For unitary circuits, inverse() is equivalent to adjoint().
   *
   * @return New circuit representing the inverse unitary.
   *
   * @throws std::invalid_argument If adjoint() would reject the circuit.
   */
  [[nodiscard]] QuantumCircuit inverse() const;

  /**
   * @brief Estimates circuit depth from qubit usage.
   *
   * Instructions touching disjoint qubits can occupy the same layer. Barriers
   * force a new layer across all qubits.
   *
   * @return Estimated circuit depth.
   */
  [[nodiscard]] std::size_t depth() const;

  /**
   * @brief Counts gate instructions.
   *
   * @return Number of GateInstruction entries in the circuit.
   */
  [[nodiscard]] std::size_t gate_count() const;

  /**
   * @brief Counts gate instructions by name.
   *
   * @return Map from gate name to count.
   */
  [[nodiscard]] std::unordered_map<std::string, std::size_t> gate_counts() const;

  /**
   * @brief Returns qubits that appear in measurement-like instructions.
   *
   * This includes standard measurement, custom-basis measurement, and
   * post-selection.
   *
   * @return Sorted measured qubit indices.
   */
  [[nodiscard]] std::vector<types::qubit> measured_qubits() const;

  /**
   * @brief Returns qubits that do not appear in measurement-like instructions.
   *
   * @return Sorted unmeasured qubit indices.
   */
  [[nodiscard]] std::vector<types::qubit> unmeasured_qubits() const;

  /**
   * @brief Checks whether the circuit satisfies its core invariants.
   *
   * This verifies instruction indices, duplicate gate targets, duplicate
   * measurement targets, classical conditions, and custom measurement bases.
   *
   * @return True when the circuit is internally valid, false otherwise.
   */
  [[nodiscard]] bool is_valid() const;

private:
  std::size_t num_qubits_{};
  std::size_t num_clbits_{};
  std::vector<Instruction> instructions_;

  // circuit helper function
  void validate_qubit(types::qubit qubit, const char* caller) const;
  void validate_clbit(types::clbit clbit, const char* caller) const;
  void validate_qubits(const std::vector<types::qubit>& qubits,
                       const char* caller) const;
  void validate_clbits(const std::vector<types::clbit>& clbits,
                       const char* caller) const;
  void validate_gate(const QuantumGate& gate, const char* caller) const;
  void validate_measurement_lists(const std::vector<types::qubit>& qubits,
                                  const std::vector<types::clbit>& clbits,
                                  const char* caller) const;
};

}  // namespace quira
