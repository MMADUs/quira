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

#include "quira/circuit.hpp"
#include "quira/output.hpp"
#include "quira/register.hpp"
#include "quira/types.hpp"

#include <random>
#include <vector>

namespace quira {

/**
 * @brief Pure-state backend represented by a state vector.
 *
 * StateVector stores the complete computational-basis amplitude vector for an
 * n-qubit pure state. The vector dimension is always 2^n, and the state is kept
 * normalized after state-setting, gate application, and measurement collapse.
 */
class StateVector {
public:
  /**
   * @brief Creates the |0...0> state for the requested qubit count.
   *
   * @param num_qubits Number of qubits in the state.
   */
  explicit StateVector(Index num_qubits = 0);

  /**
   * @brief Creates a state vector from explicit amplitudes.
   *
   * @param amplitudes Complex amplitude vector. Its length must be a power of
   * two.
   *
   * @throws std::invalid_argument If the amplitude count is not a power of two.
   * @throws std::runtime_error If the state cannot be normalized.
   */
  explicit StateVector(Vector amplitudes);

  /**
   * @brief Returns the number of qubits represented by this state.
   */
  [[nodiscard]] Index num_qubits() const noexcept;

  /**
   * @brief Returns the Hilbert-space dimension.
   */
  [[nodiscard]] Index dimension() const noexcept;

  /**
   * @brief Returns the full amplitude vector.
   *
   * @return Read-only reference to the internal state vector.
   */
  [[nodiscard]] const Vector& state() const noexcept;

  /**
   * @brief Resets the backend to the |0...0> state.
   *
   * @param num_qubits Number of qubits in the new state.
   */
  void reset(Index num_qubits);

  /**
   * @brief Replaces the full state vector.
   *
   * @param amplitudes Complex amplitude vector. Its length must be a power of
   * two.
   *
   * @throws std::invalid_argument If the amplitude count is not a power of two.
   * @throws std::runtime_error If the state cannot be normalized.
   */
  void set_state(Vector amplitudes);

  /**
   * @brief Returns the probability of a computational basis state.
   *
   * @param basis_state Basis-state index.
   *
   * @throws std::out_of_range If basis_state is outside the state dimension.
   */
  [[nodiscard]] double probability(Index basis_state) const;

  /**
   * @brief Returns the Euclidean norm of the state vector.
   */
  [[nodiscard]] double norm() const;

  /**
   * @brief Normalizes the state vector in place.
   *
   * @throws std::runtime_error If the state has zero norm.
   */
  void normalize();

  /**
   * @brief Applies a local unitary to target qubits.
   *
   * @param unitary Matrix with dimension 2^k by 2^k, where k is targets.size().
   * @param targets Qubits acted on by the unitary.
   *
   * @throws std::invalid_argument If targets are empty, duplicated, or the
   * matrix dimension does not match the target count.
   * @throws std::out_of_range If any target qubit is outside the state.
   */
  void apply(const Matrix& unitary, const std::vector<Qubit>& targets);

  /**
   * @brief Measures one qubit and collapses the state.
   *
   * @param qubit Qubit to measure.
   * @param rng Random generator used for sampling.
   * @return Measured Boolean outcome.
   *
   * @throws std::out_of_range If qubit is outside the state.
   */
  bool measure(Qubit qubit, std::mt19937_64& rng);

  /**
   * @brief Measures all qubits and collapses to one computational basis state.
   *
   * @param rng Random generator used for sampling.
   * @return Measured outcomes in qubit index order.
   */
  std::vector<bool> measure_all(std::mt19937_64& rng);

  /**
   * @brief Resets one qubit to |0>.
   *
   * The reset is implemented as measurement followed by a Pauli-X correction
   * when the measured outcome is 1.
   *
   * @param qubit Qubit to reset.
   * @param rng Random generator used for measurement sampling.
   *
   * @throws std::out_of_range If qubit is outside the state.
   */
  void reset_qubit(Qubit qubit, std::mt19937_64& rng);

private:
  Index num_qubits_{};
  Vector state_;

  void validate_qubit(Qubit qubit) const;
  void validate_targets(const std::vector<Qubit>& targets) const;
};

/**
 * @brief State-vector circuit simulator.
 *
 * StateVectorSimulator executes QuantumCircuit instructions by creating a fresh
 * StateVector for each shot, applying gates in insertion order, and collecting
 * classical measurement outcomes into a RunResult.
 */
class StateVectorSimulator {
public:
  /**
   * @brief Creates a simulator with a non-deterministic random seed.
   */
  StateVectorSimulator();

  /**
   * @brief Creates a simulator with a deterministic random seed.
   *
   * @param seed Seed used by the simulator random generator.
   */
  explicit StateVectorSimulator(std::uint64_t seed);

  /**
   * @brief Runs a circuit for the requested number of shots.
   *
   * @param circuit Circuit to execute.
   * @param shots Number of independent shots.
   * @return Aggregated run result.
   *
   * @throws std::invalid_argument If circuit instructions contain invalid gate
   * dimensions or targets.
   * @throws std::out_of_range If an instruction references an invalid qubit or
   * classical bit.
   */
  SimulationOutput run(const QuantumCircuit& circuit, Index shots = 1024);

private:
  std::mt19937_64 rng_;

  ClassicalRegister run_single_shot(const QuantumCircuit& circuit);
};

}  // namespace quira
