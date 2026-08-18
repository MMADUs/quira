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

#include "quira/circuit.hpp"
#include "quira/output.hpp"
#include "quira/random.hpp"
#include "quira/register.hpp"
#include "quira/types.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace quira {

/**
 * @brief Pure-state backend represented by a state vector.
 *
 * StateVector stores the full computational-basis amplitude vector for an
 * n-qubit pure state. Qubit indexing follows Quira's little-endian convention:
 * qubit 0 corresponds to the least-significant bit of a basis index.
 */
class StateVector {
public:
  /**
   * @brief Creates the |0...0> state for the requested qubit count.
   *
   * @param num_qubits Number of qubits in the state.
   *
   * @throws std::out_of_range If the state dimension cannot fit in Index.
   */
  explicit StateVector(std::size_t num_qubits = 0);

  /**
   * @brief Creates a state vector from explicit amplitudes.
   *
   * @param amplitudes Complex amplitude vector. Its length must be a power of
   * two.
   *
   * @throws std::invalid_argument If the amplitude count is not a power of two.
   * @throws std::runtime_error If the state cannot be normalized.
   */
  explicit StateVector(types::ket amplitudes);

  /**
   * @brief Returns the number of qubits represented by this state.
   */
  [[nodiscard]] std::size_t num_qubits() const noexcept;

  /**
   * @brief Returns the Hilbert-space dimension.
   */
  [[nodiscard]] std::size_t dimension() const noexcept;

  /**
   * @brief Returns the full amplitude vector.
   *
   * @return Read-only reference to the internal state vector.
   */
  [[nodiscard]] const types::ket& state() const noexcept;

  /**
   * @brief Resets the backend to the |0...0> state.
   *
   * @param num_qubits Number of qubits in the new state.
   *
   * @throws std::out_of_range If the state dimension cannot fit in Index.
   */
  void reset(std::size_t num_qubits);

  /**
   * @brief Replaces the full state vector.
   *
   * @param amplitudes Complex amplitude vector. Its length must be a power of
   * two.
   *
   * @throws std::invalid_argument If the amplitude count is not a power of two.
   * @throws std::runtime_error If the state cannot be normalized.
   */
  void set_state(types::ket amplitudes);

  /**
   * @brief Returns the probability of one computational basis state.
   *
   * @param basis_state Basis-state index.
   * @return Probability equal to the squared magnitude of the amplitude.
   *
   * @throws std::out_of_range If basis_state is outside the state dimension.
   */
  [[nodiscard]] types::real_n probability(std::size_t basis_state) const;

  /**
   * @brief Returns marginal probabilities over selected qubits.
   *
   * @param targets Qubits to inspect.
   * @return Probability vector of size 2^targets.size().
   *
   * @note The returned probability index uses targets[0] as local bit 0.
   *
   * @throws std::invalid_argument If targets are empty or duplicated.
   * @throws std::out_of_range If any target qubit is outside the state.
   */
  [[nodiscard]] std::vector<types::real_n> probabilities(
      const std::vector<types::qubit>& targets) const;

  /**
   * @brief Returns the Euclidean norm of the state vector.
   */
  [[nodiscard]] types::real_n norm() const;

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
   * @note targets[0] maps to local bit 0 in the unitary matrix.
   *
   * @throws std::invalid_argument If targets are empty, duplicated, or the
   * matrix dimension does not match the target count.
   * @throws std::out_of_range If any target qubit is outside the state.
   */
  void apply(const types::c_mat& unitary, const std::vector<types::qubit>& targets);

  /**
   * @brief Samples selected qubits without collapsing the state.
   *
   * @param targets Qubits to sample.
   * @param rng Random engine used for measurement sampling.
   * @return Sampled Boolean outcomes in target order.
   *
   * @note This does not mutate the state vector.
   *
   * @throws std::invalid_argument If targets are empty or duplicated.
   * @throws std::out_of_range If any target qubit is outside the state.
   */
  [[nodiscard]] std::vector<bool> sample(const std::vector<types::qubit>& targets,
                                         RandomEngine& rng) const;

  /**
   * @brief Measures one qubit in the standard Z basis and collapses the state.
   *
   * @param qubit Qubit to measure.
   * @param rng Random engine used for measurement sampling.
   * @return Measured Boolean outcome.
   *
   * @throws std::out_of_range If qubit is outside the state.
   */
  bool measure(types::qubit qubit, RandomEngine& rng);

  /**
   * @brief Measures selected qubits in the standard Z basis.
   *
   * @param qubits Qubits to measure.
   * @param rng Random engine used for measurement sampling.
   * @return Measured Boolean outcomes in qubit order.
   *
   * @throws std::invalid_argument If qubits are empty or duplicated.
   * @throws std::out_of_range If any qubit is outside the state.
   */
  std::vector<bool> measure(const std::vector<types::qubit>& qubits, RandomEngine& rng);

  /**
   * @brief Measures selected qubits in a custom basis.
   *
   * @param basis Unitary basis matrix. Columns are measurement basis vectors.
   * @param qubits Qubits to measure.
   * @param rng Random engine used for measurement sampling.
   * @return Measured Boolean outcomes in qubit order.
   *
   * @throws std::invalid_argument If qubits are empty, duplicated, basis has an
   * invalid dimension, or basis is not unitary.
   * @throws std::out_of_range If any qubit is outside the state.
   */
  std::vector<bool> measure_basis(const types::c_mat& basis,
                                  const std::vector<types::qubit>& qubits,
                                  RandomEngine& rng);

  /**
   * @brief Measures all qubits in the standard Z basis.
   *
   * @param rng Random engine used for measurement sampling.
   * @return Measured outcomes in qubit index order.
   */
  std::vector<bool> measure_all(RandomEngine& rng);

  /**
   * @brief Resets one qubit to |0>.
   *
   * Reset is implemented as measurement followed by a Pauli-X correction when
   * the measured outcome is 1.
   *
   * @param qubit Qubit to reset.
   * @param rng Random engine used for measurement sampling.
   *
   * @throws std::out_of_range If qubit is outside the state.
   */
  void reset_qubit(types::qubit qubit, RandomEngine& rng);

  /**
   * @brief Resets selected qubits to |0>.
   *
   * @param qubits Qubits to reset.
   * @param rng Random engine used for measurement sampling.
   *
   * @throws std::invalid_argument If qubits are empty or duplicated.
   * @throws std::out_of_range If any qubit is outside the state.
   */
  void reset_qubits(const std::vector<types::qubit>& qubits, RandomEngine& rng);

private:
  std::size_t num_qubits_{};
  types::ket state_;

  void collapse(const std::vector<types::qubit>& targets,
                const std::vector<bool>& outcomes);
  void validate_qubit(types::qubit qubit, const char* caller = nullptr) const;
  void validate_targets(const std::vector<types::qubit>& targets,
                        const char* caller = nullptr) const;
};

/**
 * @brief State-vector circuit simulator.
 *
 * StateVectorSimulator executes QuantumCircuit instructions by creating a fresh
 * StateVector for each requested shot. Post-selection may reject shots, in
 * which case the rejected shot is counted but no ClassicalRegister is stored
 * for it.
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
   * @param shots Number of requested shots.
   * @return Aggregated simulation output.
   *
   * @throws std::invalid_argument If circuit instructions contain invalid gate,
   * measurement, or post-selection data.
   * @throws std::out_of_range If an instruction references an invalid qubit or
   * classical bit.
   */
  SimulationOutput run(const QuantumCircuit& circuit, std::size_t shots = 1);

  /**
   * @brief Evolves a circuit into a final pure state without sampling.
   *
   * Only unitary gate instructions and barriers are supported. Measurement,
   * reset, post-selection, and conditional gates are rejected.
   */
  [[nodiscard]] StateVector evolve(const QuantumCircuit& circuit) const;

private:
  RandomEngine rng_;

  std::optional<ClassicalRegister> run_single_shot(const QuantumCircuit& circuit);
};

}  // namespace quira
