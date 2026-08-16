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

#include "quira/types.hpp"

#include <cstdint>

namespace quira::states {

/**
 * @brief Bell-state variant.
 *
 * The computational basis follows Quira's little-endian basis-index convention.
 */
enum class BellState : std::uint8_t {
  PhiPlus,
  PhiMinus,
  PsiPlus,
  PsiMinus,
};

/**
 * @brief Returns the computational zero state |0...0>.
 *
 * @param num_qubits Number of qubits.
 * @return State vector with amplitude 1 at basis index 0.
 *
 * @throws std::out_of_range If the state dimension cannot fit in Index.
 */
[[nodiscard]] types::ket zero_state(std::size_t num_qubits);

/**
 * @brief Returns the computational one state |1...1>.
 *
 * @param num_qubits Number of qubits.
 * @return State vector with amplitude 1 at the final basis index.
 *
 * @throws std::out_of_range If the state dimension cannot fit in Index.
 */
[[nodiscard]] types::ket one_state(std::size_t num_qubits);

/**
 * @brief Returns a computational basis state.
 *
 * @param num_qubits Number of qubits.
 * @param basis_index Computational basis index.
 * @return State vector with amplitude 1 at basis_index.
 *
 * @throws std::out_of_range If basis_index is outside the state dimension.
 */
[[nodiscard]] types::ket basis_state(std::size_t num_qubits, std::size_t basis_index);

/**
 * @brief Returns the equal-superposition state |+>^n.
 *
 * @param num_qubits Number of qubits.
 * @return State vector with equal amplitudes over all basis states.
 *
 * @throws std::out_of_range If the state dimension cannot fit in Index.
 */
[[nodiscard]] types::ket plus_state(std::size_t num_qubits);

/**
 * @brief Returns the repeated minus state |->^n.
 *
 * @param num_qubits Number of qubits.
 * @return State vector whose basis amplitudes have signs from basis parity.
 *
 * @throws std::out_of_range If the state dimension cannot fit in Index.
 */
[[nodiscard]] types::ket minus_state(std::size_t num_qubits);

/**
 * @brief Returns a single-qubit X-basis eigenstate.
 *
 * @param value False returns |+>, true returns |->.
 * @return X-basis state vector.
 */
[[nodiscard]] types::ket x_basis_state(bool value);

/**
 * @brief Returns a single-qubit Y-basis eigenstate.
 *
 * @param value False returns |y+>, true returns |y->.
 * @return Y-basis state vector.
 */
[[nodiscard]] types::ket y_basis_state(bool value);

/**
 * @brief Returns a two-qubit Bell state (default = (|00> + |11>) / sqrt(2))
 *
 * @param state Bell-state variant.
 * @return Bell state vector.
 */
[[nodiscard]] types::ket bell_state(BellState state = BellState::PhiPlus);

/**
 * @brief Returns the n-qubit GHZ state.
 *
 * The GHZ state is (|0...0> + |1...1>) / sqrt(2).
 *
 * @param num_qubits Number of qubits.
 * @return GHZ state vector.
 *
 * @throws std::invalid_argument If num_qubits is less than 2.
 * @throws std::out_of_range If the state dimension cannot fit in Index.
 */
[[nodiscard]] types::ket ghz_state(std::size_t num_qubits);

/**
 * @brief Returns the n-qubit W state.
 *
 * The W state is the equal superposition of all basis states with exactly one
 * qubit in the |1> state.
 *
 * @param num_qubits Number of qubits.
 * @return W state vector.
 *
 * @throws std::invalid_argument If num_qubits is zero.
 * @throws std::out_of_range If the state dimension cannot fit in Index.
 */
[[nodiscard]] types::ket w_state(std::size_t num_qubits);

/**
 * @brief Returns the projector for the computational zero state.
 *
 * @param num_qubits Number of qubits.
 * @return Projector |0...0><0...0|.
 */
[[nodiscard]] types::c_mat zero_projector(std::size_t num_qubits);

/**
 * @brief Returns the projector for the computational one state.
 *
 * @param num_qubits Number of qubits.
 * @return Projector |1...1><1...1|.
 */
[[nodiscard]] types::c_mat one_projector(std::size_t num_qubits);

/**
 * @brief Returns the projector for a computational basis state.
 *
 * @param num_qubits Number of qubits.
 * @param basis_index Computational basis index.
 * @return Projector |basis_index><basis_index|.
 */
[[nodiscard]] types::c_mat basis_projector(std::size_t num_qubits,
                                           std::size_t basis_index);

/**
 * @brief Returns the projector for |+>^n.
 *
 * @param num_qubits Number of qubits.
 * @return Projector onto the plus state.
 */
[[nodiscard]] types::c_mat plus_projector(std::size_t num_qubits);

/**
 * @brief Returns the projector for |->^n.
 *
 * @param num_qubits Number of qubits.
 * @return Projector onto the minus state.
 */
[[nodiscard]] types::c_mat minus_projector(std::size_t num_qubits);

/**
 * @brief Returns the projector for a single-qubit X-basis state.
 *
 * @param value False returns |+><+|, true returns |-><-|.
 * @return X-basis projector.
 */
[[nodiscard]] types::c_mat x_basis_projector(bool value);

/**
 * @brief Returns the projector for a single-qubit Y-basis state.
 *
 * @param value False returns |y+><y+|, true returns |y-><y-|.
 * @return Y-basis projector.
 */
[[nodiscard]] types::c_mat y_basis_projector(bool value);

/**
 * @brief Returns the projector for a Bell state.
 *
 * @param state Bell-state variant.
 * @return Bell-state projector.
 */
[[nodiscard]] types::c_mat bell_projector(BellState state = BellState::PhiPlus);

/**
 * @brief Returns the projector for the GHZ state.
 *
 * @param num_qubits Number of qubits.
 * @return GHZ-state projector.
 */
[[nodiscard]] types::c_mat ghz_projector(std::size_t num_qubits);

/**
 * @brief Returns the projector for the W state.
 *
 * @param num_qubits Number of qubits.
 * @return W-state projector.
 */
[[nodiscard]] types::c_mat w_projector(std::size_t num_qubits);

}  // namespace quira::states
