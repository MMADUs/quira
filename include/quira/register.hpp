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

#include "quira/types.hpp"

#include <optional>
#include <string>
#include <vector>

namespace quira {

/**
 * @brief Classical bit used to store a measurement outcome.
 *
 * A ClassicalBit has a fixed index and an optional Boolean value. The value is
 * empty until a measurement writes an outcome into the bit.
 */
class ClassicalBit {
public:
  /**
   * @brief Creates an unset classical bit.
   *
   * @param index Position of this bit in its classical register.
   */
  explicit ClassicalBit(Clbit index);

  /**
   * @brief Returns this bit's register index.
   */
  Clbit index() const noexcept;

  /**
   * @brief Checks whether this bit contains a measurement outcome.
   */
  bool has_value() const noexcept;

  /**
   * @brief Returns the stored measurement outcome.
   *
   * @return Stored Boolean measurement value.
   *
   * @throws std::runtime_error If the bit has no stored outcome.
   */
  bool value() const;

  /**
   * @brief Returns the optional stored measurement outcome.
   *
   * @return Empty optional when the bit has not been measured.
   */
  std::optional<bool> raw_value() const noexcept;

  /**
   * @brief Stores a measurement outcome.
   *
   * @param value Boolean measurement value.
   */
  void set(bool value);

  /**
   * @brief Clears the stored measurement outcome.
   */
  void clear() noexcept;

private:
  Clbit index_{};
  std::optional<bool> value_;
};

/**
 * @brief Classical register used for per-shot measurement storage.
 *
 * ClassicalRegister owns a fixed-size sequence of classical bits. Simulators
 * write measurement outcomes into this register while executing a circuit shot.
 */
class ClassicalRegister {
public:
  /**
   * @brief Creates a classical register with the requested bit count.
   *
   * @param size Number of classical bits in the register.
   */
  explicit ClassicalRegister(Index size = 0);

  /**
   * @brief Returns the number of classical bits.
   */
  Index size() const noexcept;

  /**
   * @brief Checks whether a classical bit contains a measurement outcome.
   *
   * @param clbit Classical bit index.
   *
   * @throws std::out_of_range If clbit is outside the register.
   */
  bool has_value(Clbit clbit) const;

  /**
   * @brief Returns a measured classical bit value.
   *
   * @param clbit Classical bit index.
   * @return Stored Boolean measurement value.
   *
   * @throws std::out_of_range If clbit is outside the register.
   * @throws std::runtime_error If the bit has no stored outcome.
   */
  bool get(Clbit clbit) const;

  /**
   * @brief Returns an optional measured classical bit value.
   *
   * @param clbit Classical bit index.
   * @return Empty optional when the bit has not been measured.
   *
   * @throws std::out_of_range If clbit is outside the register.
   */
  std::optional<bool> get_raw(Clbit clbit) const;

  /**
   * @brief Stores a measurement outcome in a classical bit.
   *
   * @param clbit Classical bit index.
   * @param value Boolean measurement value.
   *
   * @throws std::out_of_range If clbit is outside the register.
   */
  void set(Clbit clbit, bool value);

  /**
   * @brief Clears all stored measurement outcomes.
   */
  void clear();

  /**
   * @brief Converts measured bits into a bit string.
   *
   * Unset bits are skipped. The string follows classical bit index order.
   *
   * @return Bit string containing measured outcomes.
   */
  std::string bit_string() const;

  /**
   * @brief Returns all classical bits in this register.
   *
   * @return Read-only bit storage.
   */
  const std::vector<ClassicalBit>& bits() const noexcept;

private:
  std::vector<ClassicalBit> bits_;

  void validate_clbit(Clbit clbit) const;
};

/**
 * @brief Indexed qubit handle.
 *
 * QuantumBit identifies a qubit by its index inside a quantum register or
 * circuit. It is a lightweight value object and does not own quantum state.
 */
class QuantumBit {
public:
  /**
   * @brief Creates a qubit handle.
   *
   * @param index Qubit index.
   */
  explicit QuantumBit(Qubit index);

  /**
   * @brief Returns the qubit index.
   */
  Qubit index() const noexcept;

  /**
   * @brief Returns the optional qubit label.
   */
  const std::optional<std::string>& label() const noexcept;

  /**
   * @brief Sets a human-readable qubit label.
   *
   * @param label Label assigned to this qubit.
   */
  void set_label(std::string label);

private:
  Qubit index_{};
  std::optional<std::string> label_;
};

/**
 * @brief Register containing indexed qubit handles.
 *
 * QuantumRegister owns qubit handles only. It does not own amplitudes,
 * density matrices, or simulator state.
 */
class QuantumRegister {
public:
  /**
   * @brief Creates a quantum register with indices [0, size).
   *
   * @param size Number of qubits in the register.
   */
  explicit QuantumRegister(Index size = 0);

  /**
   * @brief Returns the number of qubits.
   */
  Index size() const noexcept;

  /**
   * @brief Returns true when the register is empty.
   */
  bool empty() const noexcept;

  /**
   * @brief Returns the optional register label.
   */
  const std::optional<std::string>& label() const noexcept;

  /**
   * @brief Sets a human-readable register label.
   *
   * @param label Label assigned to this register.
   */
  void set_label(std::string label);

  /**
   * @brief Returns a qubit by index.
   *
   * @param qubit Qubit index.
   * @return Read-only qubit handle.
   *
   * @throws std::out_of_range If qubit is outside the register.
   */
  const QuantumBit& at(Qubit qubit) const;

  /**
   * @brief Returns a mutable qubit by index.
   *
   * @param qubit Qubit index.
   * @return Mutable qubit handle.
   *
   * @throws std::out_of_range If qubit is outside the register.
   */
  QuantumBit& at(Qubit qubit);

  /**
   * @brief Returns a qubit by index.
   *
   * @note This performs the same validation as at().
   */
  const QuantumBit& operator[](Qubit qubit) const;

  /**
   * @brief Returns a mutable qubit by index.
   *
   * @note This performs the same validation as at().
   */
  QuantumBit& operator[](Qubit qubit);

  /**
   * @brief Appends one qubit and returns its handle.
   *
   * The new qubit receives the next available index.
   */
  QuantumBit& add_qubit();

  /**
   * @brief Appends multiple qubits.
   *
   * @param count Number of qubits to append.
   */
  void add_qubits(Index count);

  /**
   * @brief Returns all qubit handles.
   */
  const std::vector<QuantumBit>& qubits() const noexcept;

private:
  std::vector<QuantumBit> qubits_;
  std::optional<std::string> label_;

  void validate_qubit(Qubit qubit) const;
};

}  // namespace quira
