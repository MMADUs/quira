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

#include <memory>
#include <string>
#include <vector>

namespace quira {

/**
 * @brief Base interface for gates that can be recorded in a quantum circuit.
 *
 * Concrete gate types own their parameters and qubit metadata. Simulators use
 * this interface to retrieve the gate name, unitary matrix, and target qubits
 * without depending on the concrete gate type.
 */
class QuantumGate {
public:
  QuantumGate() = default;
  virtual ~QuantumGate() = default;

  QuantumGate(const QuantumGate&) = default;
  QuantumGate& operator=(const QuantumGate&) = default;

  QuantumGate(QuantumGate&&) noexcept = default;
  QuantumGate& operator=(QuantumGate&&) noexcept = default;

  /**
   * @brief Returns the human-readable gate name.
   */
  virtual std::string name() const = 0;

  /**
   * @brief Returns the unitary matrix representation of the gate.
   */
  [[nodiscard]] virtual types::c_mat unitary() const = 0;

  /**
   * @brief Returns the qubits acted on by the gate.
   *
   * @return Target qubits in the order expected by the gate matrix.
   */
  [[nodiscard]] virtual std::vector<types::qubit> targets() const = 0;

  /**
   * @brief Returns a heap-allocated copy of the concrete gate.
   *
   * QuantumCircuit stores gates through QuantumGate pointers while preserving
   * value-like ownership. Cloning lets each concrete gate copy itself without
   * object slicing.
   */
  [[nodiscard]] virtual std::unique_ptr<QuantumGate> clone() const = 0;
};

/**
 * @brief CRTP helper for cloneable quantum gates.
 *
 * @tparam Derived Concrete gate type.
 *
 * This helper implements QuantumGate::clone() for concrete gate classes that
 * are copy-constructible.
 */
template<typename Derived>
class CloneableGate : public QuantumGate {
public:
  [[nodiscard]] std::unique_ptr<QuantumGate> clone() const override {
    return std::make_unique<Derived>(static_cast<const Derived&>(*this));
  }
};

/**
 * @brief Base class for cloneable single-qubit gates.
 *
 * @tparam Derived Concrete gate type.
 *
 * @note target() expose qubit acted on by the gate
 */
template<typename Derived>
class SingleQubit : public CloneableGate<Derived> {
public:
  explicit SingleQubit(types::qubit target) : target_(target) {}

  [[nodiscard]] std::vector<types::qubit> targets() const override { return {target_}; }

  types::qubit target() const noexcept { return target_; }

private:
  types::qubit target_{};
};

/**
 * @brief Base class for cloneable two-qubit gates.
 *
 * @tparam Derived Concrete gate type.
 *
 * @note Controlled gates may expose semantic accessors such as control() and
 * target() in their concrete classes.
 */
template<typename Derived>
class TwoQubit : public CloneableGate<Derived> {
public:
  TwoQubit(types::qubit first, types::qubit second) : first_(first), second_(second) {}

  [[nodiscard]] std::vector<types::qubit> targets() const override {
    return {first_, second_};
  }

  types::qubit first() const noexcept { return first_; }
  types::qubit second() const noexcept { return second_; }

private:
  types::qubit first_{};
  types::qubit second_{};
};

/**
 * @brief Base class for cloneable multi-gates acting on three or more qubits.
 *
 * @tparam Derived Concrete gate type.
 *
 * Stores an ordered target-qubit list and implements QuantumGate::targets().
 * Concrete gates can expose semantic accessors such as control(), target(),
 * first_control(), or second_target() when those names make the API clearer.
 */
template<typename Derived>
class MultiQubit : public CloneableGate<Derived> {
public:
  explicit MultiQubit(std::vector<types::qubit> targets)
      : targets_(std::move(targets)) {}

  [[nodiscard]] std::vector<types::qubit> targets() const override { return targets_; }

protected:
  [[nodiscard]] const std::vector<types::qubit>& target_list() const noexcept {
    return targets_;
  }

private:
  std::vector<types::qubit> targets_;
};

}  // namespace quira
