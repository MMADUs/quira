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
  [[nodiscard]] virtual Matrix unitary() const = 0;

  /**
   * @brief Returns the qubits acted on by the gate.
   *
   * @return Target qubits in the order expected by the gate matrix.
   */
  [[nodiscard]] virtual std::vector<Qubit> targets() const = 0;

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
  explicit SingleQubit(Qubit target) : target_(target) {}

  [[nodiscard]] std::vector<Qubit> targets() const override { return {target_}; }

  Qubit target() const noexcept { return target_; }

private:
  Qubit target_{};
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
  TwoQubit(Qubit first, Qubit second) : first_(first), second_(second) {}

  [[nodiscard]] std::vector<Qubit> targets() const override {
    return {first_, second_};
  }

  Qubit first() const noexcept { return first_; }
  Qubit second() const noexcept { return second_; }

private:
  Qubit first_{};
  Qubit second_{};
};

/**
 * @brief Base class for cloneable three-qubit gates.
 *
 * @tparam Derived Concrete gate type.
 *
 * @note Concrete gates may expose semantic accessors such as first_control(),
 * second_control(), control(), target(), first_target(), or second_target().
 */
template<typename Derived>
class ThreeQubit : public CloneableGate<Derived> {
public:
  ThreeQubit(Qubit first, Qubit second, Qubit third)
      : first_(first), second_(second), third_(third) {}

  [[nodiscard]] std::vector<Qubit> targets() const override {
    return {first_, second_, third_};
  }

  Qubit first() const noexcept { return first_; }
  Qubit second() const noexcept { return second_; }
  Qubit third() const noexcept { return third_; }

private:
  Qubit first_{};
  Qubit second_{};
  Qubit third_{};
};

}  // namespace quira
