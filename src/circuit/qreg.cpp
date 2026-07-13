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

#include "quira/register.hpp"

#include <sstream>
#include <stdexcept>

namespace quira {

/**
 * Classical bit.
 *
 * Holds one optional measurement outcome.
 */
QuantumBit::QuantumBit(Qubit index) : index_(index) {
}

Qubit QuantumBit::index() const noexcept {
  return index_;
}

const std::optional<std::string>& QuantumBit::label() const noexcept {
  return label_;
}

void QuantumBit::set_label(std::string label) {
  label_ = std::move(label);
}

/**
 * Classical register.
 *
 * Stores measurement outcomes in classical bit index order.
 */
QuantumRegister::QuantumRegister(Index size) {
  qubits_.reserve(size);
  for (Qubit qubit = 0; qubit < size; ++qubit) {
    qubits_.emplace_back(qubit);
  }
}

Index QuantumRegister::size() const noexcept {
  return qubits_.size();
}

bool QuantumRegister::empty() const noexcept {
  return qubits_.empty();
}

const std::optional<std::string>& QuantumRegister::label() const noexcept {
  return label_;
}

void QuantumRegister::set_label(std::string label) {
  label_ = std::move(label);
}

const QuantumBit& QuantumRegister::at(Qubit qubit) const {
  validate_qubit(qubit);
  return qubits_[qubit];
}

QuantumBit& QuantumRegister::at(Qubit qubit) {
  validate_qubit(qubit);
  return qubits_[qubit];
}

const QuantumBit& QuantumRegister::operator[](Qubit qubit) const {
  return at(qubit);
}

QuantumBit& QuantumRegister::operator[](Qubit qubit) {
  return at(qubit);
}

QuantumBit& QuantumRegister::add_qubit() {
  const Qubit index = qubits_.size();
  qubits_.emplace_back(index);
  return qubits_.back();
}

void QuantumRegister::add_qubits(Index count) {
  qubits_.reserve(qubits_.size() + count);
  for (Index i = 0; i < count; ++i) {
    add_qubit();
  }
}

const std::vector<QuantumBit>& QuantumRegister::qubits() const noexcept {
  return qubits_;
}

void QuantumRegister::validate_qubit(Qubit qubit) const {
  if (qubit >= qubits_.size()) {
    std::ostringstream oss;
    oss << "Qubit index " << qubit << " is out of range for register with "
        << qubits_.size() << " qubits";
    throw std::out_of_range(oss.str());
  }
}

}  // namespace quira
