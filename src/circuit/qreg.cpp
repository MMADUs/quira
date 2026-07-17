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

#include "quira/exception.hpp"
#include "quira/register.hpp"

#include <sstream>

namespace quira {

QuantumBit::QuantumBit(types::qubit index) : index_(index) {
}

types::qubit QuantumBit::index() const noexcept {
  return index_;
}

const std::optional<std::string>& QuantumBit::label() const noexcept {
  return label_;
}

void QuantumBit::set_label(std::string label) {
  label_ = std::move(label);
}

QuantumRegister::QuantumRegister(std::size_t size) {
  qubits_.reserve(size);

  for (std::size_t i = 0; i < size; ++i) {
    qubits_.emplace_back(static_cast<types::qubit>(i));
  }
}

std::size_t QuantumRegister::size() const noexcept {
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

const QuantumBit& QuantumRegister::at(types::qubit qubit) const {
  validate_qubit(qubit);
  return qubits_[qubit];
}

QuantumBit& QuantumRegister::at(types::qubit qubit) {
  validate_qubit(qubit);
  return qubits_[qubit];
}

const QuantumBit& QuantumRegister::operator[](types::qubit qubit) const {
  return at(qubit);
}

QuantumBit& QuantumRegister::operator[](types::qubit qubit) {
  return at(qubit);
}

QuantumBit& QuantumRegister::add_qubit() {
  const types::qubit index = qubits_.size();
  qubits_.emplace_back(index);
  return qubits_.back();
}

void QuantumRegister::add_qubits(std::size_t count) {
  qubits_.reserve(qubits_.size() + count);

  for (std::size_t i = 0; i < count; ++i) {
    add_qubit();
  }
}

const std::vector<QuantumBit>& QuantumRegister::qubits() const noexcept {
  return qubits_;
}

void QuantumRegister::validate_qubit(types::qubit qubit) const {
  if (qubit >= qubits_.size()) {
    std::ostringstream oss;
    oss << "qubit index " << qubit << " is out of range for register with "
        << qubits_.size() << " qubits";

    throw exception::OutOfRange("QuantumRegister::validate_qubit()", oss.str());
  }
}

}  // namespace quira
