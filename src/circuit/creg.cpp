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

#include "quira/register.hpp"

#include <sstream>
#include <stdexcept>

namespace quira {

/**
 * Classical bit.
 *
 * Holds one optional measurement outcome.
 */
ClassicalBit::ClassicalBit(Clbit index) : index_(index) {
}

Clbit ClassicalBit::index() const noexcept {
  return index_;
}

bool ClassicalBit::has_value() const noexcept {
  return value_.has_value();
}

bool ClassicalBit::value() const {
  if (!value_) {
    std::ostringstream oss;
    oss << "Classical bit " << index_ << " has no measurement outcome";
    throw std::runtime_error(oss.str());
  }
  return *value_;
}

std::optional<bool> ClassicalBit::raw_value() const noexcept {
  return value_;
}

void ClassicalBit::set(bool value) {
  value_ = value;
}

void ClassicalBit::clear() noexcept {
  value_.reset();
}

/**
 * Classical register.
 *
 * Stores measurement outcomes in classical bit index order.
 */
ClassicalRegister::ClassicalRegister(Index size) {
  bits_.reserve(size);
  for (Clbit bit = 0; bit < size; ++bit) {
    bits_.emplace_back(bit);
  }
}

Index ClassicalRegister::size() const noexcept {
  return bits_.size();
}

bool ClassicalRegister::has_value(Clbit clbit) const {
  validate_clbit(clbit);
  return bits_[clbit].has_value();
}

bool ClassicalRegister::get(Clbit clbit) const {
  validate_clbit(clbit);
  return bits_[clbit].value();
}

std::optional<bool> ClassicalRegister::get_raw(Clbit clbit) const {
  validate_clbit(clbit);
  return bits_[clbit].raw_value();
}

void ClassicalRegister::set(Clbit clbit, bool value) {
  validate_clbit(clbit);
  bits_[clbit].set(value);
}

void ClassicalRegister::clear() {
  for (ClassicalBit& bit : bits_) {
    bit.clear();
  }
}

std::string ClassicalRegister::bit_string() const {
  std::string result;
  result.reserve(bits_.size());

  // Match Qiskit display order: highest classical bit on the left, c0 on the
  // right. Unset bits are still skipped for the current output behavior.
  for (Index index = bits_.size(); index > 0; --index) {
    const ClassicalBit& bit = bits_[index - 1];
    if (bit.has_value()) {
      result.push_back(bit.value() ? '1' : '0');
    }
  }

  return result;
}

const std::vector<ClassicalBit>& ClassicalRegister::bits() const noexcept {
  return bits_;
}

void ClassicalRegister::validate_clbit(Clbit clbit) const {
  if (clbit >= bits_.size()) {
    std::ostringstream oss;
    oss << "Classical bit index " << clbit << " is out of range for register with "
        << bits_.size() << " bits";
    throw std::out_of_range(oss.str());
  }
}

}  // namespace quira
