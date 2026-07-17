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

ClassicalBit::ClassicalBit(types::clbit index) : index_(index) {
}

types::clbit ClassicalBit::index() const noexcept {
  return index_;
}

bool ClassicalBit::has_value() const noexcept {
  return value_.has_value();
}

bool ClassicalBit::value() const {
  if (!value_) {
    std::ostringstream oss;
    oss << "classical bit " << index_ << " has no measurement outcome";
    throw exception::InvalidClbit("ClassicalBit::value()", oss.str());
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

ClassicalRegister::ClassicalRegister(std::size_t size) {
  bits_.reserve(size);

  for (std::size_t i = 0; i < size; ++i) {
    bits_.emplace_back(static_cast<types::clbit>(i));
  }
}

std::size_t ClassicalRegister::size() const noexcept {
  return bits_.size();
}

bool ClassicalRegister::has_value(types::clbit clbit) const {
  validate_clbit(clbit);
  return bits_[clbit].has_value();
}

bool ClassicalRegister::get(types::clbit clbit) const {
  validate_clbit(clbit);
  return bits_[clbit].value();
}

std::optional<bool> ClassicalRegister::get_raw(types::clbit clbit) const {
  validate_clbit(clbit);
  return bits_[clbit].raw_value();
}

void ClassicalRegister::set(types::clbit clbit, bool value) {
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
  for (std::size_t index = bits_.size(); index > 0; --index) {
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

void ClassicalRegister::validate_clbit(types::clbit clbit) const {
  if (clbit >= bits_.size()) {
    std::ostringstream oss;
    oss << "classical bit index " << clbit << " is out of range for register with "
        << bits_.size() << " bits";

    throw exception::OutOfRange("ClassicalRegister::validate_clbit()", oss.str());
  }
}

}  // namespace quira
