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

#include "quira/gates/single.hpp"

#include "quira/constants.hpp"

#include <cmath>
#include <sstream>

namespace quira {

namespace {

std::string single_gate_name(const char* gate_name, types::qubit target) {
  std::ostringstream oss;
  oss << gate_name << "(target=" << target << ")";
  return oss.str();
}

}  // namespace

I::I(types::qubit target) : SingleQubit<I>(target) {
}

std::string I::name() const {
  return single_gate_name("I", target());
}

types::c_mat I::unitary() const {
  return types::c_mat::Identity(2, 2);
}

H::H(types::qubit target) : SingleQubit<H>(target) {
}

std::string H::name() const {
  return single_gate_name("H", target());
}

types::c_mat H::unitary() const {
  types::c_mat matrix(2, 2);
  matrix << constants::INV_SQRT_2, constants::INV_SQRT_2, constants::INV_SQRT_2,
      -constants::INV_SQRT_2;
  return matrix;
}

X::X(types::qubit target) : SingleQubit<X>(target) {
}

std::string X::name() const {
  return single_gate_name("X", target());
}

types::c_mat X::unitary() const {
  types::c_mat matrix(2, 2);
  matrix << constants::ZERO, constants::ONE, constants::ONE, constants::ZERO;
  return matrix;
}

Y::Y(types::qubit target) : SingleQubit<Y>(target) {
}

std::string Y::name() const {
  return single_gate_name("Y", target());
}

types::c_mat Y::unitary() const {
  types::c_mat matrix(2, 2);
  matrix << constants::ZERO, -constants::IM, constants::IM, constants::ZERO;
  return matrix;
}

Z::Z(types::qubit target) : SingleQubit<Z>(target) {
}

std::string Z::name() const {
  return single_gate_name("Z", target());
}

types::c_mat Z::unitary() const {
  types::c_mat matrix(2, 2);
  matrix << constants::ONE, constants::ZERO, constants::ZERO, -constants::ONE;
  return matrix;
}

S::S(types::qubit target) : SingleQubit<S>(target) {
}

std::string S::name() const {
  return single_gate_name("S", target());
}

types::c_mat S::unitary() const {
  types::c_mat matrix(2, 2);
  matrix << constants::ONE, constants::ZERO, constants::ZERO, constants::IM;
  return matrix;
}

T::T(types::qubit target) : SingleQubit<T>(target) {
}

std::string T::name() const {
  return single_gate_name("T", target());
}

types::c_mat T::unitary() const {
  const types::cplx_n phase =
      std::exp(constants::IM * (constants::PI / types::real_n{4}));

  types::c_mat matrix(2, 2);
  matrix << constants::ONE, constants::ZERO, constants::ZERO, phase;
  return matrix;
}

}  // namespace quira
