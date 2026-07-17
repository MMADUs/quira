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

#include "quira/gates/param.hpp"

#include "quira/constants.hpp"

#include <cmath>
#include <sstream>

namespace quira {

namespace {

std::string rotation_gate_name(const char* gate_name, types::qubit target,
                               types::real_n theta) {
  std::ostringstream oss;
  oss << gate_name << "(target=" << target << ", theta=" << theta << ")";
  return oss.str();
}

}  // namespace

RX::RX(types::qubit target, types::real_n theta)
    : SingleQubit<RX>(target), theta_(theta) {
}

std::string RX::name() const {
  return rotation_gate_name("RX", target(), theta_);
}

types::c_mat RX::unitary() const {
  const types::real_n c = std::cos(theta_ / 2.0);
  const types::real_n s = std::sin(theta_ / 2.0);

  types::c_mat matrix(2, 2);
  matrix << types::cplx_n{c, 0.0}, types::cplx_n{0.0, -s}, types::cplx_n{0.0, -s},
      types::cplx_n{c, 0.0};
  return matrix;
}

types::real_n RX::theta() const noexcept {
  return theta_;
}

RY::RY(types::qubit target, types::real_n theta)
    : SingleQubit<RY>(target), theta_(theta) {
}

std::string RY::name() const {
  return rotation_gate_name("RY", target(), theta_);
}

types::c_mat RY::unitary() const {
  const types::real_n c = std::cos(theta_ / 2.0);
  const types::real_n s = std::sin(theta_ / 2.0);

  types::c_mat matrix(2, 2);
  matrix << types::cplx_n{c, 0.0}, types::cplx_n{-s, 0.0}, types::cplx_n{s, 0.0},
      types::cplx_n{c, 0.0};
  return matrix;
}

types::real_n RY::theta() const noexcept {
  return theta_;
}

RZ::RZ(types::qubit target, types::real_n theta)
    : SingleQubit<RZ>(target), theta_(theta) {
}

std::string RZ::name() const {
  return rotation_gate_name("RZ", target(), theta_);
}

types::c_mat RZ::unitary() const {
  const types::cplx_n negative_phase = std::exp(-constants::IM * (theta_ / 2.0));
  const types::cplx_n positive_phase = std::exp(constants::IM * (theta_ / 2.0));

  types::c_mat matrix(2, 2);
  matrix << negative_phase, types::cplx_n{0.0, 0.0}, types::cplx_n{0.0, 0.0},
      positive_phase;
  return matrix;
}

types::real_n RZ::theta() const noexcept {
  return theta_;
}

}  // namespace quira
