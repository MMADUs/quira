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

#include <cmath>
#include <complex>
#include <sstream>

namespace quira {
namespace {

const Complex I{0.0, 1.0};

std::string rotation_gate_name(const char* gate_name, Qubit target, double theta) {
  std::ostringstream oss;
  oss << gate_name << "(target=" << target << ", theta=" << theta << ")";
  return oss.str();
}

}  // namespace

RotateX::RotateX(Qubit target, double theta)
    : SingleQubit<RotateX>(target), theta_(theta) {
}

std::string RotateX::name() const {
  return rotation_gate_name("RX", target_, theta_);
}

Matrix RotateX::unitary() const {
  const double c = std::cos(theta_ / 2.0);
  const double s = std::sin(theta_ / 2.0);

  Matrix matrix(2, 2);
  matrix << Complex{c, 0.0}, Complex{0.0, -s}, Complex{0.0, -s}, Complex{c, 0.0};
  return matrix;
}

double RotateX::theta() const noexcept {
  return theta_;
}

RotateY::RotateY(Qubit target, double theta)
    : SingleQubit<RotateY>(target), theta_(theta) {
}

std::string RotateY::name() const {
  return rotation_gate_name("RY", target_, theta_);
}

Matrix RotateY::unitary() const {
  const double c = std::cos(theta_ / 2.0);
  const double s = std::sin(theta_ / 2.0);

  Matrix matrix(2, 2);
  matrix << Complex{c, 0.0}, Complex{-s, 0.0}, Complex{s, 0.0}, Complex{c, 0.0};
  return matrix;
}

double RotateY::theta() const noexcept {
  return theta_;
}

RotateZ::RotateZ(Qubit target, double theta)
    : SingleQubit<RotateZ>(target), theta_(theta) {
}

std::string RotateZ::name() const {
  return rotation_gate_name("RZ", target_, theta_);
}

Matrix RotateZ::unitary() const {
  const Complex negative_phase = std::exp(-I * (theta_ / 2.0));
  const Complex positive_phase = std::exp(I * (theta_ / 2.0));

  Matrix matrix(2, 2);
  matrix << negative_phase, Complex{0.0, 0.0}, Complex{0.0, 0.0}, positive_phase;
  return matrix;
}

double RotateZ::theta() const noexcept {
  return theta_;
}

}  // namespace quira
