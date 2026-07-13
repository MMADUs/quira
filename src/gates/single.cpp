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

#include "quira/gates/single.hpp"

#include "quira/constants.hpp"

#include <cmath>
#include <complex>
#include <sstream>

namespace quira {
namespace {

const Complex I{0.0, 1.0};

std::string single_gate_name(const char* gate_name, Qubit target) {
  std::ostringstream oss;
  oss << gate_name << "(target=" << target << ")";
  return oss.str();
}

}  // namespace

/**
 * Identity.
 *
 * Leaves the target qubit unchanged.
 */
Identity::Identity(Qubit target) : SingleQubit<Identity>(target) {
}

std::string Identity::name() const {
  return single_gate_name("I", target_);
}

Matrix Identity::unitary() const {
  return Matrix::Identity(2, 2);
}

/**
 * Hadamard.
 *
 * Creates an equal superposition from computational basis states.
 */
Hadamard::Hadamard(Qubit target) : SingleQubit<Hadamard>(target) {
}

std::string Hadamard::name() const {
  return single_gate_name("H", target_);
}

Matrix Hadamard::unitary() const {
  Matrix matrix(2, 2);
  matrix << INV_SQRT_2, INV_SQRT_2, INV_SQRT_2, -INV_SQRT_2;
  return matrix;
}

/**
 * Pauli X.
 *
 * Bit-flip gate that maps |0> to |1> and |1> to |0>.
 */
PauliX::PauliX(Qubit target) : SingleQubit<PauliX>(target) {
}

std::string PauliX::name() const {
  return single_gate_name("X", target_);
}

Matrix PauliX::unitary() const {
  Matrix matrix(2, 2);
  matrix << 0.0, 1.0, 1.0, 0.0;
  return matrix;
}

/**
 * Pauli Y.
 *
 * Bit-and-phase-flip gate with imaginary off-diagonal entries.
 */
PauliY::PauliY(Qubit target) : SingleQubit<PauliY>(target) {
}

std::string PauliY::name() const {
  return single_gate_name("Y", target_);
}

Matrix PauliY::unitary() const {
  Matrix matrix(2, 2);
  matrix << Complex{0.0, 0.0}, -I, I, Complex{0.0, 0.0};
  return matrix;
}

/**
 * Pauli Z.
 *
 * Phase-flip gate that changes the sign of the |1> amplitude.
 */
PauliZ::PauliZ(Qubit target) : SingleQubit<PauliZ>(target) {
}

std::string PauliZ::name() const {
  return single_gate_name("Z", target_);
}

Matrix PauliZ::unitary() const {
  Matrix matrix(2, 2);
  matrix << 1.0, 0.0, 0.0, -1.0;
  return matrix;
}

/**
 * S phase.
 *
 * Applies a pi/2 phase shift to the |1> state.
 */
SGate::SGate(Qubit target) : SingleQubit<SGate>(target) {
}

std::string SGate::name() const {
  return single_gate_name("S", target_);
}

Matrix SGate::unitary() const {
  Matrix matrix(2, 2);
  matrix << Complex{1.0, 0.0}, Complex{0.0, 0.0}, Complex{0.0, 0.0}, I;
  return matrix;
}

/**
 * T phase.
 *
 * Applies a pi/4 phase shift to the |1> state.
 */
TGate::TGate(Qubit target) : SingleQubit<TGate>(target) {
}

std::string TGate::name() const {
  return single_gate_name("T", target_);
}

Matrix TGate::unitary() const {
  const Complex phase = std::exp(I * (PI / 4.0));

  Matrix matrix(2, 2);
  matrix << Complex{1.0, 0.0}, Complex{0.0, 0.0}, Complex{0.0, 0.0}, phase;
  return matrix;
}

}  // namespace quira
