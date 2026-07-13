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

#include "quira/gate.hpp"

namespace quira {
/**
 * Identity gate.
 */
class Identity final : public SingleQubit<Identity> {
public:
  explicit Identity(Qubit target);

  std::string name() const override;
  Matrix unitary() const override;
};

/**
 * Hadamard gate.
 */
class Hadamard final : public SingleQubit<Hadamard> {
public:
  explicit Hadamard(Qubit target);

  std::string name() const override;

  Matrix unitary() const override;
};

/**
 * Pauli X gate.
 */
class PauliX final : public SingleQubit<PauliX> {
public:
  explicit PauliX(Qubit target);

  std::string name() const override;

  Matrix unitary() const override;
};

/**
 * Pauli Y gate.
 */
class PauliY final : public SingleQubit<PauliY> {
public:
  explicit PauliY(Qubit target);

  std::string name() const override;

  Matrix unitary() const override;
};

/**
 * Pauli Z gate.
 */
class PauliZ final : public SingleQubit<PauliZ> {
public:
  explicit PauliZ(Qubit target);

  std::string name() const override;

  Matrix unitary() const override;
};

/**
 * S phase gate: diag(1, i).
 */
class SGate final : public SingleQubit<SGate> {
public:
  explicit SGate(Qubit target);

  std::string name() const override;

  Matrix unitary() const override;
};

/**
 * T phase gate: diag(1, exp(i*pi/4)).
 */
class TGate final : public SingleQubit<TGate> {
public:
  explicit TGate(Qubit target);

  std::string name() const override;

  Matrix unitary() const override;
};

}  // namespace quira
