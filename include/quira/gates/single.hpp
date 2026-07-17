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

#pragma once

#include "quira/gate.hpp"

namespace quira {
/**
 * Identity gate.
 */
class I final : public SingleQubit<I> {
public:
  explicit I(types::qubit target);

  std::string name() const override;
  types::c_mat unitary() const override;
};

/**
 * Hadamard gate.
 */
class H final : public SingleQubit<H> {
public:
  explicit H(types::qubit target);

  std::string name() const override;

  types::c_mat unitary() const override;
};

/**
 * Pauli X gate.
 */
class X final : public SingleQubit<X> {
public:
  explicit X(types::qubit target);

  std::string name() const override;

  types::c_mat unitary() const override;
};

/**
 * Pauli Y gate.
 */
class Y final : public SingleQubit<Y> {
public:
  explicit Y(types::qubit target);

  std::string name() const override;

  types::c_mat unitary() const override;
};

/**
 * Pauli Z gate.
 */
class Z final : public SingleQubit<Z> {
public:
  explicit Z(types::qubit target);

  std::string name() const override;

  types::c_mat unitary() const override;
};

/**
 * S phase gate: diag(1, i).
 */
class S final : public SingleQubit<S> {
public:
  explicit S(types::qubit target);

  std::string name() const override;

  types::c_mat unitary() const override;
};

/**
 * T phase gate: diag(1, exp(i*pi/4)).
 */
class T final : public SingleQubit<T> {
public:
  explicit T(types::qubit target);

  std::string name() const override;

  types::c_mat unitary() const override;
};

}  // namespace quira
