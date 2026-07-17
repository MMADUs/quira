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
 * Controlled-NOT gate.
 */
class CNOT final : public TwoQubit<CNOT> {
public:
  CNOT(types::qubit control, types::qubit target);

  std::string name() const override;

  types::c_mat unitary() const override;
  types::qubit control() const noexcept;
  types::qubit target() const noexcept;
};

/**
 * Controlled-Z gate.
 */
class CZ final : public TwoQubit<CZ> {
public:
  CZ(types::qubit control, types::qubit target);

  std::string name() const override;

  types::c_mat unitary() const override;
  types::qubit control() const noexcept;
  types::qubit target() const noexcept;
};

/**
 * Swap gate.
 */
class Swap final : public TwoQubit<Swap> {
public:
  Swap(types::qubit control, types::qubit target);

  std::string name() const override;

  types::c_mat unitary() const override;
};

}  // namespace quira
