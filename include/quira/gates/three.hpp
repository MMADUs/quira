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
 * Toffoli gate, also known as controlled-controlled-X.
 */
class Toffoli final : public ThreeQubit<Toffoli> {
public:
  Toffoli(Qubit first_control, Qubit second_control, Qubit target);

  std::string name() const override;

  Matrix unitary() const override;
  Qubit first_control() const noexcept;
  Qubit second_control() const noexcept;
  Qubit target() const noexcept;
};

/**
 * Fredkin gate, also known as controlled-SWAP.
 */
class Fredkin final : public ThreeQubit<Fredkin> {
public:
  Fredkin(Qubit control, Qubit first_target, Qubit second_target);

  std::string name() const override;

  Matrix unitary() const override;
  Qubit control() const noexcept;
  Qubit first_target() const noexcept;
  Qubit second_target() const noexcept;
};

}  // namespace quira
