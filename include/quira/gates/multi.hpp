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
 * @brief Toffoli gate, also known as controlled-controlled-X.
 */
class Toff final : public MultiQubit<Toff> {
public:
  Toff(types::qubit first_control, types::qubit second_control, types::qubit target);

  std::string name() const override;

  types::c_mat unitary() const override;

  types::qubit first_control() const noexcept;
  types::qubit second_control() const noexcept;
  types::qubit target() const noexcept;
};

/**
 * @brief Fredkin gate, also known as controlled-SWAP.
 */
class Fred final : public MultiQubit<Fred> {
public:
  Fred(types::qubit control, types::qubit first_target, types::qubit second_target);

  std::string name() const override;

  types::c_mat unitary() const override;

  types::qubit control() const noexcept;
  types::qubit first_target() const noexcept;
  types::qubit second_target() const noexcept;
};

}  // namespace quira
