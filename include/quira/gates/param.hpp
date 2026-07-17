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
#include "quira/types.hpp"

namespace quira {
/**
 * Rotate X gate.
 */
class RX final : public SingleQubit<RX> {
public:
  RX(types::qubit target, types::real_n theta);

  std::string name() const override;

  types::c_mat unitary() const override;
  types::real_n theta() const noexcept;

private:
  types::real_n theta_{};
};

/**
 * Rotate Y gate.
 */
class RY final : public SingleQubit<RY> {
public:
  RY(types::qubit target, types::real_n theta);

  std::string name() const override;

  types::c_mat unitary() const override;
  types::real_n theta() const noexcept;

private:
  types::real_n theta_{};
};

/**
 * Rotate Z gate.
 */
class RZ final : public SingleQubit<RZ> {
public:
  RZ(types::qubit target, types::real_n theta);

  std::string name() const override;

  types::c_mat unitary() const override;
  types::real_n theta() const noexcept;

private:
  types::real_n theta_{};
};

}  // namespace quira
