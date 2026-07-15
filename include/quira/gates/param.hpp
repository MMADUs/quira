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
 * Rotate X gate.
 */
class RotateX final : public SingleQubit<RotateX> {
public:
  RotateX(Qubit target, double theta);

  std::string name() const override;

  Matrix unitary() const override;
  double theta() const noexcept;

private:
  double theta_{};
};

/**
 * Rotate Y gate.
 */
class RotateY final : public SingleQubit<RotateY> {
public:
  RotateY(Qubit target, double theta);

  std::string name() const override;

  Matrix unitary() const override;
  double theta() const noexcept;

private:
  double theta_{};
};

/**
 * Rotate Z gate.
 */
class RotateZ final : public SingleQubit<RotateZ> {
public:
  RotateZ(Qubit target, double theta);

  std::string name() const override;

  Matrix unitary() const override;
  double theta() const noexcept;

private:
  double theta_{};
};

}  // namespace quira
