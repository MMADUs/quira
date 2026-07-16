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

#include <numbers>
#include <quira/types.hpp>

namespace quira::constants {

/**
 * @brief Represent the constant pi = 3.141592653589793238462643383279502884...
 */
inline constexpr types::real_n PI = std::numbers::pi;

/**
 * @brief Represent the sqrt of 2 = 1.414213562373095048801688724209698079...
 */
inline constexpr types::real_n SQRT_2 = std::numbers::sqrt2;

/**
 * @brief Represent the inverse sqrt of 2 = 0.707106781186547524400844362104849...
 */
inline constexpr types::real_n INV_SQRT_2 = 1.0 / SQRT_2;

/**
 * @brief Represent the constant epsilon 10^{-12} for thresholding
 */
inline constexpr types::real_n EPS = 1e-12;

/**
 * @brief Represent the complex number of 0 = 0.0 + 0.0i
 */
inline constexpr types::cplx_n ZERO = types::cplx_n{0.0, 0.0};

/**
 * @brief Represent the complex number of 1 = 1.0 + 0.0i
 */
inline constexpr types::cplx_n ONE = types::cplx_n{1.0, 0.0};

/**
 * @brief Represent the imaginary i = 0.0 + 1.0i or sqrt(-1)
 */
inline constexpr types::cplx_n IM = types::cplx_n{0.0, 1.0};

}  // namespace quira::constants
