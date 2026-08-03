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

#include "quira/types.hpp"

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <random>
#include <vector>

namespace quira {

/**
 * @brief Central random-number engine used by Quira components.
 *
 * RandomEngine owns the pseudo-random generator used by simulators, samplers,
 * measurements, resets, and future randomized state-vector utilities. Keeping
 * randomness explicit makes simulations reproducible and avoids hidden global
 * mutable state.
 */
class RandomEngine {
public:
  using engine_type = std::mt19937_64;

  /**
   * @brief Creates an engine with non-deterministic seeding.
   */
  RandomEngine();

  /**
   * @brief Creates an engine with deterministic seeding.
   *
   * @param seed Seed value for reproducible random draws.
   */
  explicit RandomEngine(std::uint64_t seed);

  /**
   * @brief Reseeds the engine.
   *
   * @param seed Seed value for reproducible random draws.
   */
  void seed(std::uint64_t seed);

  /**
   * @brief Returns the underlying standard-library engine.
   *
   * This is intended for internal algorithms that need direct interoperability
   * with standard algorithms such as std::shuffle.
   */
  [[nodiscard]] engine_type& native() noexcept;

  /**
   * @brief Returns the underlying standard-library engine.
   */
  [[nodiscard]] const engine_type& native() const noexcept;

  /**
   * @brief Saves the pseudo-random generator state.
   *
   * @param output Stream receiving the engine state.
   * @return Reference to output.
   */
  std::ostream& save(std::ostream& output) const;

  /**
   * @brief Loads the pseudo-random generator state.
   *
   * @param input Stream containing a previously saved engine state.
   * @return Reference to input.
   */
  std::istream& load(std::istream& input);

  /**
   * @brief Draws a uniformly distributed real number in [a, b).
   *
   * @throws exception::OutOfRange If a >= b.
   */
  [[nodiscard]] types::real_n uniform(types::real_n a = 0.0, types::real_n b = 1.0);

  /**
   * @brief Draws a uniformly distributed index in [a, b].
   *
   * @throws exception::OutOfRange If a > b.
   */
  [[nodiscard]] std::size_t uniform_index(std::size_t a, std::size_t b);

  /**
   * @brief Draws a normally distributed real number.
   *
   * @throws exception::OutOfRange If sigma <= 0.
   */
  [[nodiscard]] types::real_n normal(types::real_n mean = 0.0,
                                     types::real_n sigma = 1.0);

  /**
   * @brief Draws a Bernoulli random boolean.
   *
   * @param probability Probability of returning true.
   *
   * @throws exception::OutOfRange If probability is outside [0, 1].
   */
  [[nodiscard]] bool bernoulli(types::real_n probability = 0.5);

  /**
   * @brief Samples one index from non-negative weights or probabilities.
   *
   * The values do not need to be normalized. This method is the central helper
   * used by state-vector measurement sampling.
   *
   * @param probabilities Non-empty non-negative weights.
   * @return Sampled index.
   *
   * @throws exception::InvalidArgument If probabilities is empty, contains a
   * negative/NaN value, or sums to zero.
   */
  [[nodiscard]] std::size_t sample_index(
      const std::vector<types::real_n>& probabilities);

  /**
   * @brief Generates a normalized random state vector for num_qubits qubits.
   */
  [[nodiscard]] types::ket random_state(std::size_t num_qubits, RandomEngine& rng);

  /**
   * @brief Generates a random unitary matrix.
   */
  [[nodiscard]] types::c_mat random_unitary(std::size_t dimension, RandomEngine& rng);

  /**
   * @brief Generates a random Hermitian matrix.
   */
  [[nodiscard]] types::c_mat random_hermitian(std::size_t dimension, RandomEngine& rng);

  /**
   * @brief Generates a normalized random probability vector.
   */
  [[nodiscard]] std::vector<types::real_n> random_probability_vector(std::size_t size,
                                                                     RandomEngine& rng);

  /**
   * @brief Generates a random permutation of [0, size).
   */
  [[nodiscard]] std::vector<std::size_t> random_permutation(std::size_t size,
                                                            RandomEngine& rng);

private:
  engine_type engine_;
};

}  // namespace quira
