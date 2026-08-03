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

#include "quira/random.hpp"

#include "quira/constants.hpp"
#include "quira/exception.hpp"
#include "quira/linalg.hpp"

#include <cmath>
#include <istream>
#include <ostream>
#include <random>
#include <sstream>

namespace quira {

namespace {

std::mt19937_64 make_nondeterministic_engine() {
  std::random_device device;

  std::seed_seq seed{
      device(), device(), device(), device(), device(), device(), device(), device(),
  };

  return std::mt19937_64{seed};
}

void validate_finite(types::real_n value, const char* name, const char* caller) {
  if (!std::isfinite(value)) {
    std::ostringstream oss;
    oss << name << " must be finite";

    throw exception::InvalidArgument(caller, oss.str());
  }
}

}  // namespace

RandomEngine::RandomEngine() : engine_(make_nondeterministic_engine()) {
}

RandomEngine::RandomEngine(std::uint64_t seed_value) : engine_(seed_value) {
}

void RandomEngine::seed(std::uint64_t seed_value) {
  engine_.seed(seed_value);
}

RandomEngine::engine_type& RandomEngine::native() noexcept {
  return engine_;
}

const RandomEngine::engine_type& RandomEngine::native() const noexcept {
  return engine_;
}

std::ostream& RandomEngine::save(std::ostream& output) const {
  return output << engine_;
}

std::istream& RandomEngine::load(std::istream& input) {
  return input >> engine_;
}

types::real_n RandomEngine::uniform(types::real_n a, types::real_n b) {
  validate_finite(a, "a", "RandomEngine::uniform()");
  validate_finite(b, "b", "RandomEngine::uniform()");

  if (!(a < b)) {
    throw exception::OutOfRange("RandomEngine::uniform()", "requires a < b");
  }

  std::uniform_real_distribution<types::real_n> distribution(a, b);
  return distribution(engine_);
}

std::size_t RandomEngine::uniform_index(std::size_t a, std::size_t b) {
  if (a > b) {
    throw exception::OutOfRange("RandomEngine::uniform_index()", "requires a <= b");
  }

  std::uniform_int_distribution<std::size_t> distribution(a, b);
  return distribution(engine_);
}

types::real_n RandomEngine::normal(types::real_n mean, types::real_n sigma) {
  validate_finite(mean, "mean", "RandomEngine::normal()");
  validate_finite(sigma, "sigma", "RandomEngine::normal()");

  if (!(sigma > types::real_n{0.0})) {
    throw exception::OutOfRange("RandomEngine::normal()", "requires sigma > 0");
  }

  std::normal_distribution<types::real_n> distribution(mean, sigma);
  return distribution(engine_);
}

bool RandomEngine::bernoulli(types::real_n probability) {
  validate_finite(probability, "probability", "RandomEngine::bernoulli()");

  if (probability < types::real_n{0.0} || probability > types::real_n{1.0}) {
    throw exception::OutOfRange("RandomEngine::bernoulli()",
                                "requires probability in [0, 1]");
  }

  std::bernoulli_distribution distribution(probability);
  return distribution(engine_);
}

std::size_t RandomEngine::sample_index(
    const std::vector<types::real_n>& probabilities) {
  if (probabilities.empty()) {
    throw exception::InvalidArgument("RandomEngine::sample_index()",
                                     "probability vector cannot be empty");
  }

  types::real_n total = 0.0;
  std::size_t last_positive = probabilities.size() - 1;

  for (std::size_t index = 0; index < probabilities.size(); ++index) {
    const types::real_n value = probabilities[index];

    if (!std::isfinite(value) || value < types::real_n{0.0}) {
      std::ostringstream oss;
      oss << "probability at index " << index << " must be finite and non-negative";

      throw exception::InvalidArgument("RandomEngine::sample_index()", oss.str());
    }

    if (value > types::real_n{0.0}) {
      last_positive = index;
    }

    total += value;
  }

  if (!(total > types::real_n{0.0})) {
    throw exception::InvalidArgument("RandomEngine::sample_index()",
                                     "probabilities must sum to a positive value");
  }

  const types::real_n draw = uniform(types::real_n{0.0}, total);
  types::real_n cumulative = 0.0;

  for (std::size_t index = 0; index < probabilities.size(); ++index) {
    cumulative += probabilities[index];

    if (draw < cumulative) {
      return index;
    }
  }

  return last_positive;
}

types::ket random_state(std::size_t num_qubits, RandomEngine& rng) {
  const std::size_t dimension = linalg::basis_size(num_qubits);

  types::ket state(static_cast<Eigen::Index>(dimension));
  for (Eigen::Index i = 0; i < state.size(); ++i) {
    state(i) = types::cplx_n{rng.normal(), rng.normal()};
  }

  return linalg::normalize(state);
}

types::c_mat random_unitary(std::size_t dimension, RandomEngine& rng) {
  if (dimension == 0) {
    throw exception::InvalidArgument("random_unitary()",
                                     "dimension must be greater than zero");
  }

  types::c_mat matrix(static_cast<Eigen::Index>(dimension),
                      static_cast<Eigen::Index>(dimension));

  const types::real_n scale = types::real_n{1.0} / std::sqrt(types::real_n{2.0});
  for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
    for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
      matrix(row, col) = scale * types::cplx_n{rng.normal(), rng.normal()};
    }
  }

  Eigen::HouseholderQR<types::c_mat> qr(matrix);
  types::c_mat q =
      qr.householderQ() * types::c_mat::Identity(matrix.rows(), matrix.cols());
  types::c_mat r = qr.matrixQR().template triangularView<Eigen::Upper>();

  for (Eigen::Index i = 0; i < r.rows(); ++i) {
    const types::cplx_n diagonal = r(i, i);
    const types::real_n magnitude = std::abs(diagonal);

    if (magnitude > constants::EPS) {
      q.col(i) *= diagonal / magnitude;
    }
  }

  return q;
}

types::c_mat random_hermitian(std::size_t dimension, RandomEngine& rng) {
  if (dimension == 0) {
    throw exception::InvalidArgument("random_hermitian()",
                                     "dimension must be greater than zero");
  }

  types::c_mat matrix(static_cast<Eigen::Index>(dimension),
                      static_cast<Eigen::Index>(dimension));

  for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
    for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
      matrix(row, col) = types::cplx_n{rng.normal(), rng.normal()};
    }
  }

  return (matrix + matrix.adjoint()) / types::real_n{2.0};
}

std::vector<types::real_n> random_probability_vector(std::size_t size,
                                                     RandomEngine& rng) {
  if (size == 0) {
    throw exception::InvalidArgument("random_probability_vector()",
                                     "size must be greater than zero");
  }

  std::exponential_distribution<types::real_n> distribution(1.0);

  std::vector<types::real_n> probabilities(size);
  types::real_n total = 0.0;

  for (types::real_n& value : probabilities) {
    value = distribution(rng.native());
    total += value;
  }

  for (types::real_n& value : probabilities) {
    value /= total;
  }

  return probabilities;
}

std::vector<std::size_t> random_permutation(std::size_t size, RandomEngine& rng) {
  if (size == 0) {
    throw exception::InvalidArgument("random_permutation()",
                                     "size must be greater than zero");
  }

  std::vector<std::size_t> permutation(size);
  std::iota(permutation.begin(), permutation.end(), 0);
  std::shuffle(permutation.begin(), permutation.end(), rng.native());

  return permutation;
}

}  // namespace quira
