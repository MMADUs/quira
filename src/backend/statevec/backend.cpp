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

#include "quira/backend/statevec.hpp"

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace quira {

namespace {

Index basis_size(Index num_qubits) {
  if (num_qubits >= std::numeric_limits<Index>::digits) {
    throw std::out_of_range("Number of qubits is too large for Index");
  }

  // A pure n-qubit state has 2^n computational-basis amplitudes.
  return Index{1} << num_qubits;
}

bool is_power_of_two(Index value) {
  return value != 0 && (value & (value - 1)) == 0;
}

Index log2_exact(Index value) {
  Index result = 0;
  while (value > 1) {
    value >>= 1;
    ++result;
  }
  return result;
}

}  // namespace

/**
 * State vector.
 *
 * Initializes to |0...0> for the requested qubit count.
 */
StateVector::StateVector(Index num_qubits) {
  reset(num_qubits);
}

StateVector::StateVector(Vector amplitudes) {
  set_state(std::move(amplitudes));
}

Index StateVector::num_qubits() const noexcept {
  return num_qubits_;
}

Index StateVector::dimension() const noexcept {
  return static_cast<Index>(state_.size());
}

const Vector& StateVector::state() const noexcept {
  return state_;
}

void StateVector::reset(Index num_qubits) {
  num_qubits_ = num_qubits;
  state_ = Vector::Zero(static_cast<Eigen::Index>(basis_size(num_qubits)));

  // |0...0> is represented by basis index 0 with amplitude 1.
  state_[0] = Complex{1.0, 0.0};
}

void StateVector::set_state(Vector amplitudes) {
  const Index size = static_cast<Index>(amplitudes.size());
  if (!is_power_of_two(size)) {
    throw std::invalid_argument("State vector dimension must be a power of two");
  }

  state_ = std::move(amplitudes);
  num_qubits_ = log2_exact(size);
  normalize();
}

double StateVector::probability(Index basis_state) const {
  if (basis_state >= dimension()) {
    throw std::out_of_range("Basis state index is out of range");
  }

  // Born rule: probability of a basis state is |amplitude|^2.
  return std::norm(state_[static_cast<Eigen::Index>(basis_state)]);
}

double StateVector::norm() const {
  double total = 0.0;
  for (Eigen::Index i = 0; i < state_.size(); ++i) {
    total += std::norm(state_[i]);
  }
  return std::sqrt(total);
}

void StateVector::normalize() {
  const double current_norm = norm();
  if (current_norm == 0.0) {
    throw std::runtime_error("Cannot normalize a zero state vector");
  }
  state_ /= current_norm;
}

void StateVector::apply(const Matrix& unitary, const std::vector<Qubit>& targets) {
  validate_targets(targets);

  const Index target_count = targets.size();
  const Index local_dimension = basis_size(target_count);

  if (unitary.rows() != static_cast<Eigen::Index>(local_dimension) ||
      unitary.cols() != static_cast<Eigen::Index>(local_dimension)) {
    throw std::invalid_argument("Gate unitary dimension does not match targets");
  }

  // Qiskit-style little-endian convention:
  //
  //   global basis bit q stores qubit q
  //   qubit 0 is the least-significant bit of the basis index
  //
  // Example for 3 qubits:
  //
  //   basis index 5 = binary 101
  //   q0 = 1, q1 = 0, q2 = 1
  //
  // When written as a ket in the usual display order, this is |q2 q1 q0>.
  Vector next = Vector::Zero(state_.size());

  for (Index basis = 0; basis < dimension(); ++basis) {
    Index local_col = 0;

    // Extract the target-qubit bits from the global basis index.
    //
    // The local gate matrix also uses little endian:
    //
    //   targets[0] -> local bit 0
    //   targets[1] -> local bit 1
    //   targets[k] -> local bit k
    //
    // local_col is the input column of the small gate unitary.
    for (Index target_index = 0; target_index < target_count; ++target_index) {
      if ((basis >> targets[target_index]) & Index{1}) {
        local_col |= Index{1} << target_index;
      }
    }

    for (Index local_row = 0; local_row < local_dimension; ++local_row) {
      Index destination = basis;

      // Replace only the target-qubit bits in the global basis index.
      //
      // Non-target qubits stay unchanged. The chosen local_row supplies the
      // output bits for the target qubits after applying the gate.
      for (Index target_index = 0; target_index < target_count; ++target_index) {
        const Index mask = Index{1} << targets[target_index];
        if ((local_row >> target_index) & Index{1}) {
          destination |= mask;
        } else {
          destination &= ~mask;
        }
      }

      // Matrix-vector multiplication over the affected subspace:
      //
      //   next[destination] += U[local_row, local_col] * state[basis]
      //
      // This avoids building a full 2^n by 2^n matrix for every gate.
      next[static_cast<Eigen::Index>(destination)] +=
          unitary(static_cast<Eigen::Index>(local_row),
                  static_cast<Eigen::Index>(local_col)) *
          state_[static_cast<Eigen::Index>(basis)];
    }
  }

  state_ = std::move(next);
  normalize();
}

bool StateVector::measure(Qubit qubit, std::mt19937_64& rng) {
  validate_qubit(qubit);

  double probability_one = 0.0;

  // Sum the probabilities of every basis state whose measured qubit is 1.
  for (Index basis = 0; basis < dimension(); ++basis) {
    if ((basis >> qubit) & Index{1}) {
      probability_one += probability(basis);
    }
  }

  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  const bool outcome = distribution(rng) < probability_one;

  // Collapse the state by deleting amplitudes inconsistent with the sampled
  // outcome, then renormalize the remaining amplitudes.
  for (Index basis = 0; basis < dimension(); ++basis) {
    const bool bit_is_one = ((basis >> qubit) & Index{1}) != 0;
    if (bit_is_one != outcome) {
      state_[static_cast<Eigen::Index>(basis)] = Complex{0.0, 0.0};
    }
  }

  normalize();
  return outcome;
}

std::vector<bool> StateVector::measure_all(std::mt19937_64& rng) {
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  const double sample = distribution(rng);

  double cumulative = 0.0;
  Index measured_basis = dimension() - 1;

  // Sample one full computational-basis state from the state distribution.
  for (Index basis = 0; basis < dimension(); ++basis) {
    cumulative += probability(basis);
    if (sample <= cumulative) {
      measured_basis = basis;
      break;
    }
  }

  // Full-register measurement collapses the state to the sampled basis vector.
  state_ = Vector::Zero(state_.size());
  state_[static_cast<Eigen::Index>(measured_basis)] = Complex{1.0, 0.0};

  std::vector<bool> outcomes;
  outcomes.reserve(num_qubits_);

  // Return outcomes in qubit-index order: q0, q1, q2, ...
  for (Qubit qubit = 0; qubit < num_qubits_; ++qubit) {
    outcomes.push_back(((measured_basis >> qubit) & Index{1}) != 0);
  }

  return outcomes;
}

void StateVector::reset_qubit(Qubit qubit, std::mt19937_64& rng) {
  // Quantum reset is modeled as measure followed by X when the outcome is 1.
  const bool outcome = measure(qubit, rng);
  if (!outcome) {
    return;
  }

  Matrix x(2, 2);
  x << 0.0, 1.0, 1.0, 0.0;
  apply(x, {qubit});
}

void StateVector::validate_qubit(Qubit qubit) const {
  if (qubit >= num_qubits_) {
    std::ostringstream oss;
    oss << "Qubit index " << qubit << " is out of range for " << num_qubits_
        << " qubits";
    throw std::out_of_range(oss.str());
  }
}

void StateVector::validate_targets(const std::vector<Qubit>& targets) const {
  if (targets.empty()) {
    throw std::invalid_argument("StateVector::apply requires at least one target");
  }

  std::unordered_set<Qubit> seen;
  for (Qubit target : targets) {
    validate_qubit(target);
    if (!seen.insert(target).second) {
      throw std::invalid_argument("Gate target list contains duplicate qubits");
    }
  }
}

}  // namespace quira
