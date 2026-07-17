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

#include "quira/unitary.hpp"

#include "quira/exception.hpp"
#include "quira/linalg.hpp"

#include <sstream>
#include <unordered_set>
#include <utility>

namespace quira {

namespace {

void validate_targets(const std::vector<types::qubit>& targets) {
  if (targets.empty()) {
    throw exception::InvalidArgument("unitary requires at least one target");
  }

  std::unordered_set<types::qubit> seen;

  for (types::qubit target : targets) {
    if (!seen.insert(target).second) {
      throw exception::InvalidQubit("unitary target qubits contains duplicates");
    }
  }
}

}  // namespace

Unitary::Unitary(std::string name, types::c_mat unitary,
                 std::vector<types::qubit> targets)
    : name_(std::move(name)),
      unitary_(std::move(unitary)),
      targets_(std::move(targets)) {
  validate_targets(targets_);

  const std::size_t dimension = linalg::basis_size(targets_.size());

  if (unitary_.rows() != static_cast<Eigen::Index>(dimension) ||
      unitary_.cols() != static_cast<Eigen::Index>(dimension)) {
    std::ostringstream oss;
    oss << "unitary matrix dimension does not match target qubits, got shape ("
        << unitary_.rows() << ", " << unitary_.cols() << ") for expected basis size "
        << dimension;

    throw exception::DimensionMismatch("Unitary::Unitary()", oss.str());
  }

  if (!linalg::is_unitary(unitary_)) {
    throw exception::InvalidArgument("Unitary::Unitary()",
                                     "unitary gate requires a unitary matrix");
  }
}

std::string Unitary::name() const {
  return name_;
}

types::c_mat Unitary::unitary() const {
  return unitary_;
}

std::vector<types::qubit> Unitary::targets() const {
  return targets_;
}

}  // namespace quira
