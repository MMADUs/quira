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

#include "quira/states.hpp"

#include "quira/constants.hpp"
#include "quira/exception.hpp"
#include "quira/linalg.hpp"

#include <cmath>

namespace quira::states {

namespace {

void validate_non_zero_qubit(std::size_t num_qubits, const char* caller) {
  if (num_qubits == 0) {
    throw exception::InvalidQubit(caller, "requires at least one qubit");
  }
}

bool has_odd_basis_parity(std::size_t basis_index) noexcept {
  bool odd = false;

  while (basis_index != 0) {
    odd = !odd;
    basis_index &= basis_index - 1;
  }

  return odd;
}

types::c_mat projector_from_state(const types::ket& state) {
  return linalg::projector(state);
}

}  // namespace

types::ket zero_state(std::size_t num_qubits) {
  const std::size_t dimension = linalg::basis_size(num_qubits);

  types::ket state = types::ket::Zero(static_cast<Eigen::Index>(dimension));
  state[0] = constants::ONE;

  return state;
}

types::ket one_state(std::size_t num_qubits) {
  const std::size_t dimension = linalg::basis_size(num_qubits);

  types::ket state = types::ket::Zero(static_cast<Eigen::Index>(dimension));
  state[static_cast<Eigen::Index>(dimension - 1)] = constants::ONE;

  return state;
}

types::ket basis_state(std::size_t num_qubits, std::size_t basis_index) {
  const std::size_t dimension = linalg::basis_size(num_qubits);

  if (basis_index >= dimension) {
    throw exception::OutOfRange("states::basis_state()",
                                "basis index is outside the state dimension");
  }

  types::ket state = types::ket::Zero(static_cast<Eigen::Index>(dimension));
  state[static_cast<Eigen::Index>(basis_index)] = constants::ONE;

  return state;
}

types::ket plus_state(std::size_t num_qubits) {
  const std::size_t dimension = linalg::basis_size(num_qubits);
  const types::real_n amplitude =
      1.0 / std::sqrt(static_cast<types::real_n>(dimension));

  types::ket state(static_cast<Eigen::Index>(dimension));
  state.setConstant(types::cplx_n{amplitude, 0.0});

  return state;
}

types::ket minus_state(std::size_t num_qubits) {
  const std::size_t dimension = linalg::basis_size(num_qubits);
  const types::real_n amplitude =
      1.0 / std::sqrt(static_cast<types::real_n>(dimension));

  types::ket state(static_cast<Eigen::Index>(dimension));

  // |->^n gives each basis state the sign (-1)^(number of one bits).
  // This matches Quira's little-endian basis indexing because only the bit
  // pattern parity matters, not the printed bit-string order.
  for (std::size_t basis = 0; basis < dimension; ++basis) {
    const types::real_n phase =
        has_odd_basis_parity(basis) ? -types::real_n{1} : types::real_n{1};

    state[static_cast<Eigen::Index>(basis)] = types::cplx_n{phase * amplitude, 0.0};
  }

  return state;
}

types::ket x_basis_state(bool value) {
  types::ket state(2);

  state[0] = types::cplx_n{constants::INV_SQRT_2, 0.0};
  state[1] = value ? types::cplx_n{-constants::INV_SQRT_2, 0.0}
                   : types::cplx_n{constants::INV_SQRT_2, 0.0};

  return state;
}

types::ket y_basis_state(bool value) {
  types::ket state(2);

  state[0] = types::cplx_n{constants::INV_SQRT_2, 0.0};
  state[1] = value ? -constants::IM * constants::INV_SQRT_2
                   : constants::IM * constants::INV_SQRT_2;

  return state;
}

types::ket bell_state(BellState state) {
  types::ket result = types::ket::Zero(4);

  switch (state) {
    case BellState::PhiPlus:
      result[0] = types::cplx_n{constants::INV_SQRT_2, 0.0};
      result[3] = types::cplx_n{constants::INV_SQRT_2, 0.0};
      break;

    case BellState::PhiMinus:
      result[0] = types::cplx_n{constants::INV_SQRT_2, 0.0};
      result[3] = types::cplx_n{-constants::INV_SQRT_2, 0.0};
      break;

    case BellState::PsiPlus:
      result[1] = types::cplx_n{constants::INV_SQRT_2, 0.0};
      result[2] = types::cplx_n{constants::INV_SQRT_2, 0.0};
      break;

    case BellState::PsiMinus:
      result[1] = types::cplx_n{constants::INV_SQRT_2, 0.0};
      result[2] = types::cplx_n{-constants::INV_SQRT_2, 0.0};
      break;
  }

  return result;
}

types::ket ghz_state(std::size_t num_qubits) {
  if (num_qubits < 2) {
    throw exception::InvalidArgument("states::ghz_state()",
                                     "GHZ state requires at least two qubits");
  }

  const std::size_t dimension = linalg::basis_size(num_qubits);

  types::ket state = types::ket::Zero(static_cast<Eigen::Index>(dimension));

  state[0] = types::cplx_n{constants::INV_SQRT_2, 0.0};
  state[static_cast<Eigen::Index>(dimension - 1)] =
      types::cplx_n{constants::INV_SQRT_2, 0.0};

  return state;
}

types::ket w_state(std::size_t num_qubits) {
  validate_non_zero_qubit(num_qubits, "w_state");

  const std::size_t dimension = linalg::basis_size(num_qubits);
  const types::real_n amplitude =
      1.0 / std::sqrt(static_cast<types::real_n>(num_qubits));

  types::ket state = types::ket::Zero(static_cast<Eigen::Index>(dimension));

  // In little-endian basis indexing, qubit q being one corresponds to the basis
  // index with bit q set, i.e. 1 << q.
  for (std::size_t i = 0; i < num_qubits; ++i) {
    const std::size_t basis = std::size_t{1} << i;
    state[static_cast<Eigen::Index>(basis)] = types::cplx_n{amplitude, 0.0};
  }

  return state;
}

types::c_mat zero_projector(std::size_t num_qubits) {
  return projector_from_state(zero_state(num_qubits));
}

types::c_mat one_projector(std::size_t num_qubits) {
  return projector_from_state(one_state(num_qubits));
}

types::c_mat basis_projector(std::size_t num_qubits, std::size_t basis_index) {
  return projector_from_state(basis_state(num_qubits, basis_index));
}

types::c_mat plus_projector(std::size_t num_qubits) {
  return projector_from_state(plus_state(num_qubits));
}

types::c_mat minus_projector(std::size_t num_qubits) {
  return projector_from_state(minus_state(num_qubits));
}

types::c_mat x_basis_projector(bool value) {
  return projector_from_state(x_basis_state(value));
}

types::c_mat y_basis_projector(bool value) {
  return projector_from_state(y_basis_state(value));
}

types::c_mat bell_projector(BellState state) {
  return projector_from_state(bell_state(state));
}

types::c_mat ghz_projector(std::size_t num_qubits) {
  return projector_from_state(ghz_state(num_qubits));
}

types::c_mat w_projector(std::size_t num_qubits) {
  return projector_from_state(w_state(num_qubits));
}

}  // namespace quira::states
