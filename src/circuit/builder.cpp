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

#include "quira/builder.hpp"

#include "quira/constants.hpp"
#include "quira/exception.hpp"
#include "quira/gates/single.hpp"
#include "quira/gates/two.hpp"
#include "quira/unitary.hpp"

#include <cmath>
#include <sstream>
#include <unordered_set>

namespace quira::algorithms {
namespace {

void validate_qubit_list(const char* caller, const std::vector<types::qubit>& qubits) {
  if (qubits.empty()) {
    throw exception::InvalidArgument(caller, "requires at least one qubit");
  }

  std::unordered_set<types::qubit> seen;
  for (types::qubit qubit : qubits) {
    if (!seen.insert(qubit).second) {
      std::ostringstream oss;
      oss << "received duplicate qubit index " << static_cast<std::size_t>(qubit);
      throw exception::InvalidArgument(caller, oss.str());
    }
  }
}

void validate_at_least_two(const char* caller,
                           const std::vector<types::qubit>& qubits) {
  validate_qubit_list(caller, qubits);

  if (qubits.size() < 2) {
    throw exception::InvalidArgument(caller, "requires at least two qubits");
  }
}

types::c_mat controlled_phase_matrix(types::real_n theta) {
  types::c_mat matrix = types::c_mat::Identity(4, 4);
  matrix(3, 3) = std::exp(constants::IM * theta);
  return matrix;
}

std::string controlled_phase_name(types::qubit control, types::qubit target,
                                  types::real_n theta) {
  std::ostringstream oss;
  oss << "CP(control=" << static_cast<std::size_t>(control)
      << ", target=" << static_cast<std::size_t>(target) << ", theta=" << theta << ")";
  return oss.str();
}

Unitary controlled_phase(types::qubit control, types::qubit target,
                         types::real_n theta) {
  return Unitary{
      controlled_phase_name(control, target, theta),
      controlled_phase_matrix(theta),
      {control, target},
  };
}

types::real_n qft_phase(std::size_t distance, bool inverse) {
  const types::real_n denominator =
      std::pow(types::real_n{2.0}, static_cast<types::real_n>(distance));
  const types::real_n sign = inverse ? types::real_n{-1.0} : types::real_n{1.0};

  return sign * types::real_n{2.0} * constants::PI / denominator;
}

}  // namespace

QuantumCircuit& qft(QuantumCircuit& circuit, const std::vector<types::qubit>& qubits,
                    bool swap_output) {
  validate_qubit_list("algorithms::qft()", qubits);

  const std::size_t count = qubits.size();

  for (std::size_t index = 0; index < count; ++index) {
    circuit.add(H{qubits[index]});

    for (std::size_t distance = 2; distance <= count - index; ++distance) {
      const types::qubit control = qubits[index + distance - 1];
      const types::qubit target = qubits[index];
      circuit.add(controlled_phase(control, target, qft_phase(distance, false)));
    }
  }

  if (swap_output) {
    for (std::size_t index = 0; index < count / 2; ++index) {
      circuit.add(Swap{qubits[index], qubits[count - index - 1]});
    }
  }

  return circuit;
}

QuantumCircuit& inverse_qft(QuantumCircuit& circuit,
                            const std::vector<types::qubit>& qubits, bool swap_output) {
  validate_qubit_list("algorithms::inverse_qft()", qubits);

  const std::size_t count = qubits.size();

  if (swap_output) {
    for (std::size_t index = 0; index < count / 2; ++index) {
      circuit.add(Swap{qubits[index], qubits[count - index - 1]});
    }
  }

  for (std::size_t reverse = count; reverse > 0; --reverse) {
    const std::size_t index = reverse - 1;

    for (std::size_t distance = count - index; distance > 1; --distance) {
      const types::qubit control = qubits[index + distance - 1];
      const types::qubit target = qubits[index];
      circuit.add(controlled_phase(control, target, qft_phase(distance, true)));
    }

    circuit.add(H{qubits[index]});
  }

  return circuit;
}

QuantumCircuit& ghz(QuantumCircuit& circuit, const std::vector<types::qubit>& qubits) {
  validate_at_least_two("algorithms::ghz()", qubits);

  circuit.add(H{qubits.front()});

  for (std::size_t index = 1; index < qubits.size(); ++index) {
    circuit.add(CNOT{qubits.front(), qubits[index]});
  }

  return circuit;
}

QuantumCircuit& bell(QuantumCircuit& circuit, types::qubit first, types::qubit second) {
  if (first == second) {
    throw exception::InvalidArgument("algorithms::bell()", "requires distinct qubits");
  }

  circuit.add(H{first});
  circuit.add(CNOT{first, second});

  return circuit;
}

}  // namespace quira::algorithms
