// Copyright 2026 Quira contributors
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

#include "quira/output.hpp"

#include <stdexcept>

namespace quira {

SimulationOutput::SimulationOutput(std::vector<ClassicalRegister> shots,
                                   Index shot_count)
    : shot_registers_(std::move(shots)), shots_(shot_count) {
}

Index SimulationOutput::shots() const noexcept {
  return shots_;
}

const std::vector<ClassicalRegister>& SimulationOutput::shot_registers()
    const noexcept {
  return shot_registers_;
}

const ClassicalRegister& SimulationOutput::shot(Index index) const {
  if (index >= shot_registers_.size()) {
    throw std::out_of_range("RunResult shot index is out of range");
  }
  return shot_registers_[index];
}

std::unordered_map<std::string, Index> SimulationOutput::get_counts() const {
  std::unordered_map<std::string, Index> counts;

  for (const ClassicalRegister& reg : shot_registers_) {
    const std::string key = reg.bit_string();
    if (!key.empty()) {
      ++counts[key];
    }
  }

  return counts;
}

std::unordered_map<std::string, double> SimulationOutput::get_probabilities() const {
  std::unordered_map<std::string, double> probabilities;
  const auto counts = get_counts();

  if (shots_ == 0) {
    return probabilities;
  }

  for (const auto& [outcome, count] : counts) {
    probabilities[outcome] = static_cast<double>(count) / static_cast<double>(shots_);
  }

  return probabilities;
}

std::optional<std::pair<std::string, Index>> SimulationOutput::get_most_frequent()
    const {
  const auto counts = get_counts();

  std::optional<std::pair<std::string, Index>> best;
  for (const auto& entry : counts) {
    if (!best || entry.second > best->second) {
      best = entry;
    }
  }

  return best;
}

}  // namespace quira
