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

#include "quira/instruction.hpp"

#include "quira/exception.hpp"

namespace quira {

GateInstruction::GateInstruction(std::unique_ptr<QuantumGate> gate)
    : gate_(std::move(gate)) {
  if (!gate_) {
    throw exception::InvalidArgument("GateInstruction::GateInstruction()",
                                     "requires a non-null gate");
  }
}

GateInstruction::GateInstruction(std::unique_ptr<QuantumGate> gate,
                                 ClassicalCondition condition)
    : gate_(std::move(gate)), condition_(condition) {
  if (!gate_) {
    throw exception::InvalidArgument("GateInstruction::GateInstruction()",
                                     "requires a non-null gate");
  }
}

GateInstruction::GateInstruction(const GateInstruction& other)
    : gate_(other.gate_ ? other.gate_->clone() : nullptr),
      condition_(other.condition_) {
}

GateInstruction& GateInstruction::operator=(const GateInstruction& other) {
  if (this == &other) {
    return *this;
  }

  gate_ = other.gate_ ? other.gate_->clone() : nullptr;
  condition_ = other.condition_;

  return *this;
}

const QuantumGate& GateInstruction::gate() const noexcept {
  return *gate_;
}

const std::optional<ClassicalCondition>& GateInstruction::condition() const noexcept {
  return condition_;
}

}  // namespace quira
