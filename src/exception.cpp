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

#include "quira/exception.hpp"

#include <utility>

namespace quira {

const char* to_string(ErrorCode code) noexcept {
  switch (code) {
    case ErrorCode::Unknown:
      return "Unknown";
    case ErrorCode::InvalidArgument:
      return "InvalidArgument";
    case ErrorCode::OutOfRange:
      return "OutOfRange";
    case ErrorCode::InvalidQubit:
      return "InvalidQubit";
    case ErrorCode::InvalidClbit:
      return "InvalidClbit";
    case ErrorCode::InvalidGate:
      return "InvalidGate";
    case ErrorCode::InvalidCircuit:
      return "InvalidCircuit";
    case ErrorCode::InvalidState:
      return "InvalidState";
    case ErrorCode::DimensionMismatch:
      return "DimensionMismatch";
    case ErrorCode::DuplicateQubit:
      return "DuplicateQubit";
    case ErrorCode::NullPointer:
      return "NullPointer";
    case ErrorCode::SimulationError:
      return "SimulationError";
    case ErrorCode::NotImplemented:
      return "NotImplemented";
  }

  return "Unknown";
}

QuiraException::QuiraException(ErrorCode code, std::string message)
    : QuiraException(code, "", std::move(message)) {
}

QuiraException::QuiraException(ErrorCode code, std::string where, std::string message)
    : std::runtime_error(format_message(code, where, message)),
      code_(code),
      where_(std::move(where)),
      message_(std::move(message)) {
}

ErrorCode QuiraException::code() const noexcept {
  return code_;
}

const std::string& QuiraException::where() const noexcept {
  return where_;
}

const std::string& QuiraException::message() const noexcept {
  return message_;
}

std::string QuiraException::format_message(ErrorCode code, const std::string& where,
                                           const std::string& message) {
  std::string formatted = "Quira ";
  formatted += to_string(code);

  if (!where.empty()) {
    formatted += " in ";
    formatted += where;
  }

  formatted += ": ";
  formatted += message;
  return formatted;
}

InvalidArgument::InvalidArgument(std::string message)
    : QuiraException(ErrorCode::InvalidArgument, std::move(message)) {
}

InvalidArgument::InvalidArgument(std::string where, std::string message)
    : QuiraException(ErrorCode::InvalidArgument, std::move(where), std::move(message)) {
}

OutOfRange::OutOfRange(std::string message)
    : QuiraException(ErrorCode::OutOfRange, std::move(message)) {
}

OutOfRange::OutOfRange(std::string where, std::string message)
    : QuiraException(ErrorCode::OutOfRange, std::move(where), std::move(message)) {
}

InvalidQubit::InvalidQubit(std::string message)
    : QuiraException(ErrorCode::InvalidQubit, std::move(message)) {
}

InvalidQubit::InvalidQubit(std::string where, std::string message)
    : QuiraException(ErrorCode::InvalidQubit, std::move(where), std::move(message)) {
}

InvalidClbit::InvalidClbit(std::string message)
    : QuiraException(ErrorCode::InvalidClbit, std::move(message)) {
}

InvalidClbit::InvalidClbit(std::string where, std::string message)
    : QuiraException(ErrorCode::InvalidClbit, std::move(where), std::move(message)) {
}

InvalidGate::InvalidGate(std::string message)
    : QuiraException(ErrorCode::InvalidGate, std::move(message)) {
}

InvalidGate::InvalidGate(std::string where, std::string message)
    : QuiraException(ErrorCode::InvalidGate, std::move(where), std::move(message)) {
}

InvalidCircuit::InvalidCircuit(std::string message)
    : QuiraException(ErrorCode::InvalidCircuit, std::move(message)) {
}

InvalidCircuit::InvalidCircuit(std::string where, std::string message)
    : QuiraException(ErrorCode::InvalidCircuit, std::move(where), std::move(message)) {
}

InvalidState::InvalidState(std::string message)
    : QuiraException(ErrorCode::InvalidState, std::move(message)) {
}

InvalidState::InvalidState(std::string where, std::string message)
    : QuiraException(ErrorCode::InvalidState, std::move(where), std::move(message)) {
}

DimensionMismatch::DimensionMismatch(std::string message)
    : QuiraException(ErrorCode::DimensionMismatch, std::move(message)) {
}

DimensionMismatch::DimensionMismatch(std::string where, std::string message)
    : QuiraException(ErrorCode::DimensionMismatch, std::move(where),
                     std::move(message)) {
}

DuplicateQubit::DuplicateQubit(std::string message)
    : QuiraException(ErrorCode::DuplicateQubit, std::move(message)) {
}

DuplicateQubit::DuplicateQubit(std::string where, std::string message)
    : QuiraException(ErrorCode::DuplicateQubit, std::move(where), std::move(message)) {
}

NullPointer::NullPointer(std::string message)
    : QuiraException(ErrorCode::NullPointer, std::move(message)) {
}

NullPointer::NullPointer(std::string where, std::string message)
    : QuiraException(ErrorCode::NullPointer, std::move(where), std::move(message)) {
}

SimulationError::SimulationError(std::string message)
    : QuiraException(ErrorCode::SimulationError, std::move(message)) {
}

SimulationError::SimulationError(std::string where, std::string message)
    : QuiraException(ErrorCode::SimulationError, std::move(where), std::move(message)) {
}

NotImplemented::NotImplemented(std::string message)
    : QuiraException(ErrorCode::NotImplemented, std::move(message)) {
}

NotImplemented::NotImplemented(std::string where, std::string message)
    : QuiraException(ErrorCode::NotImplemented, std::move(where), std::move(message)) {
}

}  // namespace quira
