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

#pragma once

#include <stdexcept>
#include <string>

namespace quira {

/**
 * @brief Stable error categories used by Quira exceptions.
 *
 * ErrorCode gives users a lightweight way to branch on failure categories
 * without parsing exception messages.
 */
enum class ErrorCode {
  Unknown,
  InvalidArgument,
  OutOfRange,
  InvalidQubit,
  InvalidClbit,
  InvalidGate,
  InvalidCircuit,
  InvalidState,
  DimensionMismatch,
  DuplicateQubit,
  NullPointer,
  SimulationError,
  NotImplemented
};

/**
 * @brief Returns a stable string name for an error code.
 *
 * @param code Error category.
 * @return String representation of the error category.
 */
const char* to_string(ErrorCode code) noexcept;

/**
 * @brief Base exception type for Quira library errors.
 *
 * All Quira-specific exceptions derive from this class, so users can catch
 * `QuiraException` to handle library errors without depending on every
 * concrete exception type.
 */
class QuiraException : public std::runtime_error {
public:
  /**
   * @brief Creates a Quira exception without a source-location label.
   *
   * @param code Stable error category.
   * @param message Human-readable error message.
   */
  QuiraException(ErrorCode code, std::string message);

  /**
   * @brief Creates a Quira exception with a source-location label.
   *
   * @param code Stable error category.
   * @param where Name of the API or subsystem where the error occurred.
   * @param message Human-readable error message.
   */
  QuiraException(ErrorCode code, std::string where, std::string message);

  /**
   * @brief Returns the stable error category.
   *
   * @return Error code associated with the exception.
   */
  ErrorCode code() const noexcept;

  /**
   * @brief Returns the source-location label, if one was provided.
   *
   * @return API or subsystem name associated with the failure.
   */
  const std::string& where() const noexcept;

  /**
   * @brief Returns the original human-readable message.
   *
   * @return Message provided when constructing the exception.
   */
  const std::string& message() const noexcept;

private:
  ErrorCode code_{ErrorCode::Unknown};
  std::string where_;
  std::string message_;

  static std::string format_message(ErrorCode code, const std::string& where,
                                    const std::string& message);
};

/**
 * @brief Generic invalid argument error.
 */
class InvalidArgument : public QuiraException {
public:
  explicit InvalidArgument(std::string message);
  InvalidArgument(std::string where, std::string message);
};

/**
 * @brief Generic out-of-range error.
 */
class OutOfRange : public QuiraException {
public:
  explicit OutOfRange(std::string message);
  OutOfRange(std::string where, std::string message);
};

/**
 * @brief Error for invalid qubit indices or qubit-related inputs.
 */
class InvalidQubit : public QuiraException {
public:
  explicit InvalidQubit(std::string message);
  InvalidQubit(std::string where, std::string message);
};

/**
 * @brief Error for invalid classical-bit indices or classical inputs.
 */
class InvalidClbit : public QuiraException {
public:
  explicit InvalidClbit(std::string message);
  InvalidClbit(std::string where, std::string message);
};

/**
 * @brief Error for invalid gate objects, gate targets, or gate matrices.
 */
class InvalidGate : public QuiraException {
public:
  explicit InvalidGate(std::string message);
  InvalidGate(std::string where, std::string message);
};

/**
 * @brief Error for invalid circuit structure or circuit operations.
 */
class InvalidCircuit : public QuiraException {
public:
  explicit InvalidCircuit(std::string message);
  InvalidCircuit(std::string where, std::string message);
};

/**
 * @brief Error for invalid quantum state vectors or backend state.
 */
class InvalidState : public QuiraException {
public:
  explicit InvalidState(std::string message);
  InvalidState(std::string where, std::string message);
};

/**
 * @brief Error for incompatible matrix, vector, register, or subsystem sizes.
 */
class DimensionMismatch : public QuiraException {
public:
  explicit DimensionMismatch(std::string message);
  DimensionMismatch(std::string where, std::string message);
};

/**
 * @brief Error for duplicate qubit indices where distinct qubits are required.
 */
class DuplicateQubit : public QuiraException {
public:
  explicit DuplicateQubit(std::string message);
  DuplicateQubit(std::string where, std::string message);
};

/**
 * @brief Error for null owning pointers passed into public APIs.
 */
class NullPointer : public QuiraException {
public:
  explicit NullPointer(std::string message);
  NullPointer(std::string where, std::string message);
};

/**
 * @brief Error raised while executing a circuit on a simulator backend.
 */
class SimulationError : public QuiraException {
public:
  explicit SimulationError(std::string message);
  SimulationError(std::string where, std::string message);
};

/**
 * @brief Error for intentionally unavailable functionality.
 */
class NotImplemented : public QuiraException {
public:
  explicit NotImplemented(std::string message);
  NotImplemented(std::string where, std::string message);
};

}  // namespace quira
