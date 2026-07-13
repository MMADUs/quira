# Contributing to Quira

Quira is a C++ quantum computing library. The current implementation focuses on
a state-vector simulator and a small public API that can grow toward additional
backends later.

## Build

Use CMake from the project root:

```bash
./scripts/build.sh clean
```

For normal rebuilds:

```bash 
./scripts/build.sh
```

Run tests with:

```bash
ctest --test-dir build --output-on-failure
```

## Code Style

- Use C++20.
- Put public headers under include/quira/.
- Put implementation files under src/.
- Put library code in the quira namespace.
- Use PascalCase for types.
- Use snake_case for functions and variables.
- Use UPPER_CASE for constants.
- Prefer value types where ownership is simple.
- Use std::unique_ptr for exclusive polymorphic ownership.
- Avoid raw owning pointers.
- Avoid macros unless they are clearly necessary.
- Avoid hidden global mutable state.

## Formatting

Format C++ files with clang-format using the repository `.clang-format.`

## Static Analysis

Use `clang-tidy` as guidance. Some warnings may be adjusted while the API is
still stabilizing.

## Comments

Public headers should use Doxygen-style comments for API contracts.
Implementation files should use normal comments for non-obvious logic, especially
quantum math, indexing conventions, validation, and performance-sensitive code.

## Quantum Conventions

Quira currently follows a Qiskit-like little-endian convention for state-vector
basis indexing:

```txt 
qubit 0 -> least significant bit
targets[0] -> local bit 0
```
Multi-qubit gate matrices must match this convention.

## Exceptions

A Quira exception hierarchy exists, but wiring may happen gradually while the
core API is still stabilizing. New validation should prefer clear, intentional
error messages.

## Pull Request Checklist

Before submitting a change:
- Build succeeds.
- Tests pass.
- Public headers are documented.
- New quantum indexing behavior is explained in comments or tests.
- No generated build artifacts are committed.
