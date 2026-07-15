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

#pragma once

/**
 * @file quira/quira.hpp
 * @brief Main public header for the Quira C++ API.
 *
 * Include this header when using the full Quira public API. Projects that care
 * about compile times may include narrower headers directly.
 */

#include "quira/backend/statevec.hpp"
#include "quira/circuit.hpp"
#include "quira/constants.hpp"
#include "quira/exception.hpp"
#include "quira/gate.hpp"
#include "quira/gates/param.hpp"
#include "quira/gates/single.hpp"
#include "quira/gates/three.hpp"
#include "quira/gates/two.hpp"
#include "quira/output.hpp"
#include "quira/register.hpp"
#include "quira/types.hpp"
