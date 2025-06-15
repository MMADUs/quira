//! Copyright (c) 2024-2025 Quira, Inc.
//!
//! This file is part of Quira
//!
//! This program is free software: you can redistribute it and/or modify
//! it under the terms of the GNU Affero General Public License as published by
//! the Free Software Foundation, either version 3 of the License, or
//! (at your option) any later version.
//!
//! This program is distributed in the hope that it will be useful
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//! GNU Affero General Public License for more details.
//!
//! You should have received a copy of the GNU Affero General Public License
//! along with this program.  If not, see <http://www.gnu.org/licenses/>.

use core::f64;

use crate::types::Complex;

/// Boltzmann constant ('k') = 1.380649 × 10^-23 J/K  
/// Relates the average kinetic energy of particles in a gas with temperature  
pub const BOLTZMANN: f64 = 1.380649e-23; // J/K

/// Planck constant ('h') = 6.62607015 × 10^-34 J⋅s  
/// Fundamental constant relating energy and frequency: E = h·f  
pub const PLANCK: f64 = 6.62607015e-34; // J⋅s

/// mathematical constant of pi ('π') = 3.14159265358979323846264338327950288_f64
pub const PI: f64 = std::f64::consts::PI;

/// Mathematical constant for Euler's number ('e') = 2.71828182845904523536028747135266250_f64
pub const E: f64 = std::f64::consts::E;

/// Mathematical constant for regular Epsilon ('ε') = 1e-10 = 1 * 10^-10 = 0.0000000001
pub const EPSILON: f64 = 1e-10;

/// Mathematical constant for machine Epsilon ('ε') = 2.220446049250313e-16
pub const M_EPSILON: f64 = f64::EPSILON;

/// The largest finite representable number (before it overflows to infinity) = 1.7976931348623157E+308f64
pub const MAX_FLOAT: f64 = f64::MAX;

/// Mathematical constant for positive infinity (`+∞`), a value exceeds the largest representable finite number.
pub const INF: f64 = f64::INFINITY;

/// Mathematical constant for negative infinity (`-∞`), a value is below the smallest representable finite number.
pub const NEG_INF: f64 = f64::NEG_INFINITY;

/// The imaginary unit (i), where i² = -1
pub const IM: Complex = Complex::new(0.0, 1.0);

/// The complex number 1 + 0i (the real number 1 as a complex number)
pub const C_ONE: Complex = Complex::new(1.0, 0.0);

/// The complex number 0 + 0i (the real number 0 as a complex number)
pub const C_ZERO: Complex = Complex::new(0.0, 0.0);
