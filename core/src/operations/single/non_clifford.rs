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

use super::SingleQubit;
use crate::constant::PI;
use crate::operations::QuantumGate;
use crate::types::{Complex, Matrix, Qubit};
use ndarray::array;

#[derive(Debug, Clone)]
/// Represents the T gate, which is a pi/4 phase rotation gate.
///
/// The T gate applies a π/4 phase shift to the |1⟩ state while leaving the |0⟩ state unchanged.
/// It is a special case of the phase shift gate and is commonly used in quantum circuits.
///
/// The matrix form is:
///
/// T = [ [ 1, 0        ],
///       [ 0, e^(iπ/4) ] ]
///   = [ [ 1, 0                     ],
///       [ 0, cos(π/4) + i*sin(π/4) ] ]
///   = [ [ 1, 0           ],
///       [ 0, 1/√2 + i/√2 ] ]
///
/// The T gate is not a Clifford gate but is needed along with Clifford gates
/// to achieve universal quantum computation. It's also known as the π/8 gate
/// (because its effect is a rotation of π/4 which is π/8 radians).
pub struct TGate {
    qubit: Qubit,
}

impl TGate {
    pub fn new(qubit: Qubit) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for TGate {
    /// construct the 2x2 unitary matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [
                Complex::new(0.0, 0.0),
                Complex::new((PI / 4.0).cos(), (PI / 4.0).sin())
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        String::from("T")
    }

    /// returns the index of the qubit this gate operates on.
    fn target_qubit(&self) -> Qubit {
        self.qubit
    }
}

impl SingleQubit for TGate {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        (PI / 8.0).cos()
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        (-1.0) * (PI / 8.0).sin()
    }

    /// returns the real part of beta
    fn beta_re(&self) -> f64 {
        0.0
    }

    /// returns the imaginary part of beta
    fn beta_im(&self) -> f64 {
        0.0
    }

    /// returns the global phase
    fn global_phase(&self) -> f64 {
        PI / 8.0
    }
}

#[derive(Debug, Clone)]
/// Represents the inverse T gate (T†), which is a -pi/4 phase rotation gate.
///
/// The inverse T gate applies a -π/4 phase shift to the |1⟩ state while leaving
/// the |0⟩ state unchanged. It is the Hermitian conjugate (inverse) of the T gate.
///
/// The matrix form is:
///
/// T† = [ [ 1, 0         ],
///        [ 0, e^(-iπ/4) ] ]
///    = [ [ 1, 0                     ],
///        [ 0, cos(π/4) - i*sin(π/4) ] ]
///    = [ [ 1, 0           ],
///        [ 0, 1/√2 - i/√2 ] ]
///
/// The T† gate is used to undo the effect of the T gate and is also an important
/// component in many quantum algorithms and error correction protocols.
pub struct InvTGate {
    qubit: Qubit,
}

impl InvTGate {
    pub fn new(qubit: Qubit) -> Self {
        Self { qubit }
    }
}

impl QuantumGate for InvTGate {
    /// construct the 2x2 unitary matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [
                Complex::new(0.0, 0.0),
                Complex::new((PI / 4.0).cos(), -1.0 * (PI / 4.0).sin())
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        String::from("Inv-T")
    }

    /// returns the index of the qubit this gate operates on.
    fn target_qubit(&self) -> Qubit {
        self.qubit
    }
}

impl SingleQubit for InvTGate {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        (PI / 8.0).cos()
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        (PI / 8.0).sin()
    }

    /// returns the real part of beta
    fn beta_re(&self) -> f64 {
        0.0
    }

    /// returns the imaginary part of beta
    fn beta_im(&self) -> f64 {
        0.0
    }

    /// returns the global phase
    fn global_phase(&self) -> f64 {
        (-1.0) * PI / 8.0
    }
}
