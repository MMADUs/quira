/*
Copyright (c) 2024-2025 Quira, Inc.

This file is part of Quira

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

use ndarray::array;

use crate::operations::QuantumGate;
use crate::operations::single_qubit::{SingleQubitGate, SingleQubitType};
use crate::{Complex, GateType, Matrix, Qubit, constant::PI};

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
///
///   = [ [ 1, 0                     ],
///       [ 0, cos(π/4) + i*sin(π/4) ] ]
///
///   = [ [ 1, 0           ],
///       [ 0, 1/√2 + i/√2 ] ]
///
/// The T gate is not a Clifford gate but is needed along with Clifford gates
/// to achieve universal quantum computation. It's also known as the π/8 gate
/// (because its effect is a rotation of π/4 which is π/8 radians).
pub struct TGate {
    target: Qubit,
}

impl TGate {
    pub fn new(target: Qubit) -> Self {
        Self { target }
    }
}

impl QuantumGate for TGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [
                Complex::new(0.0, 0.0),
                Complex::new((PI / 4.0).cos(), (PI / 4.0).sin())
            ]
        ]
    }

    fn name(&self) -> String {
        format!("T(target={})", self.target)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::TGate(Self::new(self.target)))
    }
}

impl SingleQubitGate for TGate {
    fn alpha_re(&self) -> f64 {
        (PI / 8.0).cos()
    }

    fn alpha_im(&self) -> f64 {
        (-1.0) * (PI / 8.0).sin()
    }

    fn beta_re(&self) -> f64 {
        0.0
    }

    fn beta_im(&self) -> f64 {
        0.0
    }

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
///
///    = [ [ 1, 0                     ],
///        [ 0, cos(π/4) - i*sin(π/4) ] ]
///
///    = [ [ 1, 0           ],
///        [ 0, 1/√2 - i/√2 ] ]
///
/// The T† gate is used to undo the effect of the T gate and is also an important
/// component in many quantum algorithms and error correction protocols.
pub struct InvTGate {
    target: Qubit,
}

impl InvTGate {
    pub fn new(target: Qubit) -> Self {
        Self { target }
    }
}

impl QuantumGate for InvTGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [
                Complex::new(0.0, 0.0),
                Complex::new((PI / 4.0).cos(), -1.0 * (PI / 4.0).sin())
            ]
        ]
    }

    fn name(&self) -> String {
        format!("Inv-T(target={})", self.target)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::InvTGate(Self::new(self.target)))
    }
}

impl SingleQubitGate for InvTGate {
    fn alpha_re(&self) -> f64 {
        (PI / 8.0).cos()
    }

    fn alpha_im(&self) -> f64 {
        (PI / 8.0).sin()
    }

    fn beta_re(&self) -> f64 {
        0.0
    }

    fn beta_im(&self) -> f64 {
        0.0
    }

    fn global_phase(&self) -> f64 {
        (-1.0 * PI) / 8.0
    }
}
