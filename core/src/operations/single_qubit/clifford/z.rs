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

use crate::bit::QuantumBit;
use crate::constant::PI;
use crate::operations::GateType;
use crate::operations::single_qubit::SingleQubitType;
use crate::ops::QuantumGate;
use crate::types::{Complex, Matrix};

#[derive(Debug, Clone)]
/// Represents the Pauli-Z gate (also known as the phase-flip gate).
///
/// This gate leaves the |0⟩ state unchanged and flips the sign of the |1⟩ state.
/// It represents a rotation around the Z-axis by π radians.
///
/// The matrix form is:
///
/// Z = [ [ 1,  0 ],
///       [ 0, -1 ] ]
pub struct PauliZ {
    target: QuantumBit,
}

impl PauliZ {
    pub fn new(target: &QuantumBit) -> Self {
        Self {
            target: target.clone(),
        }
    }

    pub fn power(&self, exponent: f64) -> ZPowGate {
        ZPowGate::new(&self.target, exponent)
    }
}

impl QuantumGate for PauliZ {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), Complex::new(-1.0, 0.0)]
        ]
    }

    fn name(&self) -> String {
        format!("Z(target={})", self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::PauliZ(Self::new(&self.target)))
    }
}

#[derive(Debug, Clone)]
/// Represents the Z^t gate (Z to the power of t).
///
/// Z^t = [ [ e^(-iπt/2), 0         ],
///         [ 0,          e^(iπt/2) ] ]
///
/// When t=1, this reduces to the standard Pauli-Z gate.
pub struct ZPowGate {
    target: QuantumBit,
    exponent: f64, // t parameter
}

impl ZPowGate {
    pub fn new(target: &QuantumBit, exponent: f64) -> Self {
        Self {
            target: target.clone(),
            exponent,
        }
    }
}

impl QuantumGate for ZPowGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let angle = PI * self.exponent / 2.0;
        let exp_neg_i_angle = Complex::new(0.0, -angle).exp();
        let exp_pos_i_angle = Complex::new(0.0, angle).exp();

        array![
            [exp_neg_i_angle, Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), exp_pos_i_angle]
        ]
    }

    fn name(&self) -> String {
        format!("Z^{:.3}(target={})", self.exponent, self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::ZPowGate(Self::new(
            &self.target,
            self.exponent,
        )))
    }
}
