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
/// Represents a rotation around the X-axis of the Bloch sphere.
///
/// This gate rotates the qubit state around the X-axis by an angle theta.
///
/// The matrix form is:
///
/// RX(θ) = [ [ cos(θ/2)   , -i*sin(θ/2) ],
///           [ -i*sin(θ/2), cos(θ/2)    ] ]
///
/// This gate is equivalent to e^(-i*θ*X/2) where X is the Pauli-X matrix.
pub struct RotateX {
    target: QuantumBit,
    theta: f64,
}

impl RotateX {
    pub fn new(target: &QuantumBit, theta: f64) -> Self {
        Self {
            target: target.clone(),
            theta,
        }
    }
}

impl QuantumGate for RotateX {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(0.0, -1.0 * s)],
            [Complex::new(0.0, -1.0 * s), Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        format!("RX(target={}, theta={:.4})", self.target, self.theta)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::RotateX(Self::new(
            &self.target,
            self.theta,
        )))
    }
}
