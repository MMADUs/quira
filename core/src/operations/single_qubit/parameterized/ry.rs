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
use crate::{Complex, GateType, Matrix, Qubit};

#[derive(Debug, Clone)]
/// Represents a rotation around the Y-axis of the Bloch sphere.
///
/// This gate rotates the qubit state around the Y-axis by an angle theta.
///
/// The matrix form is:
///
/// RY(θ) = [ [ cos(θ/2), -sin(θ/2) ],
///           [ sin(θ/2), cos(θ/2)  ] ]
///
/// This gate is equivalent to e^(-i*θ*Y/2) where Y is the Pauli-Y matrix.
pub struct RotateY {
    target: Qubit,
    theta: f64,
}

impl RotateY {
    pub fn new(target: Qubit, theta: f64) -> Self {
        Self { target, theta }
    }
}

impl QuantumGate for RotateY {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(-1.0 * s, 0.0)],
            [Complex::new(s, 0.0), Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        format!("RY(target={}, theta={:.4})", self.target, self.theta)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::RotateY(Self::new(self.target, self.theta)))
    }
}

impl SingleQubitGate for RotateY {
    fn alpha_re(&self) -> f64 {
        (self.theta / 2.0).cos()
    }

    fn alpha_im(&self) -> f64 {
        0.0
    }

    fn beta_re(&self) -> f64 {
        (self.theta / 2.0).sin()
    }

    fn beta_im(&self) -> f64 {
        0.0
    }

    fn global_phase(&self) -> f64 {
        0.0
    }
}
