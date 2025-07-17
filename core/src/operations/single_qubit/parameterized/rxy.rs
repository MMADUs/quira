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
use crate::operations::GateType;
use crate::operations::single_qubit::SingleQubitType;
use crate::ops::QuantumGate;
use crate::types::{Complex, Matrix};

#[derive(Debug, Clone)]
/// Represents a rotation in the XY-plane of the Bloch sphere.
///
/// This gate rotates the qubit state by angle theta in the direction specified by angle phi in the XY-plane.
///
/// The matrix form is:
///
/// RXY(θ, φ) = [ [ cos(θ/2)                   , -sin(θ/2)(sin(φ) + i*cos(φ)) ],
///               [ sin(θ/2)(sin(φ) - i*cos(φ)), cos(θ/2)                     ] ]
///
/// This gate combines aspects of both RX and RY rotations and allows for arbitrary rotations in the XY-plane.
pub struct RotateXY {
    target: QuantumBit,
    theta: f64,
    phi: f64,
}

impl RotateXY {
    pub fn new(target: &QuantumBit, theta: f64, phi: f64) -> Self {
        Self {
            target: target.clone(),
            theta,
            phi,
        }
    }
}

impl QuantumGate for RotateXY {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        let vx: f64 = (self.phi).cos();
        let vy: f64 = (self.phi).sin();
        array![
            [
                Complex::new(c, 0.0),
                Complex::new(-1.0 * s * vy, -1.0 * s * vx)
            ],
            [Complex::new(s * vy, -1.0 * s * vx), Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        format!(
            "RXY(target={}, theta={:.4}, phi={:.4})",
            self.target, self.theta, self.phi
        )
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::RotateXY(Self::new(
            &self.target,
            self.theta,
            self.phi,
        )))
    }
}
