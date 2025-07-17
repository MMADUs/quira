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
use crate::operations::two_qubit::TwoQubitType;
use crate::ops::QuantumGate;
use crate::types::{Complex, Matrix};

#[derive(Debug, Clone)]
/// Represents a controlled rotation in the XY-plane gate.
///
/// This gate performs a rotation in the XY-plane on the target qubit when
/// the control qubit is in the |1⟩ state. The rotation is parameterized by
/// angle θ and axis direction φ in the XY-plane.
///
/// The matrix form involves rotations around the axis (cos(φ), sin(φ), 0):
/// CRXY = [ [ 1, 0,               0,                 0 ],
///          [ 0, 1,               0,                 0 ],
///          [ 0, 0,        cos(θ/2), -sin(θ/2)·e^(-iφ) ],
///          [ 0, 0, sin(θ/2)·e^(iφ),          cos(θ/2) ] ]
///
/// # Parameters
/// - `control`: The control qubit index
/// - `target`: The target qubit index
/// - `theta`: The rotation angle θ in radians
/// - `phi`: The axis angle φ in the XY-plane
pub struct ControlledRotateXY {
    control: QuantumBit,
    target: QuantumBit,
    theta: f64,
    phi: f64,
}

impl ControlledRotateXY {
    pub fn new(control: &QuantumBit, target: &QuantumBit, theta: f64, phi: f64) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
            theta,
            phi,
        }
    }
}

impl QuantumGate for ControlledRotateXY {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        let vx: f64 = (self.phi).cos();
        let vy: f64 = (self.phi).sin();
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(c, 0.0),
                Complex::new(-1.0 * s * vy, -1.0 * s * vx),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(s * vy, -1.0 * s * vx),
                Complex::new(c, 0.0),
            ]
        ]
    }

    fn name(&self) -> String {
        format!(
            "CRXY(control={}, target={}, theta={}, phi={})",
            self.control, self.target, self.theta, self.phi
        )
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::ControlledRotateXY(Self::new(
            &self.control,
            &self.target,
            self.theta,
            self.phi,
        )))
    }
}
