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
/// Represents a generalized π (Pi) pulse quantum gate.
///
/// This gate flips the state of a qubit around an axis in the XY plane of the Bloch sphere.
/// The angle theta determines the orientation of the rotation axis in the XY plane.
///
/// The matrix form is:
///
/// GPi(θ) = [ [ 0              , cos(θ)-i*sin(θ) ],
///            [ cos(θ)+i*sin(θ), 0               ] ]
///
/// When θ=0, this becomes equivalent to the Pauli-X gate.
pub struct GPi {
    target: QuantumBit,
    theta: f64,
}

impl GPi {
    pub fn new(target: &QuantumBit, theta: f64) -> Self {
        Self {
            target: target.clone(),
            theta,
        }
    }
}

impl QuantumGate for GPi {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta).cos();
        let s: f64 = (self.theta).sin();
        array![
            [Complex::new(0.0, 0.0), Complex::new(c, -1.0 * s)],
            [Complex::new(c, s), Complex::new(0.0, 0.0)]
        ]
    }

    fn name(&self) -> String {
        format!("GPi(target={}, theta={:.4})", self.target, self.theta)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::GPi(Self::new(&self.target, self.theta)))
    }
}

#[derive(Debug, Clone)]
/// Represents a generalized π/2 (Pi/2) pulse quantum gate.
///
/// This gate performs a 90-degree rotation around an axis in the XY plane of the Bloch sphere.
/// The angle theta determines the orientation of the rotation axis in the XY plane.
///
/// The matrix form is:
///
/// GPi2(θ) = (1/√2) * [ [ 1              , -sin(θ)-i*cos(θ) ],
///                      [ sin(θ)-i*cos(θ), 1                ] ]
///
/// When θ=0, this is equivalent to a 90-degree rotation around the X-axis (√X gate without the global phase).
pub struct GPi2 {
    target: QuantumBit,
    theta: f64,
}

impl GPi2 {
    pub fn new(target: &QuantumBit, theta: f64) -> Self {
        Self {
            target: target.clone(),
            theta,
        }
    }
}

impl QuantumGate for GPi2 {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta).cos();
        let s: f64 = (self.theta).sin();
        array![
            [Complex::new(1.0, 0.0), Complex::new(-1.0 * s, -1.0 * c)],
            [Complex::new(s, -1.0 * c), Complex::new(1.0, 0.0)]
        ] / 2.0_f64.sqrt()
    }

    fn name(&self) -> String {
        format!("GPi2(target={}, theta={:.4})", self.target, self.theta)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::GPi2(Self::new(&self.target, self.theta)))
    }
}
