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
/// Represents the Pauli-Y gate.
///
/// This gate performs a rotation around the Y-axis by π radians,
/// with an additional phase factor.
///
/// The matrix form is:
///
/// Y = [ [ 0, -i ],
///       [ i,  0 ] ]
pub struct PauliY {
    target: QuantumBit,
}

impl PauliY {
    pub fn new(target: &QuantumBit) -> Self {
        Self {
            target: target.clone(),
        }
    }

    pub fn sqrt(&self) -> SqrtPauliY {
        SqrtPauliY::new(&self.target)
    }

    pub fn inv_sqrt(&self) -> InvSqrtPauliY {
        InvSqrtPauliY::new(&self.target)
    }

    pub fn power(&self, exponent: f64) -> YPowGate {
        YPowGate::new(&self.target, exponent)
    }
}

impl QuantumGate for PauliY {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(0.0, 0.0), Complex::new(0.0, -1.0)],
            [Complex::new(0.0, 1.0), Complex::new(0.0, 0.0)]
        ]
    }

    fn name(&self) -> String {
        format!("Y(target={})", self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::PauliY(Self::new(&self.target)))
    }
}

#[derive(Debug, Clone)]
/// Represents the square root of the Pauli-Y gate.
///
/// This gate performs half the action of the Pauli-Y gate, equivalent to
/// a 90-degree rotation around the Y-axis.
///
/// Applying this gate twice is equivalent to applying the Y gate once.
pub struct SqrtPauliY {
    target: QuantumBit,
}

impl SqrtPauliY {
    pub fn new(target: &QuantumBit) -> Self {
        Self {
            target: target.clone(),
        }
    }
}

impl QuantumGate for SqrtPauliY {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = PI / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(-1.0 * s, 0.0)],
            [Complex::new(s, 0.0), Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        format!("Sqrt-Y(target={})", self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::SqrtPauliY(Self::new(&self.target)))
    }
}

#[derive(Debug, Clone)]
/// Represents the inverse of the square root of the Pauli-Y gate.
///
/// This gate is the adjoint (inverse) of the √Y gate, equivalent to
/// a -90-degree rotation around the Y-axis.
///
/// The matrix form is the complex conjugate transpose of the √Y matrix.
///
/// Applying this gate twice is equivalent to applying the Y gate once,
/// and applying it after the √Y gate results in the identity gate.
pub struct InvSqrtPauliY {
    target: QuantumBit,
}

impl InvSqrtPauliY {
    pub fn new(target: &QuantumBit) -> Self {
        Self {
            target: target.clone(),
        }
    }
}

impl QuantumGate for InvSqrtPauliY {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = (-1.0 * PI) / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(-1.0 * s, 0.0)],
            [Complex::new(s, 0.0), Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        format!("Inv-Sqrt-Y(target={})", self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::InvSqrtPauliY(Self::new(&self.target)))
    }
}

#[derive(Debug, Clone)]
/// Represents the Y^t gate (Y to the power of t).
///
/// Y^t = [ [ cos(πt/2) - i*sin(πt/2), -sin(πt/2) - i*cos(πt/2) ],
///         [ sin(πt/2) - i*cos(πt/2),  cos(πt/2) + i*sin(πt/2) ] ]
///
/// When t=1, this reduces to the standard Pauli-Y gate.
pub struct YPowGate {
    target: QuantumBit,
    exponent: f64, // t parameter
}

impl YPowGate {
    pub fn new(target: &QuantumBit, exponent: f64) -> Self {
        Self {
            target: target.clone(),
            exponent,
        }
    }
}

impl QuantumGate for YPowGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let angle = PI * self.exponent / 2.0;
        let cos_angle = angle.cos();
        let sin_angle = angle.sin();

        array![
            [
                Complex::new(cos_angle, -sin_angle),
                Complex::new(-sin_angle, -cos_angle)
            ],
            [
                Complex::new(sin_angle, -cos_angle),
                Complex::new(cos_angle, sin_angle)
            ]
        ]
    }

    fn name(&self) -> String {
        format!("Y^{:.3}(target={})", self.exponent, self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::YPowGate(Self::new(
            &self.target,
            self.exponent,
        )))
    }
}
