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
/// Represents the Pauli-X gate (also known as the bit-flip or NOT gate).
///
/// This gate flips the state of a qubit from |0⟩ to |1⟩ and vice versa.
/// It is the quantum equivalent of the classical NOT gate.
///
/// The matrix form is:
///
/// X = [ [ 0, 1 ],
///       [ 1, 0 ] ]
pub struct PauliX {
    target: QuantumBit,
}

impl PauliX {
    pub fn new(target: &QuantumBit) -> Self {
        Self {
            target: target.clone(),
        }
    }

    pub fn sqrt(&self) -> SqrtPauliX {
        SqrtPauliX::new(&self.target)
    }

    pub fn inv_sqrt(&self) -> InvSqrtPauliX {
        InvSqrtPauliX::new(&self.target)
    }

    pub fn power(&self, exponent: f64) -> XPowGate {
        XPowGate::new(&self.target, exponent)
    }
}

impl QuantumGate for PauliX {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)]
        ]
    }

    fn name(&self) -> String {
        format!("X(target={})", self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::PauliX(Self::new(&self.target)))
    }
}

#[derive(Debug, Clone)]
/// Represents the square root of the Pauli-X gate.
///
/// This gate performs half the action of the Pauli-X gate, equivalent to
/// a 90-degree rotation around the X-axis.
///
/// The matrix form is:
///
/// √X = [ [ (1+i)/2, (1-i)/2 ],
///        [ (1-i)/2, (1+i)/2 ] ]
///
/// Applying this gate twice is equivalent to applying the X gate once.
pub struct SqrtPauliX {
    target: QuantumBit,
}

impl SqrtPauliX {
    pub fn new(target: &QuantumBit) -> Self {
        Self {
            target: target.clone(),
        }
    }
}

impl QuantumGate for SqrtPauliX {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = PI / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(0.0, -1.0 * s)],
            [Complex::new(0.0, -1.0 * s), Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        format!("Sqrt-X(target={})", self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::SqrtPauliX(Self::new(&self.target)))
    }
}

#[derive(Debug, Clone)]
/// Represents the inverse of the square root of the Pauli-X gate.
///
/// This gate is the adjoint (inverse) of the √X gate, equivalent to
/// a -90-degree rotation around the X-axis.
///
/// The matrix form is the complex conjugate transpose of the √X matrix.
///
/// Applying this gate twice is equivalent to applying the X gate once,
/// and applying it after the √X gate results in the identity gate.
pub struct InvSqrtPauliX {
    target: QuantumBit,
}

impl InvSqrtPauliX {
    pub fn new(target: &QuantumBit) -> Self {
        Self {
            target: target.clone(),
        }
    }
}

impl QuantumGate for InvSqrtPauliX {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = PI / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(0.0, 1.0 * s)],
            [Complex::new(0.0, 1.0 * s), Complex::new(c, 0.0)]
        ]
    }

    fn name(&self) -> String {
        format!("Inv-Sqrt-X(target={})", self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::InvSqrtPauliX(Self::new(&self.target)))
    }
}

#[derive(Debug, Clone)]
/// Represents the X^t gate (X to the power of t).
///
/// X^t = [ [ cos(πt/2) - i*sin(πt/2), -i*sin(πt/2) ],
///         [ -i*sin(πt/2),           cos(πt/2) + i*sin(πt/2) ] ]
///
/// When t=1, this reduces to the standard Pauli-X gate.
pub struct XPowGate {
    target: QuantumBit,
    exponent: f64, // t parameter
}

impl XPowGate {
    pub fn new(target: &QuantumBit, exponent: f64) -> Self {
        Self {
            target: target.clone(),
            exponent,
        }
    }

    pub fn phased(&self, phase: f64) -> PhasedXPowGate {
        PhasedXPowGate::new(&self.target, self.exponent, phase)
    }
}

impl QuantumGate for XPowGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let angle = PI * self.exponent / 2.0;
        let cos_angle = angle.cos();
        let sin_angle = angle.sin();

        array![
            [
                Complex::new(cos_angle, -sin_angle),
                Complex::new(0.0, -sin_angle)
            ],
            [
                Complex::new(0.0, -sin_angle),
                Complex::new(cos_angle, sin_angle)
            ]
        ]
    }

    fn name(&self) -> String {
        format!("X^{:.3}(target={})", self.exponent, self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::XPowGate(Self::new(
            &self.target,
            self.exponent,
        )))
    }
}

#[derive(Debug, Clone)]
/// Represents the PhasedXPow gate.
///
/// This gate applies an X^t rotation with an additional phase offset.
/// PhasedXPow(t, φ) = Z^(-φ/π) * X^t * Z^(φ/π)
///
/// The matrix form is:
/// [ [ cos(πt/2) - i*sin(πt/2)*cos(φ), -i*sin(πt/2)*sin(φ) - sin(πt/2)*cos(φ) ],
///   [ -i*sin(πt/2)*sin(φ) + sin(πt/2)*cos(φ), cos(πt/2) + i*sin(πt/2)*cos(φ) ] ]
pub struct PhasedXPowGate {
    target: QuantumBit,
    exponent: f64,       // t parameter
    phase_exponent: f64, // φ parameter
}

impl PhasedXPowGate {
    pub fn new(target: &QuantumBit, exponent: f64, phase_exponent: f64) -> Self {
        Self {
            target: target.clone(),
            exponent,
            phase_exponent,
        }
    }
}

impl QuantumGate for PhasedXPowGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let t_angle = PI * self.exponent / 2.0;
        let cos_t = t_angle.cos();
        let sin_t = t_angle.sin();
        let cos_phi = self.phase_exponent.cos();
        let sin_phi = self.phase_exponent.sin();

        array![
            [
                Complex::new(cos_t, -sin_t * cos_phi),
                Complex::new(-sin_t * cos_phi, -sin_t * sin_phi)
            ],
            [
                Complex::new(sin_t * cos_phi, -sin_t * sin_phi),
                Complex::new(cos_t, sin_t * cos_phi)
            ]
        ]
    }

    fn name(&self) -> String {
        format!(
            "PhasedX^{:.3}(target={}, φ={:.3})",
            self.exponent, self.target, self.phase_exponent
        )
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::PhasedXPowGate(Self::new(
            &self.target,
            self.exponent,
            self.phase_exponent,
        )))
    }
}
