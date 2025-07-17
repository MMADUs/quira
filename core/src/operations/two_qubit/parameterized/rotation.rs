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
/// Represents a Givens rotation gate.
///
/// This gate implements a Givens rotation in the computational basis,
/// which rotates the |01⟩ and |10⟩ subspace by angle θ with an additional
/// phase parameter φ. It's useful for quantum chemistry and optimization algorithms.
///
/// The matrix form is:
/// GR = [ [ 1,              0,      0,      0 ],
///        [ 0,  cos(θ)·e^(iφ), sin(θ),      0 ],
///        [ 0, -sin(θ)·e^(iφ), cos(θ),      0 ],
///        [ 0,              0,      0, e^(iφ) ] ]
///
/// # Parameters
/// - `control`: The first qubit index
/// - `target`: The second qubit index
/// - `theta`: The rotation angle θ in radians
pub struct GivensRotation {
    control: QuantumBit,
    target: QuantumBit,
    theta: f64,
    phi: f64,
}

impl GivensRotation {
    pub fn new(control: &QuantumBit, target: &QuantumBit, theta: f64, phi: f64) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
            theta,
            phi,
        }
    }
}

impl QuantumGate for GivensRotation {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let ct: f64 = (self.theta).cos();
        let st: f64 = (self.theta).sin();
        let cp: f64 = (self.phi).cos();
        let sp: f64 = (self.phi).sin();
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(ct * cp, ct * sp),
                Complex::new(st, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(-1.0 * st * cp, -1.0 * st * sp),
                Complex::new(ct, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(cp, sp),
            ]
        ]
    }

    fn name(&self) -> String {
        format!(
            "GR(control={}, target={}, theta={}, phi={})",
            self.control, self.target, self.theta, self.phi
        )
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::GivensRotation(Self::new(
            &self.control,
            &self.target,
            self.theta,
            self.phi,
        )))
    }
}

#[derive(Debug, Clone)]
/// Represents a Givens rotation gate with little-endian qubit ordering.
///
/// This is similar to the standard Givens rotation but uses little-endian
/// qubit ordering conventions. The rotation acts on the |01⟩ and |10⟩ subspace
/// with different matrix structure to account for the bit ordering.
///
/// The matrix form is:
/// GRLE = [ [ 1,              0,             0,      0 ],
///          [ 0,         cos(θ),        sin(θ),      0 ],
///          [ 0, -sin(θ)·e^(iφ), cos(θ)·e^(iφ),      0 ],
///          [ 0,              0,             0, e^(iφ) ] ]
///
/// # Parameters
/// - `control`: The first qubit index
/// - `target`: The second qubit index
/// - `theta`: The rotation angle θ in radians
/// - `phi`: The phase parameter φ in radians
pub struct GivensRotationLittleEndian {
    control: QuantumBit,
    target: QuantumBit,
    theta: f64,
    phi: f64,
}

impl GivensRotationLittleEndian {
    pub fn new(control: &QuantumBit, target: &QuantumBit, theta: f64, phi: f64) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
            theta,
            phi,
        }
    }
}

impl QuantumGate for GivensRotationLittleEndian {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let ct: f64 = (self.theta).cos();
        let st: f64 = (self.theta).sin();
        let cp: f64 = (self.phi).cos();
        let sp: f64 = (self.phi).sin();
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(ct, 0.0),
                Complex::new(st, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(-1.0 * st * cp, -1.0 * st * sp),
                Complex::new(ct * cp, ct * sp),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(cp, sp),
            ]
        ]
    }

    fn name(&self) -> String {
        format!(
            "GRLE(control={}, target={}, theta={}, phi={})",
            self.control, self.target, self.theta, self.phi
        )
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::GivensRotationLittleEndian(Self::new(
            &self.control,
            &self.target,
            self.theta,
            self.phi,
        )))
    }
}
