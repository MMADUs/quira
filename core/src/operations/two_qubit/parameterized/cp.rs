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
use crate::operations::two_qubit::TwoQubitType;
use crate::{Complex, GateType, Matrix, Qubit};

#[derive(Debug, Clone)]
/// Represents a controlled phase shift gate.
///
/// This gate applies a phase shift e^(iθ) to the |11⟩ state when both
/// control and target qubits are in the |1⟩ state, leaving other states unchanged.
///
/// The matrix form is:
/// CPS = [ [ 1, 0, 0,      0 ],
///         [ 0, 1, 0,      0 ],
///         [ 0, 0, 1,      0 ],
///         [ 0, 0, 0, e^(iθ) ] ]
///
/// # Parameters
/// - `control`: The control qubit index
/// - `target`: The target qubit index
/// - `theta`: The phase shift angle θ in radians
pub struct ControlledPhaseShift {
    control: Qubit,
    target: Qubit,
    theta: f64,
}

impl ControlledPhaseShift {
    pub fn new(control: Qubit, target: Qubit, theta: f64) -> Self {
        Self {
            control,
            target,
            theta,
        }
    }
}

impl QuantumGate for ControlledPhaseShift {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta).cos();
        let s: f64 = (self.theta).sin();
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
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(c, s),
            ]
        ]
    }

    fn name(&self) -> String {
        format!(
            "CPS(control={}, target={}, theta={})",
            self.control, self.target, self.theta
        )
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::ControlledPhaseShift(Self::new(
            self.control,
            self.target,
            self.theta,
        )))
    }
}

#[derive(Debug, Clone)]
/// Represents a phase-shifted controlled phase gate.
///
/// This gate combines a phase shift on intermediate states with a controlled
/// phase operation. It applies phase shifts e^(iφ) on the |01⟩ and |10⟩ states,
/// and e^(i(2φ+θ)) on the |11⟩ state.
///
/// The matrix form is:
/// PSCP = [ [ 1,      0,      0,           0 ],
///          [ 0, e^(iφ),      0,           0 ],
///          [ 0,      0, e^(iφ),           0 ],
///          [ 0,      0,      0, e^(i(2φ+θ)) ] ]
///
/// # Parameters
/// - `control`: The control qubit index
/// - `target`: The target qubit index
/// - `theta`: The additional phase parameter θ
/// - `phi`: The base phase shift parameter φ
pub struct PhaseShiftedControlledPhase {
    control: Qubit,
    target: Qubit,
    theta: f64,
    phi: f64,
}

impl PhaseShiftedControlledPhase {
    pub fn new(control: Qubit, target: Qubit, theta: f64, phi: f64) -> Self {
        Self {
            control,
            target,
            theta,
            phi,
        }
    }
}

impl QuantumGate for PhaseShiftedControlledPhase {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let phi: f64 = self.phi;
        let theta: f64 = self.theta;
        let c: f64 = phi.cos();
        let s: f64 = phi.sin();
        let cos: f64 = (2.0 * phi + theta).cos();
        let sin: f64 = (2.0 * phi + theta).sin();
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(c, s),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(c, s),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(cos, sin),
            ]
        ]
    }

    fn name(&self) -> String {
        format!(
            "PSCP(control={}, target={}, theta={}, phi={})",
            self.control, self.target, self.theta, self.phi
        )
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::PhaseShiftedControlledPhase(Self::new(
            self.control,
            self.target,
            self.theta,
            self.phi,
        )))
    }
}
