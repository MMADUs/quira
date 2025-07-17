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
use crate::operations::two_qubit::TwoQubitType;
use crate::ops::QuantumGate;
use crate::types::{Complex, Matrix};

#[derive(Debug, Clone)]
/// Represents the Controlled-Z (CZ) gate, applying a phase flip when both qubits are |1⟩.
///
/// The CZ gate applies a phase factor of -1 to the |11⟩ state while leaving
/// all other computational basis states unchanged. This gate is symmetric
/// with respect to both qubits and is commonly used for creating entanglement.
///
/// The gate transforms the computational basis states as:
/// - |00⟩ → |00⟩ (unchanged)
/// - |01⟩ → |01⟩ (unchanged)
/// - |10⟩ → |10⟩ (unchanged)
/// - |11⟩ → -|11⟩ (phase flip)
///
/// The 4x4 unitary matrix is:
///
/// CZ = [ [ 1,  0, 0,  0 ],
///        [ 0,  1, 0,  0 ],
///        [ 0,  0, 1,  0 ],
///        [ 0,  0, 0, -1 ] ]
///
/// The CZ gate is self-inverse and commutes with both qubits. It's equivalent
/// to a CNOT gate with Hadamard gates applied to the target qubit.
///
/// # Parameters
/// * `control` - The first qubit in the controlled phase operation
/// * `target` - The second qubit in the controlled phase operation
pub struct ControlledPauliZ {
    control: QuantumBit,
    target: QuantumBit,
}

impl ControlledPauliZ {
    pub fn new(control: &QuantumBit, target: &QuantumBit) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
        }
    }

    pub fn phase_shift(&self, phi: f64) -> PhaseShiftedControlledZ {
        PhaseShiftedControlledZ::new(&self.control, &self.target, phi)
    }
}

impl QuantumGate for ControlledPauliZ {
    fn unitary_matrix(&self) -> Matrix<Complex> {
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
                Complex::new(-1.0, 0.0)
            ]
        ]
    }

    fn name(&self) -> String {
        format!("CZ(control={}, target={})", self.control, self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::ControlledPauliZ(Self::new(
            &self.control,
            &self.target,
        )))
    }
}

#[derive(Debug, Clone)]
/// Represents a phase-shifted controlled-Z gate.
///
/// This gate applies a phase shift to the control qubit and then performs
/// a controlled-Z operation. The phase shift is applied as e^(iφ) when the
/// control qubit is in the |1⟩ state.
///
/// The matrix form is:
/// PSCZ = [ [ 1,      0,      0,           0 ],
///          [ 0, e^(iφ),      0,           0 ],
///          [ 0,      0, e^(iφ),           0 ],
///          [ 0,      0,      0, e^(i(2φ+π)) ] ]
///
/// # Parameters
/// - `control`: The control qubit index
/// - `target`: The target qubit index  
/// - `phi`: The phase shift parameter φ
pub struct PhaseShiftedControlledZ {
    control: QuantumBit,
    target: QuantumBit,
    phi: f64,
}

impl PhaseShiftedControlledZ {
    pub fn new(control: &QuantumBit, target: &QuantumBit, phi: f64) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
            phi,
        }
    }
}

impl QuantumGate for PhaseShiftedControlledZ {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let phi: f64 = self.phi;
        let c: f64 = phi.cos();
        let s: f64 = phi.sin();
        let cos: f64 = (2.0 * phi + PI).cos();
        let sin: f64 = (2.0 * phi + PI).sin();
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
            "PSCZ(control={}, target={}, phi={:.4})",
            self.control, self.target, self.phi
        )
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::PhaseShiftedControlledZ(Self::new(
            &self.control,
            &self.target,
            self.phi,
        )))
    }
}
