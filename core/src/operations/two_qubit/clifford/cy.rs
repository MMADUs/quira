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
/// Represents the Controlled-Y (CY) gate, applying Pauli-Y to the target when control is |1⟩.
///
/// The CY gate applies a Pauli-Y rotation to the target qubit if and only if
/// the control qubit is in state |1⟩. The Pauli-Y gate flips the qubit state
/// and introduces a phase factor of i.
///
/// The gate transforms the computational basis states as:
/// - |00⟩ → |00⟩ (control |0⟩: target unchanged)
/// - |01⟩ → |01⟩ (control |0⟩: target unchanged)
/// - |10⟩ → i|11⟩ (control |1⟩: Y applied to target)
/// - |11⟩ → -i|10⟩ (control |1⟩: Y applied to target)
///
/// The 4x4 unitary matrix is:
///
/// CY = [ [ 1,  0,  0,   0 ],
///        [ 0,  1,  0,  -i ],
///        [ 0,  0,  0,  -i ],
///        [ 0,  0,  i,   0 ] ]
///
/// # Arguments
/// * `control` - The control qubit that determines whether Y is applied
/// * `target` - The target qubit that receives the Pauli-Y operation
pub struct ControlledPauliY {
    control: QuantumBit,
    target: QuantumBit,
}

impl ControlledPauliY {
    pub fn new(control: &QuantumBit, target: &QuantumBit) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
        }
    }
}

impl QuantumGate for ControlledPauliY {
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
                Complex::new(0.0, -1.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -1.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 1.0),
                Complex::new(0.0, 0.0),
            ]
        ]
    }

    fn name(&self) -> String {
        format!("CY(control={}, target={})", self.control, self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::ControlledPauliY(Self::new(
            &self.control,
            &self.target,
        )))
    }
}
