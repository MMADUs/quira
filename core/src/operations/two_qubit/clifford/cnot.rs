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
/// Represents the Controlled-NOT (CNOT) gate, one of the most fundamental two-qubit gates.
///
/// The CNOT gate flips the target qubit if and only if the control qubit is in state |1⟩.
/// It's a key building block for quantum computation and can be used to create entanglement.
///
/// The gate transforms the computational basis states as:
/// - |00⟩ → |00⟩ (control |0⟩: target unchanged)
/// - |01⟩ → |01⟩ (control |0⟩: target unchanged)  
/// - |10⟩ → |11⟩ (control |1⟩: target flipped)
/// - |11⟩ → |10⟩ (control |1⟩: target flipped)
///
/// The 4x4 unitary matrix is:
///
/// CNOT = [ [ 1, 0, 0, 0 ],
///          [ 0, 1, 0, 0 ],
///          [ 0, 0, 0, 1 ],
///          [ 0, 0, 1, 0 ] ]
///
/// # Parameters
/// * `control` - The control qubit that determines whether the operation is applied
/// * `target` - The target qubit that gets flipped when control is |1⟩
pub struct ControlledNot {
    control: Qubit,
    target: Qubit,
}

impl ControlledNot {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for ControlledNot {
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
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
            ]
        ]
    }

    fn name(&self) -> String {
        format!("CNOT(control={}, target={})", self.control, self.target)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::ControlledNot(Self::new(
            self.control,
            self.target,
        )))
    }
}
