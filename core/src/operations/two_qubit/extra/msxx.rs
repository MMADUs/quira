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
/// Represents the Mølmer-Sørensen XX gate, a two-qubit entangling gate.
///
/// This gate implements the XX interaction between two qubits and is particularly
/// useful in trapped-ion quantum computing. It creates entanglement between qubits
/// by applying a rotation around the XX axis in the Bloch sphere representation.
///
/// The gate acts on the computational basis states as:
/// - |00⟩ → (|00⟩ - i|11⟩)/√2
/// - |01⟩ → (|01⟩ - i|10⟩)/√2  
/// - |10⟩ → (|10⟩ - i|01⟩)/√2
/// - |11⟩ → (|11⟩ - i|00⟩)/√2
///
/// The matrix form is:
///
/// MSXX = (1/√2) * [ [ 1,  0,  0, -i ],
///                   [ 0,  1, -i,  0 ],
///                   [ 0, -i,  1,  0 ],
///                   [-i,  0,  0,  1 ] ]
///
/// This gate is self-inverse when applied twice with appropriate phase corrections.
///
/// # Arguments
/// * `control` - The control qubit
/// * `target` - The target qubit
pub struct MolmerSorensenXX {
    control: QuantumBit,
    target: QuantumBit,
}

impl MolmerSorensenXX {
    pub fn new(control: &QuantumBit, target: &QuantumBit) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
        }
    }
}

impl QuantumGate for MolmerSorensenXX {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let f: f64 = 1.0 / ((2.0_f64).sqrt());
        array![
            [
                Complex::new(f, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -1.0 * f),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(f, 0.0),
                Complex::new(0.0, -1.0 * f),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -1.0 * f),
                Complex::new(f, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, -1.0 * f),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(f, 0.0),
            ]
        ]
    }

    fn name(&self) -> String {
        format!("MSXX(control={}, target={})", self.control, self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::MolmerSorensenXX(Self::new(
            &self.control,
            &self.target,
        )))
    }
}
