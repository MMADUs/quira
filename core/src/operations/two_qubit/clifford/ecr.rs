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
/// Represents the Echo Cross Resonance (ECR) gate, a native gate in IBM quantum processors.
///
/// The ECR gate is a two-qubit entangling gate that is natively implemented
/// in IBM's superconducting quantum processors. It creates a specific pattern
/// of entanglement and is optimized for the physical constraints of the hardware.
///
/// The gate creates a complex superposition mixing all four basis states
/// with specific amplitude and phase relationships. It's particularly useful
/// for implementing quantum algorithms on IBM quantum hardware.
///
/// The 4x4 unitary matrix is:
///
/// ECR = (1/√2) × [ [ 0,  0,  1,  i ],
///                  [ 0,  0,  i,  1 ],
///                  [ 1, -i,  0,  0 ],
///                  [-i,  1,  0,  0 ] ]
///
/// # Arguments
/// * `control` - The first qubit in the ECR operation
/// * `target` - The second qubit in the ECR operation
pub struct EchoCrossResonance {
    control: QuantumBit,
    target: QuantumBit,
}

impl EchoCrossResonance {
    pub fn new(control: &QuantumBit, target: &QuantumBit) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
        }
    }
}

impl QuantumGate for EchoCrossResonance {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let matrix: Matrix<Complex> = array![
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 1.0)
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 1.0),
                Complex::new(1.0, 0.0)
            ],
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, -1.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0)
            ],
            [
                Complex::new(0.0, -1.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0)
            ]
        ];
        matrix / 2.0_f64.sqrt()
    }

    fn name(&self) -> String {
        format!("ECR(control={}, target={})", self.control, self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::EchoCrossResonance(Self::new(
            &self.control,
            &self.target,
        )))
    }
}
