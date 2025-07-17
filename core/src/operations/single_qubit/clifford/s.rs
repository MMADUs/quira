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
/// Represents the S gate (also known as the phase gate or Z90 gate).
///
/// This gate leaves the |0⟩ state unchanged and maps |1⟩ to i|1⟩,
/// introducing a π/2 phase shift to the |1⟩ state.
///
/// The matrix form is:
///
/// S = [ [ 1, 0 ],
///       [ 0, i ] ]
///
/// Applying this gate twice is equivalent to applying the Z gate once.
pub struct SGate {
    target: QuantumBit,
}

impl SGate {
    pub fn new(target: &QuantumBit) -> Self {
        Self {
            target: target.clone(),
        }
    }

    pub fn dagger(&self) -> InvSGate {
        InvSGate::new(&self.target)
    }
}

impl QuantumGate for SGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), Complex::new(0.0, 1.0)]
        ]
    }

    fn name(&self) -> String {
        format!("S(target={})", self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::SGate(Self::new(&self.target)))
    }
}

#[derive(Debug, Clone)]
/// Represents the inverse of the S gate.
///
/// This gate leaves the |0⟩ state unchanged and maps |1⟩ to -i|1⟩,
/// introducing a -π/2 phase shift to the |1⟩ state.
///
/// The matrix form is:
///
/// S† = [ [ 1,  0 ],
///        [ 0, -i ] ]
///
/// Applying this gate after the S gate results in the identity gate.
pub struct InvSGate {
    target: QuantumBit,
}

impl InvSGate {
    pub fn new(target: &QuantumBit) -> Self {
        Self {
            target: target.clone(),
        }
    }
}

impl QuantumGate for InvSGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), Complex::new(0.0, -1.0)]
        ]
    }

    fn name(&self) -> String {
        format!("Inv-S(target={})", self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::InvSGate(Self::new(&self.target)))
    }
}
