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
/// Represents an XY gate (iSWAP-like interaction).
///
/// This gate implements an XY interaction between two qubits, creating
/// entanglement by mixing the |01⟩ and |10⟩ states with a rotation angle θ.
/// It's related to the iSWAP gate family.
///
/// The matrix form is:
/// XY = [ [ 1,          0,          0, 0 ],
///        [ 0,   cos(θ/2), i·sin(θ/2), 0 ],
///        [ 0, i·sin(θ/2),   cos(θ/2), 0 ],
///        [ 0,          0,          0, 1 ] ]
///
/// # Parameters
/// - `control`: The first qubit index
/// - `target`: The second qubit index
/// - `theta`: The interaction strength θ in radians
pub struct XY {
    control: Qubit,
    target: Qubit,
    theta: f64,
}

impl XY {
    pub fn new(control: Qubit, target: Qubit, theta: f64) -> Self {
        Self {
            control,
            target,
            theta,
        }
    }
}

impl QuantumGate for XY {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(c, 0.0),
                Complex::new(0.0, s),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, s),
                Complex::new(c, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ]
        ]
    }

    fn name(&self) -> String {
        format!(
            "XY(control={}, target={}, theta={})",
            self.control, self.target, self.theta
        )
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::XY(Self::new(
            self.control,
            self.target,
            self.theta,
        )))
    }
}
