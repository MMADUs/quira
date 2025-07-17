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
/// Represents a variable Mølmer-Sørensen XX gate.
///
/// This gate implements a parametric version of the Mølmer-Sørensen gate,
/// which creates entanglement through XX interactions. It's commonly used
/// in trapped ion quantum computing.
///
/// The matrix form is:
/// VMSXX = [ [    cos(θ/2),           0,           0, -i·sin(θ/2) ],
///           [           0,    cos(θ/2), -i·sin(θ/2),           0 ],
///           [           0, -i·sin(θ/2),    cos(θ/2),           0 ],
///           [ -i·sin(θ/2),           0,           0,    cos(θ/2) ] ]
///
/// # Parameters
/// - `control`: The first qubit index
/// - `target`: The second qubit index
/// - `theta`: The interaction strength θ in radians
pub struct VariableMSXX {
    control: QuantumBit,
    target: QuantumBit,
    theta: f64,
}

impl VariableMSXX {
    pub fn new(control: &QuantumBit, target: &QuantumBit, theta: f64) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
            theta,
        }
    }
}

impl QuantumGate for VariableMSXX {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        array![
            [
                Complex::new(c, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -1.0 * s),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(c, 0.0),
                Complex::new(0.0, -1.0 * s),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -1.0 * s),
                Complex::new(c, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, -1.0 * s),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(c, 0.0),
            ]
        ]
    }

    fn name(&self) -> String {
        format!(
            "VMSXX(control={}, target={}, theta={})",
            self.control, self.target, self.theta
        )
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::VariableMSXX(Self::new(
            &self.control,
            &self.target,
            self.theta,
        )))
    }
}
