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
/// Represents the Fsim (Fermionic simulation) gate, used for simulating fermionic systems.
///
/// The Fsim gate is designed to simulate fermionic interactions in quantum systems,
/// particularly useful in quantum chemistry and condensed matter physics simulations.
/// It implements a combination of iswap-like operations with controlled phase rotations.
///
/// The gate acts on the computational basis as:
/// - |00⟩ → cos(δ)|00⟩ + i*sin(δ)|11⟩
/// - |01⟩ → -i*sin(hs)|01⟩ + cos(hs)|10⟩
/// - |10⟩ → cos(hs)|01⟩ - i*sin(hs)|10⟩
/// - |11⟩ → complex combination involving all parameters
///
/// # Arguments
/// * `control` - The control qubit
/// * `target` - The target qubit
/// * `hs` - Hopping strength parameter controlling the iswap-like interaction (in radians)
/// * `is` - Interaction strength parameter for the phase rotation (in radians)
/// * `delta` - Phase parameter for additional control (in radians)
pub struct Fsim {
    control: QuantumBit,
    target: QuantumBit,
    hs: f64,
    is: f64,
    delta: f64,
}

impl Fsim {
    pub fn new(control: &QuantumBit, target: &QuantumBit, hs: f64, is: f64, delta: f64) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
            hs,
            is,
            delta,
        }
    }
}

impl QuantumGate for Fsim {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let h: f64 = self.hs;
        let i: f64 = self.is;
        let d: f64 = self.delta;
        array![
            [
                Complex::new(d.cos(), 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, d.sin()),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -1.0 * h.sin()),
                Complex::new(h.cos(), 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(h.cos(), 0.0),
                Complex::new(0.0, -1.0 * h.sin()),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(-1.0 * d.sin() * i.sin(), -1.0 * d.sin() * i.cos()),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(-1.0 * d.cos() * i.cos(), d.cos() * i.sin())
            ]
        ]
    }

    fn name(&self) -> String {
        format!(
            "Fsim(control={}, target={}, hs={}, is={}, delta={})",
            self.control, self.target, self.hs, self.is, self.delta
        )
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::Fsim(Self::new(
            &self.control,
            &self.target,
            self.hs,
            self.is,
            self.delta,
        )))
    }
}
