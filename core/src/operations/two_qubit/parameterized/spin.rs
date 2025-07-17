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
/// Represents a general spin-spin interaction gate between two qubits.
///
/// This gate models the interaction between two quantum spins with coupling
/// strengths in all three spatial directions (x, y, z). It's commonly used
/// to simulate Heisenberg-type interactions in quantum many-body systems.
///
/// The 4x4 unitary matrix has a complex structure involving trigonometric
/// functions of the coupling parameters:
///
/// SI = [ [    cos(x-y)·e^(-iz),                 0,                 0, -i·sin(x-y)·e^(-iz) ],
///        [                   0,   cos(x+y)·e^(iz), i·sin(x+y)·e^(iz),                   0 ],
///        [                   0, i·sin(x+y)·e^(iz),   cos(x+y)·e^(iz),                   0 ],
///        [ -i·sin(x-y)·e^(-iz),                 0,                 0,    cos(x-y)·e^(-iz) ] ]
///
/// This gate is particularly useful for simulating quantum magnetism and
/// spin chain dynamics.
///
/// # Arguments
/// * `control` - The first qubit in the spin interaction
/// * `target` - The second qubit in the spin interaction
/// * `x` - Coupling strength in x-direction (σₓ ⊗ σₓ term)
/// * `y` - Coupling strength in y-direction (σᵧ ⊗ σᵧ term)
/// * `z` - Coupling strength in z-direction (σᵤ ⊗ σᵤ term)
///
/// The transformation is based on the Hamiltonian:
/// H = x·(σₓ ⊗ σₓ) + y·(σᵧ ⊗ σᵧ) + z·(σᵤ ⊗ σᵤ)
pub struct SpinInteraction {
    control: QuantumBit,
    target: QuantumBit,
    x: f64,
    y: f64,
    z: f64,
}

impl SpinInteraction {
    pub fn new(control: &QuantumBit, target: &QuantumBit, x: f64, y: f64, z: f64) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
            x,
            y,
            z,
        }
    }
}

impl QuantumGate for SpinInteraction {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let cm: f64 = (self.x - self.y).cos();
        let cp: f64 = (self.x + self.y).cos();
        let sm: f64 = (self.x - self.y).sin();
        let sp: f64 = (self.x + self.y).sin();
        let cz: f64 = self.z.cos();
        let sz: f64 = self.z.sin();
        array![
            [
                Complex::new(cm * cz, (-1.0) * cm * sz),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new((-1.0) * sm * sz, (-1.0) * sm * cz)
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(cp * cz, cp * sz),
                Complex::new(sp * sz, (-1.0) * sp * cz),
                Complex::new(0.0, 0.0)
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(sp * sz, (-1.0) * sp * cz),
                Complex::new(cp * cz, cp * sz),
                Complex::new(0.0, 0.0)
            ],
            [
                Complex::new((-1.0) * sm * sz, (-1.0) * sm * cz),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(cm * cz, (-1.0) * cm * sz)
            ],
        ]
    }

    fn name(&self) -> String {
        format!(
            "SI(control={}, target={}, x={}, y={}, z={})",
            self.control, self.target, self.x, self.y, self.z
        )
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::SpinInteraction(Self::new(
            &self.control,
            &self.target,
            self.x,
            self.y,
            self.z,
        )))
    }
}
