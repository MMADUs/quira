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
/// Represents the Qsim gate, a parameterized two-qubit gate used in quantum simulation.
///
/// The Qsim gate is a general two-qubit gate parameterized by three angles (x, y, z)
/// that allow for flexible quantum operations. It's particularly useful for simulating
/// arbitrary two-qubit interactions in quantum systems.
///
/// The matrix elements are constructed using:
/// - cm = cos(x - y), cp = cos(x + y)
/// - sm = sin(x - y), sp = sin(x + y)
/// - cz = cos(z), sz = sin(z)
///
/// The matrix form is:
///
/// Qsim = [ [  cm*e^(-iz),         0,         0, -sm*e^(-iz) ],
///          [           0, sp*e^(iz), cp*e^(iz),           0 ],
///          [           0, cp*e^(iz), sp*e^(iz),           0 ],
///          [ -sm*e^(-iz),         0,         0,  cm*e^(-iz) ] ]
///
/// This gate can simulate various two-qubit interactions by choosing appropriate parameters.
///
/// # Arguments
/// * `control` - The control qubit
/// * `target` - The target qubit
/// * `x` - First rotation parameter affecting the computational basis mixing (in radians)
/// * `y` - Second rotation parameter affecting the computational basis mixing (in radians)
/// * `z` - Phase rotation parameter (in radians)
pub struct Qsim {
    control: Qubit,
    target: Qubit,
    x: f64,
    y: f64,
    z: f64,
}

impl Qsim {
    pub fn new(control: Qubit, target: Qubit, x: f64, y: f64, z: f64) -> Self {
        Self {
            control,
            target,
            x,
            y,
            z,
        }
    }
}

impl QuantumGate for Qsim {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let cm: f64 = (self.x - self.y).cos();
        let cp: f64 = (self.x + self.y).cos();
        let sm: f64 = (self.x - self.y).sin();
        let sp: f64 = (self.x + self.y).sin();
        let cz: f64 = self.z.cos();
        let sz: f64 = self.z.sin();
        array![
            [
                Complex::new(cm * cz, -1.0 * cm * sz),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(-1.0 * sm * sz, -1.0 * sm * cz),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(sp * sz, -1.0 * sp * cz),
                Complex::new(cp * cz, cp * sz),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(cp * cz, cp * sz),
                Complex::new(sp * sz, -1.0 * sp * cz),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(-1.0 * sm * sz, -1.0 * sm * cz),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(cm * cz, -1.0 * cm * sz),
            ]
        ]
    }

    fn name(&self) -> String {
        format!(
            "Qsim(control={}, target={}, x={}, y={}, z={})",
            self.control, self.target, self.x, self.y, self.z
        )
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::Qsim(Self::new(
            self.control,
            self.target,
            self.x,
            self.y,
            self.z,
        )))
    }
}
