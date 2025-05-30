//! Copyright (c) 2024-2025 Quira, Inc.
//!
//! This file is part of Quira
//!
//! This program is free software: you can redistribute it and/or modify
//! it under the terms of the GNU Affero General Public License as published by
//! the Free Software Foundation, either version 3 of the License, or
//! (at your option) any later version.
//!
//! This program is distributed in the hope that it will be useful
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//! GNU Affero General Public License for more details.
//!
//! You should have received a copy of the GNU Affero General Public License
//! along with this program.  If not, see <http://www.gnu.org/licenses/>.

use crate::{
    operations::QuantumGate,
    types::{Complex, Matrix, Qubit},
};

use ndarray::array;

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
    control: Qubit,
    target: Qubit,
}

impl MolmerSorensenXX {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for MolmerSorensenXX {
    /// construct the 4x4 matrix representing the gate.
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

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("MSXX(control={}, target={})", self.control, self.target)
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

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
    /// construct the 4x4 matrix representing the gate.
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

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "Qsim(control={}, target={}, x={}, y={}, z={})",
            self.control, self.target, self.x, self.y, self.z
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

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
    control: Qubit,
    target: Qubit,
    hs: f64,
    is: f64,
    delta: f64,
}

impl Fsim {
    pub fn new(control: Qubit, target: Qubit, hs: f64, is: f64, delta: f64) -> Self {
        Self {
            control,
            target,
            hs,
            is,
            delta,
        }
    }
}

impl QuantumGate for Fsim {
    /// construct the 4x4 matrix representing the gate.
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

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "Fsim(control={}, target={}, hs={}, is={}, delta={})",
            self.control, self.target, self.hs, self.is, self.delta
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}
