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
/// Represents the square root of iSWAP gate, a two-qubit entangling gate.
///
/// The √iSWAP gate is the square root of the iSWAP gate, meaning that applying
/// it twice produces the full iSWAP operation. It performs a partial swap of
/// two qubits with an additional phase factor.
///
/// This gate transforms the computational basis states as:
/// - |00⟩ → |00⟩ (unchanged)
/// - |01⟩ → (|01⟩ + i|10⟩)/√2 (creates superposition with phase)
/// - |10⟩ → (i|01⟩ + |10⟩)/√2 (creates superposition with phase)
/// - |11⟩ → |11⟩ (unchanged)
///
/// The 4x4 unitary matrix is:
///
/// √iSWAP = [ [ 1,    0,    0, 0 ],
///            [ 0, 1/√2, i/√2, 0 ],
///            [ 0, i/√2, 1/√2, 0 ],
///            [ 0,    0,    0, 1 ] ]
///
/// # Arguments
/// * `control` - The first qubit to be partially swapped
/// * `target` - The second qubit to be partially swapped
pub struct SqrtISWAP {
    control: Qubit,
    target: Qubit,
}

impl SqrtISWAP {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for SqrtISWAP {
    /// construct the 4x4 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let f: f64 = 1.0 / ((2.0_f64).sqrt());
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(f, 0.0),
                Complex::new(0.0, f),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, f),
                Complex::new(f, 0.0),
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

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "Sqrt-ISWAP(control={}, target={})",
            self.control, self.target
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
/// Represents the inverse square root of iSWAP gate.
///
/// The inverse √iSWAP gate (or conjugate transpose of √iSWAP) undoes the
/// operation of the √iSWAP gate. It's the Hermitian conjugate of the √iSWAP
/// gate, obtained by taking the complex conjugate and transpose of the matrix.
///
/// This gate transforms the computational basis states as:
/// - |00⟩ → |00⟩ (unchanged)
/// - |01⟩ → (|01⟩ - i|10⟩)/√2 (creates superposition with opposite phase)
/// - |10⟩ → (-i|01⟩ + |10⟩)/√2 (creates superposition with opposite phase)
/// - |11⟩ → |11⟩ (unchanged)
///
/// The 4x4 unitary matrix is:
///
/// (√iSWAP)† = [ [ 1,     0,     0, 0 ],
///               [ 0,  1/√2, -i/√2, 0 ],
///               [ 0, -i/√2,  1/√2, 0 ],
///               [ 0,     0,     0, 1 ] ]
///
/// This gate satisfies the property: (√iSWAP)† × √iSWAP = I (identity).
///
/// # Arguments
/// * `control` - The first qubit to be operated on
/// * `target` - The second qubit to be operated on
pub struct InvSqrtISWAP {
    control: Qubit,
    target: Qubit,
}

impl InvSqrtISWAP {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for InvSqrtISWAP {
    /// construct the 4x4 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let f: f64 = 1.0 / ((2.0_f64).sqrt());
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
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
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "Inv-Sqrt-ISWAP(control={}, target={})",
            self.control, self.target
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}
