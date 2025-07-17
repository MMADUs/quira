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
/// Represents the SWAP gate, which exchanges the states of two qubits.
///
/// The SWAP gate interchanges the quantum states of two qubits completely.
/// It's equivalent to three CNOT gates and is used for qubit routing and
/// quantum circuit optimization.
///
/// The gate transforms the computational basis states as:
/// - |00⟩ → |00⟩ (unchanged)
/// - |01⟩ → |10⟩ (qubits swapped)
/// - |10⟩ → |01⟩ (qubits swapped)
/// - |11⟩ → |11⟩ (unchanged)
///
/// The 4x4 unitary matrix is:
///
/// SWAP = [ [ 1, 0, 0, 0 ],
///          [ 0, 0, 1, 0 ],
///          [ 0, 1, 0, 0 ],
///          [ 0, 0, 0, 1 ] ]
///
/// # Parameters
/// * `control` - The first qubit to be swapped
/// * `target` - The second qubit to be swapped
pub struct SWAP {
    control: QuantumBit,
    target: QuantumBit,
}

impl SWAP {
    pub fn new(control: &QuantumBit, target: &QuantumBit) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
        }
    }
}

impl QuantumGate for SWAP {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
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
        format!("SWAP(control={}, target={})", self.control, self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::SWAP(Self::new(&self.control, &self.target)))
    }
}

#[derive(Debug, Clone)]
/// Represents the iSWAP gate, which swaps qubits with an additional phase factor.
///
/// The iSWAP gate exchanges the |01⟩ and |10⟩ states while introducing an
/// imaginary phase factor of i. The iSWAP gate can be decomposed into SWAP followed by controlled-Z gates.
///
/// The gate transforms the computational basis states as:
/// - |00⟩ → |00⟩ (unchanged)
/// - |01⟩ → i|10⟩ (swapped with phase i)
/// - |10⟩ → i|01⟩ (swapped with phase i)
/// - |11⟩ → |11⟩ (unchanged)
///
/// The 4x4 unitary matrix is:
///
/// iSWAP = [ [ 1, 0, 0, 0 ],
///           [ 0, 0, i, 0 ],
///           [ 0, i, 0, 0 ],
///           [ 0, 0, 0, 1 ] ]
///
/// # Parameters
/// * `control` - The first qubit to be swapped with phase
/// * `target` - The second qubit to be swapped with phase
pub struct ISWAP {
    control: QuantumBit,
    target: QuantumBit,
}

impl ISWAP {
    pub fn new(control: &QuantumBit, target: &QuantumBit) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
        }
    }
}

impl QuantumGate for ISWAP {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 1.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 1.0),
                Complex::new(0.0, 0.0),
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
        format!("ISWAP(control={}, target={})", self.control, self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::ISWAP(Self::new(&self.control, &self.target)))
    }
}

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
    control: QuantumBit,
    target: QuantumBit,
}

impl SqrtISWAP {
    pub fn new(control: &QuantumBit, target: &QuantumBit) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
        }
    }
}

impl QuantumGate for SqrtISWAP {
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

    fn name(&self) -> String {
        format!(
            "Sqrt-ISWAP(control={}, target={})",
            self.control, self.target
        )
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::SqrtISWAP(Self::new(
            &self.control,
            &self.target,
        )))
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
    control: QuantumBit,
    target: QuantumBit,
}

impl InvSqrtISWAP {
    pub fn new(control: &QuantumBit, target: &QuantumBit) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
        }
    }
}

impl QuantumGate for InvSqrtISWAP {
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

    fn name(&self) -> String {
        format!(
            "Inv-Sqrt-ISWAP(control={}, target={})",
            self.control, self.target
        )
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::InvSqrtISWAP(Self::new(
            &self.control,
            &self.target,
        )))
    }
}

#[derive(Debug, Clone)]
/// Represents the Fredkin gate (FSWAP), a controlled SWAP gate.
///
/// The FSWAP gate swaps two qubits while introducing a phase factor of -1
/// when both qubits are in state |11⟩. This creates a conditional swap
/// operation with phase correction.
///
/// The gate transforms the computational basis states as:
/// - |00⟩ → |00⟩ (unchanged)
/// - |01⟩ → |10⟩ (swapped)
/// - |10⟩ → |01⟩ (swapped)
/// - |11⟩ → -|11⟩ (phase flip)
///
/// The 4x4 unitary matrix is:
///
/// FSWAP = [ [ 1,  0, 0,  0 ],
///           [ 0,  0, 1,  0 ],
///           [ 0,  1, 0,  0 ],
///           [ 0,  0, 0, -1 ] ]
///
/// # Parameters
/// * `control` - The first qubit in the conditional swap
/// * `target` - The second qubit in the conditional swap
pub struct FSWAP {
    control: QuantumBit,
    target: QuantumBit,
}

impl FSWAP {
    pub fn new(control: &QuantumBit, target: &QuantumBit) -> Self {
        Self {
            control: control.clone(),
            target: target.clone(),
        }
    }
}

impl QuantumGate for FSWAP {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(-1.0, 0.0)
            ]
        ]
    }

    fn name(&self) -> String {
        format!("FSWAP(control={}, target={})", self.control, self.target)
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index(), self.control.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::FSWAP(Self::new(&self.control, &self.target)))
    }
}
