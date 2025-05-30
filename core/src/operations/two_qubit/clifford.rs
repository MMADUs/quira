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
/// Represents the Controlled-NOT (CNOT) gate, one of the most fundamental two-qubit gates.
///
/// The CNOT gate flips the target qubit if and only if the control qubit is in state |1⟩.
/// It's a key building block for quantum computation and can be used to create entanglement.
///
/// The gate transforms the computational basis states as:
/// - |00⟩ → |00⟩ (control |0⟩: target unchanged)
/// - |01⟩ → |01⟩ (control |0⟩: target unchanged)  
/// - |10⟩ → |11⟩ (control |1⟩: target flipped)
/// - |11⟩ → |10⟩ (control |1⟩: target flipped)
///
/// The 4x4 unitary matrix is:
///
/// CNOT = [ [ 1, 0, 0, 0 ],
///          [ 0, 1, 0, 0 ],
///          [ 0, 0, 0, 1 ],
///          [ 0, 0, 1, 0 ] ]
///
/// # Parameters
/// * `control` - The control qubit that determines whether the operation is applied
/// * `target` - The target qubit that gets flipped when control is |1⟩
pub struct ControlledNot {
    control: Qubit,
    target: Qubit,
}

impl ControlledNot {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for ControlledNot {
    /// construct the 4x4 matrix representing the gate.
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
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("CNOT(control={}, target={})", self.control, self.target)
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

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
    control: Qubit,
    target: Qubit,
}

impl SWAP {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for SWAP {
    /// construct the 4x4 matrix representing the gate.
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

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("SWAP(control={}, target={})", self.control, self.target)
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
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
    control: Qubit,
    target: Qubit,
}

impl ISWAP {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for ISWAP {
    /// construct the 4x4 matrix representing the gate.
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

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("ISWAP(control={}, target={})", self.control, self.target)
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
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
    control: Qubit,
    target: Qubit,
}

impl FSWAP {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for FSWAP {
    /// construct the 4x4 matrix representing the gate.
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

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("FSWAP(control={}, target={})", self.control, self.target)
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
/// Represents the Controlled-Y (CY) gate, applying Pauli-Y to the target when control is |1⟩.
///
/// The CY gate applies a Pauli-Y rotation to the target qubit if and only if
/// the control qubit is in state |1⟩. The Pauli-Y gate flips the qubit state
/// and introduces a phase factor of i.
///
/// The gate transforms the computational basis states as:
/// - |00⟩ → |00⟩ (control |0⟩: target unchanged)
/// - |01⟩ → |01⟩ (control |0⟩: target unchanged)
/// - |10⟩ → i|11⟩ (control |1⟩: Y applied to target)
/// - |11⟩ → -i|10⟩ (control |1⟩: Y applied to target)
///
/// The 4x4 unitary matrix is:
///
/// CY = [ [ 1,  0,  0,   0 ],
///        [ 0,  1,  0,  -i ],
///        [ 0,  0,  0,  -i ],
///        [ 0,  0,  i,   0 ] ]
///
/// # Arguments
/// * `control` - The control qubit that determines whether Y is applied
/// * `target` - The target qubit that receives the Pauli-Y operation
pub struct ControlledPauliY {
    control: Qubit,
    target: Qubit,
}

impl ControlledPauliY {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for ControlledPauliY {
    /// construct the 4x4 matrix representing the gate.
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
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -1.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -1.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 1.0),
                Complex::new(0.0, 0.0),
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("CY(control={}, target={})", self.control, self.target)
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
/// Represents the Controlled-Z (CZ) gate, applying a phase flip when both qubits are |1⟩.
///
/// The CZ gate applies a phase factor of -1 to the |11⟩ state while leaving
/// all other computational basis states unchanged. This gate is symmetric
/// with respect to both qubits and is commonly used for creating entanglement.
///
/// The gate transforms the computational basis states as:
/// - |00⟩ → |00⟩ (unchanged)
/// - |01⟩ → |01⟩ (unchanged)
/// - |10⟩ → |10⟩ (unchanged)
/// - |11⟩ → -|11⟩ (phase flip)
///
/// The 4x4 unitary matrix is:
///
/// CZ = [ [ 1,  0, 0,  0 ],
///        [ 0,  1, 0,  0 ],
///        [ 0,  0, 1,  0 ],
///        [ 0,  0, 0, -1 ] ]
///
/// The CZ gate is self-inverse and commutes with both qubits. It's equivalent
/// to a CNOT gate with Hadamard gates applied to the target qubit.
///
/// # Parameters
/// * `control` - The first qubit in the controlled phase operation
/// * `target` - The second qubit in the controlled phase operation
pub struct ControlledPauliZ {
    control: Qubit,
    target: Qubit,
}

impl ControlledPauliZ {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for ControlledPauliZ {
    /// construct the 4x4 matrix representing the gate.
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
                Complex::new(1.0, 0.0),
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
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(-1.0, 0.0)
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("CZ(control={}, target={})", self.control, self.target)
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

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
    control: Qubit,
    target: Qubit,
}

impl EchoCrossResonance {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for EchoCrossResonance {
    /// construct the 4x4 matrix representing the gate.
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

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("ECR(control={}, target={})", self.control, self.target)
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}
