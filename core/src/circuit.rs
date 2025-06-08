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

use std::collections::HashMap;

use crate::{
    constant::{C_ONE, C_ZERO, EPSILON},
    operations::QuantumGate,
    prelude::{QuantumTokens, QubitToken},
    types::{Complex, Qubit, Vector},
};

pub struct QuantumCircuit {
    /// number of current qubits in circuit.
    current_qubits: usize,
    /// applied quantum operations.
    operations: Vec<QuantumTokens>,
    /// map of Measurements.
    /// (qubit -> classical_bit)
    registers: HashMap<Qubit, usize>,
}

impl QuantumCircuit {
    /// Create a new quantum circuit.
    pub fn new() -> Self {
        Self {
            current_qubits: 0,
            operations: Vec::new(),
            registers: HashMap::new(),
        }
    }

    /// Return the number of qubit in the circuit.
    pub fn num_qubits(&self) -> usize {
        self.current_qubits
    }

    /// Return the reference of circuit registers.
    pub(crate) fn registers_as_ref(&self) -> &HashMap<Qubit, usize> {
        &self.registers
    }

    /// Return the reference of circuit operations.
    pub(crate) fn operations_as_ref(&self) -> &Vec<QuantumTokens> {
        &self.operations
    }

    /// Create a new qubit from arbitrary iterator with complex element |?⟩.
    pub fn qb_from_iter<I>(&mut self, num_iter: I) -> Qubit
    where
        I: IntoIterator<Item = Complex>,
    {
        let vector = Vector::from_iter(num_iter);
        self.qb_from_vector(vector)
    }

    /// Create a new qubit in zero state |0⟩.
    pub fn qb_zero(&mut self) -> Qubit {
        // qubit representation in zero state as vector = [1, 0]
        let ket_zero = Vector::from_iter([C_ONE, C_ZERO]);
        self.qb_from_vector(ket_zero)
    }

    /// Create a new qubit in one state |1⟩.
    pub fn qb_one(&mut self) -> Qubit {
        // qubit representation in one state as vector = [0, 1]
        let ket_one = Vector::from_iter([C_ZERO, C_ONE]);
        self.qb_from_vector(ket_one)
    }

    /// Create a new qubit with arbitrary vector state |?⟩.
    pub fn qb_from_vector(&mut self, mut state: Vector<Complex>) -> Qubit {
        assert!(
            state.len() == 2,
            "New qubit must be a 2-dimensional vector (single qubit)"
        );

        // Normalize the vector here, since we don't know what the arbitrary input is
        let norm = state.iter().map(|x| (x * x.conj()).re).sum::<f64>().sqrt();

        if norm > EPSILON {
            state.mapv_inplace(|x| x / norm);
        } else {
            panic!("State vector norm is either zero or too small");
        }

        // Apply to circuit
        self.operations
            .push(QuantumTokens::Qubit(QubitToken::FROM(state)));
        self.current_qubits += 1;
        self.num_qubits() - 1
    }

    /// Set the qubit state to |0..0⟩.
    pub fn zeros(&mut self, num_qubits: usize) -> Vec<Qubit> {
        let start = self.current_qubits;
        let end = start + num_qubits;
        self.current_qubits = end;

        self.operations
            .push(QuantumTokens::Qubit(QubitToken::ZEROS(num_qubits)));

        (start..end).collect()
    }

    /// Set the qubit state to |1..1⟩.
    pub fn ones(&mut self, num_qubits: usize) -> Vec<Qubit> {
        let start = self.current_qubits;
        let end = start + num_qubits;
        self.current_qubits = end;

        self.operations
            .push(QuantumTokens::Qubit(QubitToken::ONES(num_qubits)));

        (start..end).collect()
    }

    /// Empty the qubit state.
    pub fn reset(&mut self) {
        self.operations
            .push(QuantumTokens::Qubit(QubitToken::RESET))
    }

    /// Conditionally add a operation to the circuit
    /// based on the bit we measured in the classical bit.
    pub fn cond_add<O>(&mut self, classical_bit: usize, if_measured: i8, operation: O)
    where
        O: QuantumGate + 'static,
    {
        if if_measured != 1 && if_measured != 0 {
            panic!("measurement condition is only between 0 and 1.")
        }
        self.operations.push(QuantumTokens::Conditional((
            classical_bit,
            if_measured != 0,
            Box::new(operation),
        )))
    }

    /// Add a operation to the circuit.
    pub fn add<O>(&mut self, operation: O)
    where
        O: QuantumGate + 'static,
    {
        self.operations
            .push(QuantumTokens::Operations(Box::new(operation)))
    }

    /// Add a measurement operation at the end of circuit.
    pub fn measure(&mut self, qubit: Qubit, classical_bit: usize) {
        if qubit > self.registers.capacity() {
            panic!("Qubit index out of range");
        }
        self.registers.insert(qubit, classical_bit);
        self.operations
            .push(QuantumTokens::Measurements((qubit, classical_bit)))
    }

    /// Add a barrier into circuit.
    pub fn barrier(&mut self) {
        self.operations.push(QuantumTokens::Barrier)
    }
}
