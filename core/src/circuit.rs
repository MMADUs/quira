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
    types::{Complex, Qubit, Vector},
};

/// Qubit operation types
pub enum QubitOperation {
    FROM(Vector<Complex>),
    ONES(usize),
    ZEROS(usize),
    RESET,
}

/// Quantum instructions type
pub enum QuantumInstructions {
    /// applying qubit to state
    Qubit(QubitOperation),
    /// apply gate operation
    Operations(Box<dyn QuantumGate>),
    /// conditionally apply gate operation
    /// 'classical_bit', 'measured', 'operations'.
    Conditional((usize, bool, Box<dyn QuantumGate>)),
    /// qubit measurements
    /// 'qubit_index', 'classical_bit'.
    Measurements((Qubit, usize)),
    /// measure all qubit
    MeasureAll,
    /// circuit barrier
    Barrier,
}

/// The Quantum Circuit.
pub struct QuantumCircuit {
    /// number of current qubits in circuit.
    curr_n_qubits: usize,
    /// applied quantum operations.
    operations: Vec<QuantumInstructions>,
    /// quantum registers capacity
    qreg: usize,
    /// classical registers capacity
    creg: usize,
}

impl QuantumCircuit {
    /// Create a new quantum circuit.
    pub fn new(qreg: usize, creg: usize) -> Self {
        Self {
            curr_n_qubits: 0,
            operations: Vec::new(),
            qreg,
            creg,
        }
    }

    /// Return the number of qubit in the circuit.
    pub fn num_qubits(&self) -> usize {
        self.curr_n_qubits
    }

    /// Return the number of classical register in the circuit.
    pub fn num_classical_register(&self) -> usize {
        self.creg
    }

    /// Return the number of quantum register in the circuit.
    pub fn num_quantum_register(&self) -> usize {
        self.qreg
    }

    /// Return the reference of circuit operations.
    pub(crate) fn operations_as_ref(&self) -> &Vec<QuantumInstructions> {
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
            .push(QuantumInstructions::Qubit(QubitOperation::FROM(state)));
        self.curr_n_qubits += 1;
        self.num_qubits() - 1
    }

    /// Set the qubit state to |0..0⟩.
    pub fn zeros(&mut self, num_qubits: usize) -> Vec<Qubit> {
        let start = self.curr_n_qubits;
        let end = start + num_qubits;
        self.curr_n_qubits = end;

        self.operations
            .push(QuantumInstructions::Qubit(QubitOperation::ZEROS(
                num_qubits,
            )));

        (start..end).collect()
    }

    /// Set the qubit state to |1..1⟩.
    pub fn ones(&mut self, num_qubits: usize) -> Vec<Qubit> {
        let start = self.curr_n_qubits;
        let end = start + num_qubits;
        self.curr_n_qubits = end;

        self.operations
            .push(QuantumInstructions::Qubit(QubitOperation::ONES(num_qubits)));

        (start..end).collect()
    }

    /// Empty the qubit state.
    pub fn reset(&mut self) {
        self.operations
            .push(QuantumInstructions::Qubit(QubitOperation::RESET))
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
        self.operations.push(QuantumInstructions::Conditional((
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
            .push(QuantumInstructions::Operations(Box::new(operation)))
    }

    /// Add a measurement operation at the end of circuit.
    pub fn measure(&mut self, qubit: Qubit, classical_bit: usize) {
        if qubit > self.qreg || classical_bit > self.creg {
            panic!("Registers out of range");
        }
        self.operations
            .push(QuantumInstructions::Measurements((qubit, classical_bit)))
    }

    /// Measure all qubits.
    pub fn measure_all(&mut self) {
        self.operations.push(QuantumInstructions::MeasureAll)
    }

    /// Add a barrier into circuit.
    pub fn barrier(&mut self) {
        self.operations.push(QuantumInstructions::Barrier)
    }
}
