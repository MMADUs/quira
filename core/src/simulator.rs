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
    QuantumCircuit, QuantumTokens, QubitToken,
    constant::{C_ONE, C_ZERO, EPSILON},
    endian::QubitIndexing,
    operations::QuantumGate,
    prelude::QubitState,
    statevec::QuantumStateVec,
    types::{Complex, Qubit, Vector},
};

pub struct QuantumSimulator {
    /// the quantum state vector.
    statevec: QuantumStateVec,
    /// map of measurements
    /// set the final outcome location of the qubit by classical bit.
    /// (qubit -> classical_bit)
    measurements: HashMap<Qubit, usize>,
    /// the classical register
    /// stores the final outcome based on the map of measurements.
    classical_register: HashMap<usize, bool>,
    /// store the result of the collapsed qubit state.
    collapsed_state: HashMap<Qubit, bool>,
    /// qubit indexing when applying to quantum state.
    qubit_indexing: QubitIndexing,
}

impl QuantumSimulator {
    /// Create a new quantum circuit.
    pub fn new(indexing: QubitIndexing) -> Self {
        Self {
            statevec: QuantumStateVec::new(),
            measurements: HashMap::new(),
            classical_register: HashMap::new(),
            collapsed_state: HashMap::new(),
            qubit_indexing: indexing,
        }
    }

    /// Return the number of qubit in the quantum statevec.
    pub fn num_qubits(&self) -> usize {
        self.statevec.num_qubits()
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

        if self.statevec.is_empty() {
            self.statevec = QuantumStateVec::from(state);
        } else {
            self.statevec = self.statevec.qubit_tensor_product(state);
        }
        self.num_qubits() - 1
    }

    /// Set and replace the entire qubit state to |0..0⟩.
    pub fn zeros(&mut self, num_qubits: usize) -> Vec<Qubit> {
        let start = self.num_qubits();
        let end = start + num_qubits;
        // qubit representation in zero state as vector = [1, 0]
        let ket_zero = Vector::from_iter([C_ONE, C_ZERO]);
        for _ in 0..num_qubits {
            if self.statevec.is_empty() {
                self.statevec = QuantumStateVec::from(ket_zero.clone());
            } else {
                self.statevec = self.statevec.qubit_tensor_product(ket_zero.clone());
            }
        }
        // return index
        (start..end).collect()
    }

    /// Set and replace the entire qubit state to |1..1⟩.
    pub fn ones(&mut self, num_qubits: usize) -> Vec<Qubit> {
        let start = self.num_qubits();
        let end = start + num_qubits;
        // qubit representation in one state as vector = [0, 1]
        let ket_one = Vector::from_iter([C_ZERO, C_ONE]);
        for _ in 0..num_qubits {
            if self.statevec.is_empty() {
                self.statevec = QuantumStateVec::from(ket_one.clone());
            } else {
                self.statevec = self.statevec.qubit_tensor_product(ket_one.clone());
            }
        }
        // return index
        (start..end).collect()
    }

    /// Empty the qubit state.
    pub fn reset(&mut self) {
        self.statevec = QuantumStateVec::new()
    }

    /// Conditionally add a operation to the circuit.
    pub fn cond_add<O>(&mut self, condition: bool, operation: O) -> Option<Vector<Complex>>
    where
        O: QuantumGate + 'static,
    {
        if condition {
            Some(self.add(operation))
        } else {
            None
        }
    }

    /// Add a operation to the circuit.
    pub fn add<O>(&mut self, operation: O) -> Vector<Complex>
    where
        O: QuantumGate + 'static,
    {
        operation.apply(&mut self.statevec, &self.qubit_indexing);
        self.statevec.amplitudes_as_ref().clone()
    }

    /// Add a measurement operation at the end of circuit.
    pub fn measure(&mut self, qubit: Qubit, classical_bit: usize) -> bool {
        if qubit > self.statevec.num_qubits() {
            panic!("Qubit index out of range");
        }
        if let Some(_) = self.collapsed_state.get(&qubit) {
            panic!("Unable to measure a collapsed qubit state");
        }
        if let Some(_) = self.measurements.get(&qubit) {
            println!("Warning, replacing existing qubit register");
        }
        self.measurements.insert(qubit, classical_bit);
        let bit = self.statevec.measure_qubit(qubit);
        self.collapsed_state.insert(qubit, bit);
        self.classical_register.insert(classical_bit, bit);
        bit
    }

    /// Do nothing in simulator
    pub fn barrier(&mut self) {}

    /// Get the current quantum state
    pub fn state_vec(&self, filter_zero: bool) -> Vec<String> {
        println!("\nVector state: ");
        for state_str in self.statevec.get_amp_state(filter_zero) {
            println!("{}", state_str);
        }
        self.statevec.get_amp_state(filter_zero)
    }

    /// Get the quantum state from specified qubit
    pub fn state_qubit(&self, qubit: Qubit) -> Option<QubitState> {
        if let Some(&measured_bit) = self.collapsed_state.get(&qubit) {
            println!("\nQubit {} state:", qubit);
            println!(
                "Qubit {} already collapsed to classical bit: {}",
                qubit, measured_bit as u8
            );
            return None;
        }
        let state = self.statevec.get_qubit_state(qubit);
        Some(state)
    }

    /// Simulate a quantum circuit.
    pub fn from_circuit(&mut self, circuit: QuantumCircuit) {
        // Apply quantum operation based on tokens.
        for circuit_ops in circuit.operations_as_ref() {
            match circuit_ops {
                QuantumTokens::Qubit(qb) => {
                    match qb {
                        QubitToken::FROM(vector) => {
                            if self.statevec.is_empty() {
                                self.statevec = QuantumStateVec::from(vector.clone());
                            } else {
                                self.statevec = self.statevec.qubit_tensor_product(vector.clone())
                            }
                        }
                        QubitToken::ONES(num_qubits) => {
                            // qubit representation in one state as vector = [0, 1]
                            let ket_one = Vector::from_iter([C_ZERO, C_ONE]);
                            for _ in 0..*num_qubits {
                                if self.statevec.is_empty() {
                                    self.statevec = QuantumStateVec::from(ket_one.clone());
                                } else {
                                    self.statevec =
                                        self.statevec.qubit_tensor_product(ket_one.clone());
                                }
                            }
                        }
                        QubitToken::ZEROS(num_qubits) => {
                            // qubit representation in zero state as vector = [1, 0]
                            let ket_zero = Vector::from_iter([C_ONE, C_ZERO]);
                            for _ in 0..*num_qubits {
                                if self.statevec.is_empty() {
                                    self.statevec = QuantumStateVec::from(ket_zero.clone());
                                } else {
                                    self.statevec =
                                        self.statevec.qubit_tensor_product(ket_zero.clone());
                                }
                            }
                        }
                        QubitToken::RESET => self.statevec = QuantumStateVec::new(),
                    }
                }
                QuantumTokens::Operations(ops) => {
                    ops.apply(&mut self.statevec, &self.qubit_indexing);
                }
                QuantumTokens::Conditional((classical_bit, if_measured, ops)) => {
                    // get classical measurement result
                    let measured = self.classical_register.get(classical_bit);
                    let measured =
                        measured.expect(&format!("classical_bit {} is empty", classical_bit));
                    // conditionally apply
                    if if_measured == measured {
                        ops.apply(&mut self.statevec, &self.qubit_indexing);
                    }
                }
                QuantumTokens::Measurements((qubit, classical_bit)) => {
                    // validate qubit size
                    if *qubit > self.num_qubits() {
                        panic!("Measurement out of index for qubit: {}", qubit);
                    }
                    if let Some(_) = self.collapsed_state.get(qubit) {
                        panic!("Unable to measure a collapsed qubit state");
                    }
                    if let Some(_) = self.measurements.get(qubit) {
                        println!("Warning, replacing existing qubit register");
                    }
                    self.measurements.insert(*qubit, *classical_bit);
                    // apply measurement
                    let bit = self.statevec.measure_qubit(*qubit);
                    self.collapsed_state.insert(*qubit, bit);
                    self.classical_register.insert(*classical_bit, bit);
                }
                QuantumTokens::Barrier => {}
            }
        }
    }
}
