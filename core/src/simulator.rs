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
    GateApplyExt, QuantumCircuit, QuantumTokens, QubitToken,
    constant::{C_ONE, C_ZERO},
    endian::QubitIndexing,
    kernel::{QuantumDebugger, QuantumState},
    types::{Qubit, Vector},
};

pub struct QuantumSimulator<T: QuantumState + QuantumDebugger + Clone> {
    /// the quantum state vector.
    kernel: T,
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

impl<T: QuantumState + QuantumDebugger + Clone> QuantumSimulator<T> {
    /// Create a new quantum circuit.
    pub fn new(state: T, indexing: QubitIndexing) -> Self {
        Self {
            kernel: state,
            measurements: HashMap::new(),
            classical_register: HashMap::new(),
            collapsed_state: HashMap::new(),
            qubit_indexing: indexing,
        }
    }

    /// Return the number of qubit in the quantum statevec.
    pub fn num_qubits(&self) -> usize {
        self.kernel.num_qubits()
    }

    /// Get the current quantum state
    pub fn entire_state(&self, filter_zero: bool) -> T::StateType {
        self.kernel.entire_state(filter_zero)
    }

    /// Get the quantum state from specified qubit
    pub fn qubit_state(&self, qubit: Qubit) -> Option<T::QubitStateType> {
        if let Some(&measured_bit) = self.collapsed_state.get(&qubit) {
            println!("\nQubit {} state:", qubit);
            println!(
                "Qubit {} already collapsed to classical bit: {}",
                qubit, measured_bit as u8
            );
            return None;
        }
        let state = self.kernel.qubit_state(qubit);
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
                            self.kernel.expand_state(vector.clone());
                        }
                        QubitToken::ONES(num_qubits) => {
                            // qubit representation in one state as vector = [0, 1]
                            let ket_one = Vector::from_iter([C_ZERO, C_ONE]);
                            for _ in 0..*num_qubits {
                                self.kernel.expand_state(ket_one.clone());
                            }
                        }
                        QubitToken::ZEROS(num_qubits) => {
                            // qubit representation in zero state as vector = [1, 0]
                            let ket_zero = Vector::from_iter([C_ONE, C_ZERO]);
                            for _ in 0..*num_qubits {
                                self.kernel.expand_state(ket_zero.clone());
                            }
                        }
                        QubitToken::RESET => self.kernel.reset_state(),
                    }
                }
                QuantumTokens::Operations(ops) => {
                    ops.apply(&mut self.kernel, &self.qubit_indexing);
                }
                QuantumTokens::Conditional((classical_bit, if_measured, ops)) => {
                    // get classical measurement result
                    let measured = self.classical_register.get(classical_bit);
                    let measured =
                        measured.expect(&format!("classical_bit {} is empty", classical_bit));
                    // conditionally apply
                    if if_measured == measured {
                        ops.apply(&mut self.kernel, &self.qubit_indexing);
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
                    let bit = self.kernel.measure_qubit(*qubit);
                    self.collapsed_state.insert(*qubit, bit);
                    self.classical_register.insert(*classical_bit, bit);
                }
                QuantumTokens::Barrier => {}
            }
        }
    }
}
