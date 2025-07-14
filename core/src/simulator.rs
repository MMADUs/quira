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
    GateApplyExt, QuantumCircuit, QuantumInstructions, QubitOperation,
    constant::{C_ONE, C_ZERO},
    endian::QubitIndexing,
    kernel::{QuantumDebugger, QuantumState},
    types::{Qubit, Vector},
};

pub struct QuantumSimulator<T: QuantumState + QuantumDebugger + Clone> {
    /// the quantum state vector.
    kernel: T,
    /// qubit indexing when applying to quantum state.
    qubit_indexing: QubitIndexing,
}

impl<T: QuantumState + QuantumDebugger + Clone> QuantumSimulator<T> {
    /// Create a new quantum circuit.
    pub fn new(state: T, indexing: QubitIndexing) -> Self {
        Self {
            kernel: state,
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

    /// Get the state probability for specific BIT
    pub fn bit_state(&self, bits: usize) -> T::BitStateType {
        let bits = bit_to_bool_slice(bits, self.num_qubits());
        self.kernel.bit_state(&bits)
    }

    /// Get the quantum state from specified qubit
    pub fn qubit_state(&self, qubit: Qubit) -> T::QubitStateType {
        self.kernel.qubit_state(qubit)
    }

    /// Simulate a quantum circuit.
    pub fn from_circuit(&mut self, circuit: QuantumCircuit) {
        let mut classical_register: Vec<Option<bool>> =
            vec![None; circuit.num_classical_register()];
        // Apply quantum operation based on tokens.
        for circuit_ops in circuit.operations_as_ref() {
            match circuit_ops {
                QuantumInstructions::Qubit(qb) => {
                    match qb {
                        QubitOperation::FROM(vector) => {
                            self.kernel.expand_state(vector.clone());
                        }
                        QubitOperation::ONES(num_qubits) => {
                            // qubit representation in one state as vector = [0, 1]
                            let ket_one = Vector::from_iter([C_ZERO, C_ONE]);
                            for _ in 0..*num_qubits {
                                self.kernel.expand_state(ket_one.clone());
                            }
                        }
                        QubitOperation::ZEROS(num_qubits) => {
                            // qubit representation in zero state as vector = [1, 0]
                            let ket_zero = Vector::from_iter([C_ONE, C_ZERO]);
                            for _ in 0..*num_qubits {
                                self.kernel.expand_state(ket_zero.clone());
                            }
                        }
                        QubitOperation::RESET => self.kernel.reset_state(),
                    }
                }
                QuantumInstructions::Operations(ops) => {
                    ops.apply(&mut self.kernel, &self.qubit_indexing);
                }
                QuantumInstructions::Conditional((classical_bit, if_measured, ops)) => {
                    // get classical measurement result
                    let measured = classical_register[*classical_bit];
                    let measured =
                        measured.expect(&format!("classical_bit {} is empty", classical_bit));
                    // conditionally apply
                    if *if_measured == measured {
                        ops.apply(&mut self.kernel, &self.qubit_indexing);
                    }
                }
                QuantumInstructions::Measurements((qubit, classical_bit)) => {
                    // validate qubit size
                    if *qubit > self.num_qubits() {
                        panic!("Measurement out of index for qubit: {}", qubit);
                    }
                    // apply measurement
                    let measured_bit = self.kernel.measure_qubit(*qubit);
                    classical_register[*classical_bit] = Some(measured_bit);
                }
                QuantumInstructions::MeasureAll => {
                    let measured = self.kernel.measure_all();
                    // map result automatically
                    for (i, measured_bit) in measured.iter().enumerate() {
                        classical_register[i] = Some(*measured_bit);
                    }
                }
                QuantumInstructions::Barrier => {}
            }
        }
    }
}

fn bit_to_bool_slice(bit: usize, width: usize) -> Vec<bool> {
    (0..width).rev().map(|i| (bit >> i) & 1 == 1).collect()
}
