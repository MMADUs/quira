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

use std::time::Instant;

use crate::{
    constant::{C_ONE, C_ZERO}, endian::QubitIndexing, kernel::QuantumState, prelude::RunResult, types::Vector, GateApplyExt, QuantumCircuit, QuantumTokens, QubitToken
};

use rayon::prelude::*;

pub struct QuantumBackend<T: QuantumState + Clone> {
    /// the selected kernel (quantum state)
    kernel: T,
    /// executed quantum circuit.
    circuits: Vec<QuantumCircuit>,
    /// qubit indexing when applying to quantum state.
    qubit_indexing: QubitIndexing,
}

impl<T: QuantumState + Clone> QuantumBackend<T> {
    /// create a new quantum backend with the given state.
    pub fn new(state: T, indexing: QubitIndexing) -> Self {
        Self {
            kernel: state,
            circuits: Vec::new(),
            qubit_indexing: indexing,
        }
    }

    /// Add multiple circuits to backend.
    pub fn from_circuits<I>(&mut self, circuits: I)
    where
        I: IntoIterator<Item = QuantumCircuit>,
    {
        self.circuits.extend(circuits);
    }

    /// Add circuit to backend.
    pub fn from_circuit(&mut self, circuit: QuantumCircuit) {
        self.circuits.push(circuit)
    }

    /// Execute the circuit and return the measurement counts and probabilities.
    pub fn run(&self, shots: usize) -> Vec<RunResult> {
        let start_time = Instant::now();
        println!(
            "\nExecuting {} quantum circuit(s) with {} shot(s) each...",
            self.circuits.len(),
            shots
        );

        // Check if any circuits have measurements
        let circuits_with_measurements = self
            .circuits
            .iter()
            .filter(|c| !c.registers_as_ref().is_empty())
            .count();

        if circuits_with_measurements == 0 {
            println!("\nWarning: No measurements defined in any circuit.");
        } else {
            println!(
                "\n{} of {} circuits have measurements defined.",
                circuits_with_measurements,
                self.circuits.len()
            );
        }

        // Run circuits in parallel, then shots in parallel for each circuit
        let results: Vec<RunResult> = self
            .circuits
            .par_iter()
            .enumerate()
            .map(|(circuit_idx, circuit)| {
                println!(
                    "Running circuit {} of {}",
                    circuit_idx + 1,
                    self.circuits.len()
                );

                // If circuit has no measurements, return empty result but still run the circuit
                if circuit.registers_as_ref().is_empty() {
                    println!(
                        "  Circuit {} has no measurements - returning empty result",
                        circuit_idx + 1
                    );
                    return RunResult {
                        classical_bit: Vec::new(),
                        shots: 0,
                    };
                }

                // Run shots in parallel for this circuit
                let shot_results: Vec<Vec<bool>> = (0..shots)
                    .into_par_iter()
                    .map(|_shot| {
                        // Prepare circuit execution
                        let mut qstate = self.kernel.clone();
                        let mut classical_register: Vec<Option<bool>> =
                            vec![None; circuit.registers_as_ref().len()];
                        // Apply quantum operation based on tokens.
                        for circuit_ops in circuit.operations_as_ref() {
                            match circuit_ops {
                                QuantumTokens::Qubit(qb) => {
                                    match qb {
                                        QubitToken::FROM(vector) => {
                                            qstate.expand_state(vector.clone());
                                        }
                                        QubitToken::ONES(num_qubits) => {
                                            // qubit representation in one state as vector = [0, 1]
                                            let ket_one = Vector::from_iter([C_ZERO, C_ONE]);
                                            for _ in 0..*num_qubits {
                                                qstate.expand_state(ket_one.clone());
                                            }
                                        }
                                        QubitToken::ZEROS(num_qubits) => {
                                            // qubit representation in zero state as vector = [1, 0]
                                            let ket_zero = Vector::from_iter([C_ONE, C_ZERO]);
                                            for _ in 0..*num_qubits {
                                                qstate.expand_state(ket_zero.clone());
                                            }
                                        }
                                        QubitToken::RESET => qstate.reset_state(),
                                    }
                                }
                                QuantumTokens::Operations(ops) => {
                                    ops.apply(&mut qstate, &self.qubit_indexing);
                                }
                                QuantumTokens::Conditional((classical_bit, if_measured, ops)) => {
                                    // get classical measurement result
                                    let measured = classical_register[*classical_bit];
                                    let measured = measured.expect(&format!(
                                        "classical_bit {} is empty",
                                        classical_bit
                                    ));
                                    // conditionally apply
                                    if *if_measured == measured {
                                        ops.apply(&mut qstate, &self.qubit_indexing);
                                    }
                                }
                                QuantumTokens::Measurements((qubit, classical_bit)) => {
                                    // validate qubit size
                                    if *qubit > qstate.num_qubits() {
                                        panic!("Measurement out of index for qubit: {}", qubit);
                                    }
                                    // apply measurement
                                    let measured = qstate.measure_qubit(*qubit);
                                    classical_register[*classical_bit] = Some(measured);
                                }
                                QuantumTokens::Barrier => {}
                            }
                        }
                        // Default to false for unmeasured bits
                        classical_register
                            .into_iter()
                            .map(|opt_bit| opt_bit.unwrap_or(false))
                            .collect::<Vec<bool>>()
                    })
                    .collect();

                RunResult {
                    classical_bit: shot_results,
                    shots,
                }
            })
            .collect();

        let duration = start_time.elapsed();
        println!("\nExecution completed in: {:.3?}", duration);
        results
    }
}
