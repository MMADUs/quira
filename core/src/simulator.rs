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

use crate::bit::QuantumBit;
use crate::circuit::{QuantumCircuit, QuantumInstructions, QubitOperation};
use crate::constant::{C_ONE, C_ZERO};
use crate::endian::QubitIndexing;
use crate::io::{QuantumDebugger, QuantumState};
use crate::noise::{NoNoise, NoiseChannelOperation, NoiseModel};
use crate::operations::single_qubit::PauliX;
use crate::ops::{GateApplyExt, QuantumGate};
use crate::register::ClassicalRegister;
use crate::types::Vector;

/// The Quantum Simulator for simulating state.
pub struct QuantumSimulator<T, N = NoNoise>
where
    T: QuantumState + QuantumDebugger + Clone,
    N: NoiseChannelOperation,
{
    /// the quantum state vector.
    kernel: T,
    /// qubit indexing when applying to quantum state.
    qubit_indexing: QubitIndexing,
    /// classical register after simulation.
    classical_register: ClassicalRegister,
    /// noise model during simulation.
    noise_model: NoiseModel<N>,
}

impl<T: QuantumState + QuantumDebugger + Clone> QuantumSimulator<T, NoNoise> {
    pub fn new(state: T, indexing: QubitIndexing) -> Self {
        Self {
            kernel: state,
            qubit_indexing: indexing,
            classical_register: ClassicalRegister::new(0),
            noise_model: NoiseModel::new(),
        }
    }
}

impl<T: QuantumState + QuantumDebugger + Clone, N: NoiseChannelOperation> QuantumSimulator<T, N> {
    pub fn with_noise(&mut self, noise_model: NoiseModel<N>) {
        self.noise_model.extend(noise_model);
    }

    /// Return the number of qubit in the quantum statevec.
    pub fn num_qubits(&self) -> usize {
        self.kernel.num_qubits()
    }

    /// The kernel backend as ref.
    pub fn backend(&self) -> &T {
        &self.kernel
    }

    /// The classical register after simulation.
    pub fn creg(&self) -> &ClassicalRegister {
        &self.classical_register
    }

    /// Simulate a quantum circuit.
    pub fn simulate(&mut self, circuit: QuantumCircuit) {
        // prepare classical register
        self.classical_register = circuit.creg_as_ref().clone();
        // apply quantum operation based on tokens.
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
                        QubitOperation::RESET(qubit_idx) => {
                            let qubit = QuantumBit::new(*qubit_idx).label("Reset qubit into |0>");
                            let pauli_x = PauliX::new(&qubit);
                            if pauli_x.construct_targets().len() > self.num_qubits() {
                                panic!(
                                    "{} Qubit Operation tries to operate on {} Qubit",
                                    pauli_x.construct_targets().len(),
                                    self.num_qubits()
                                );
                            }
                            // measure qubit
                            let measured_bit = self.kernel.measure_qubit(*qubit_idx);
                            // conditionally apply
                            if measured_bit == true {
                                pauli_x.apply(&mut self.kernel, &self.qubit_indexing);
                            }
                        }
                    }
                }
                QuantumInstructions::Operations(ops) => {
                    if ops.construct_targets().len() > self.num_qubits() {
                        panic!(
                            "{} Qubit Operation tries to operate on {} Qubit",
                            ops.construct_targets().len(),
                            self.num_qubits()
                        );
                    }
                    ops.apply(&mut self.kernel, &self.qubit_indexing);
                }
                QuantumInstructions::Conditional((classical_bit, if_measured, ops)) => {
                    if ops.construct_targets().len() > self.num_qubits() {
                        panic!(
                            "{} Qubit Operation tries to operate on {} Qubit",
                            ops.construct_targets().len(),
                            self.num_qubits()
                        );
                    }
                    // get classical measurement result
                    let measured = self.classical_register.get(*classical_bit);
                    // conditionally apply
                    if *if_measured == measured {
                        ops.apply(&mut self.kernel, &self.qubit_indexing);
                    }
                }
                QuantumInstructions::Measurements((qubit, classical_bit)) => {
                    // validate qubit size
                    if *qubit >= self.num_qubits() {
                        panic!("Measurement out of index for qubit: {}", qubit);
                    }
                    // apply measurement
                    let measured_bit = self.kernel.measure_qubit(*qubit);
                    self.classical_register.set(*classical_bit, measured_bit);
                }
                QuantumInstructions::MeasureAll => {
                    let measured = self.kernel.measure_all();
                    let mut mcreg = ClassicalRegister::new(measured.len());
                    // map result automatically
                    for (i, measured_bit) in measured.iter().enumerate() {
                        mcreg.set(i, *measured_bit);
                    }
                    self.classical_register = mcreg;
                }
                QuantumInstructions::Barrier => {}
            }
        }
    }
}
