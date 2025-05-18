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
use std::time::Instant;

use crate::{
    constant::{C_ONE, C_ZERO},
    operations::QuantumGate,
    state::QuantumState,
    types::{Complex, Qubit, Vector},
};

use rayon::prelude::*;

pub struct QuantumCircuit {
    /// state vector for the final quantum state outcome
    base_state: QuantumState,
    /// state vector that represent the current quantum state
    current_state: QuantumState,
    /// applied quantum circuit gates
    gates: Vec<Box<dyn QuantumGate>>,
    /// map of measurements (qubit_index -> classical_bit_index)
    measurements: HashMap<usize, usize>,
    /// toggle if operation directly applied to state
    direct_apply: bool,
}

impl QuantumCircuit {
    /// Create a new quantum circuit.
    pub fn new(update_flag: bool) -> Self {
        Self {
            base_state: QuantumState::new(),
            current_state: QuantumState::new(),
            gates: Vec::new(),
            measurements: HashMap::new(),
            direct_apply: update_flag,
        }
    }

    /// Return the number of qubit in the circuit.
    pub fn num_qubits(&self) -> usize {
        self.current_state.num_qubits()
    }

    /// Create a new qubit with arbitrary state |?⟩.
    pub fn qb_new(&mut self, state: Vector<Complex>) -> Qubit {
        assert!(
            state.len() == 2,
            "New qubit must be a 2-dimensional vector (single qubit)"
        );
        if self.current_state.is_empty() {
            self.current_state = QuantumState::from(state.clone());
            self.base_state = QuantumState::from(state);
        } else {
            self.current_state = self.current_state.tensor_product(state.clone());
            self.base_state = self.base_state.tensor_product(state);
        }
        self.num_qubits() - 1
    }

    /// Create a new qubit in zero state |0⟩.
    pub fn qb_zero(&mut self) -> Qubit {
        let ket_zero = Vector::from_vec(vec![C_ONE, C_ZERO]);
        self.qb_new(ket_zero)
    }

    /// Create a new qubit in one state |1⟩.
    pub fn qb_one(&mut self) -> Qubit {
        let ket_one = Vector::from_vec(vec![C_ZERO, C_ONE]);
        self.qb_new(ket_one)
    }

    /// Set and replace the entire qubit state to |0..0⟩.
    pub fn zero(&mut self, num_qubits: usize) {
        self.current_state = QuantumState::zero(num_qubits);
        self.base_state = QuantumState::zero(num_qubits)
    }

    /// Set and replace the entire qubit state to |1..1⟩.
    pub fn one(&mut self, num_qubits: usize) {
        self.current_state = QuantumState::one(num_qubits);
        self.base_state = QuantumState::one(num_qubits)
    }

    /// Empty the qubit state.
    pub fn reset(&mut self) {
        self.current_state = QuantumState::new();
        self.base_state = QuantumState::new()
    }

    /// Add a operation to the circuit.
    pub fn add_operation<O>(&mut self, operation: O) -> &mut Self
    where
        O: QuantumGate + 'static,
    {
        if self.direct_apply {
            operation.apply(&mut self.current_state);
        }
        self.gates.push(Box::new(operation));
        self
    }

    /// Add a measurement operation at the end of circuit.
    pub fn add_measurement(&mut self, qubit: usize, classical_bit: usize) -> &mut Self {
        if qubit > self.current_state.num_qubits() {
            panic!("Qubit index out of range");
        }
        self.measurements.insert(qubit, classical_bit);
        self
    }

    /// Execute the circuit and return the measurement counts and probabilities.
    pub fn execute(&self, shots: usize) -> (HashMap<String, usize>, HashMap<String, f64>) {
        let start_time = Instant::now();

        println!("Executing quantum circuit with {} shot(s)...", shots);
        println!("\nInitial state: ");
        for state_str in self.state() {
            println!("{}", state_str);
        }
        println!("\nTotal Gates to apply: {}", self.gates.len());
        println!("Measurements: {:?}", self.measurements);

        // Run shots in parallel
        let results: Vec<Vec<bool>> = (0..shots)
            .into_par_iter()
            .map(|_shot| {
                let mut state = self.base_state.clone();

                for (_i, gate) in self.gates.iter().enumerate() {
                    gate.apply(&mut state);
                }

                let num_qubits = self.num_qubits();
                let mut classical_register = vec![false; self.measurements.len()];

                if self.measurements.len() == num_qubits {
                    let measured_state = state.measure_all();
                    for (qubit, classical_bit) in &self.measurements {
                        classical_register[*classical_bit] = (measured_state & (1 << qubit)) != 0;
                    }
                } else {
                    for (qubit, classical_bit) in &self.measurements {
                        classical_register[*classical_bit] = state.measure_qubit(*qubit);
                    }
                }

                classical_register
            })
            .collect();

        let duration = start_time.elapsed();
        println!("\nExecution completed in: {:.3?}", duration);

        // Count outcomes
        let mut counts: HashMap<String, usize> = HashMap::new();

        for classical_register in &results {
            let bit_string = classical_register
                .iter()
                .rev()
                .map(|b| if *b { '1' } else { '0' })
                .collect::<String>();
            *counts.entry(bit_string).or_insert(0) += 1;
        }

        // Count |0⟩ and |1⟩ for first bit
        let mut zeros = 0;
        let mut ones = 0;
        for result in &results {
            if result[0] {
                ones += 1;
            } else {
                zeros += 1;
            }
        }

        println!("\nMeasurement results from {} iterations:", shots);
        println!(
            "Measured |0⟩: {} times ({:.2}%)",
            zeros,
            (zeros as f64 / shots as f64) * 100.0
        );
        println!(
            "Measured |1⟩: {} times ({:.2}%)",
            ones,
            (ones as f64 / shots as f64) * 100.0
        );

        // Compute probabilities
        let total = shots as f64;
        let probabilities = counts
            .iter()
            .map(|(k, v)| (k.clone(), *v as f64 / total))
            .collect::<HashMap<_, _>>();

        println!("\nFinal counts:");
        for (outcome, count) in &counts {
            println!("  {} => {}", outcome, count);
        }

        println!("\nFinal probabilities:");
        for (outcome, prob) in &probabilities {
            println!("  {} => {:.4}", outcome, prob);
        }

        (counts, probabilities)
    }

    /// Get the current quantum state
    pub fn state(&self) -> Vec<String> {
        self.current_state.get_state()
    }
}
