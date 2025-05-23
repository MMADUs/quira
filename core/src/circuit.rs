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
    constant::{C_ONE, C_ZERO, EPSILON},
    endian::QubitIndexing,
    operations::QuantumGate,
    prelude::{QubitState, RunResult},
    state::QuantumState,
    types::{Complex, Qubit, Vector},
};

use rayon::prelude::*;

pub struct QuantumCircuit {
    /// state vector for the final quantum state outcome.
    base_state: QuantumState,
    /// state vector that represent the current quantum state
    /// the state will directly updated if 'direct_apply' is enabled.
    current_state: QuantumState,
    /// applied quantum circuit gates.
    gates: Vec<Box<dyn QuantumGate>>,
    /// map of Measurements.
    /// (qubit -> classical_bit)
    measurements: HashMap<Qubit, usize>,
    /// a measurement result of the current quantum state
    /// measurement will directly updated if 'direct_apply' is enabled.
    measured_state: HashMap<Qubit, bool>,
    /// toggle if operation should be directly applied to 'current_state'.
    direct_apply: bool,
    /// qubit indexing when applying to quantum state.
    qubit_indexing: QubitIndexing,
}

impl QuantumCircuit {
    /// Create a new quantum circuit.
    pub fn new(update_flag: bool, indexing: QubitIndexing) -> Self {
        Self {
            base_state: QuantumState::new(),
            current_state: QuantumState::new(),
            gates: Vec::new(),
            measurements: HashMap::new(),
            measured_state: HashMap::new(),
            direct_apply: update_flag,
            qubit_indexing: indexing,
        }
    }

    /// Return the number of qubit in the circuit.
    pub fn num_qubits(&self) -> usize {
        self.current_state.num_qubits()
    }

    /// Create a new qubit from arbitrary iterator with complex element |?⟩.
    pub fn qb_from_iter<I>(&mut self, num_iter: I) -> Qubit
    where
        I: IntoIterator<Item = Complex>,
    {
        let vector = Vector::from_iter(num_iter);
        self.qb_from_vector(vector)
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

        if self.current_state.is_empty() {
            self.current_state = QuantumState::from(state.clone());
            self.base_state = QuantumState::from(state);
        } else {
            self.current_state = self.current_state.qubit_tensor_product(state.clone());
            self.base_state = self.base_state.qubit_tensor_product(state);
        }

        self.num_qubits() - 1
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

    /// Set and replace the entire qubit state to |0..0⟩.
    pub fn zeros(&mut self, num_qubits: usize) {
        self.current_state = QuantumState::zero(num_qubits);
        self.base_state = QuantumState::zero(num_qubits)
    }

    /// Set and replace the entire qubit state to |1..1⟩.
    pub fn ones(&mut self, num_qubits: usize) {
        self.current_state = QuantumState::one(num_qubits);
        self.base_state = QuantumState::one(num_qubits)
    }

    /// Empty the qubit state.
    pub fn reset(&mut self) {
        self.current_state = QuantumState::new();
        self.base_state = QuantumState::new()
    }

    /// Conditionally add a operation to the circuit.
    pub fn cond_add<O>(&mut self, condition: bool, operation: O) -> Option<Vector<Complex>>
    where
        O: QuantumGate + 'static,
    {
        if condition { self.add(operation) } else { None }
    }

    /// Add a operation to the circuit.
    pub fn add<O>(&mut self, operation: O) -> Option<Vector<Complex>>
    where
        O: QuantumGate + 'static,
    {
        // apply unitary to state directly
        if self.direct_apply {
            operation.apply(&mut self.current_state, &self.qubit_indexing);
        }
        // save the gate for further measurement
        self.gates.push(Box::new(operation));
        // return the current state vector
        if self.direct_apply {
            Some(self.current_state.amplitudes_as_ref().clone())
        } else {
            None
        }
    }

    /// Add a measurement operation at the end of circuit.
    pub fn measure(&mut self, qubit: Qubit, classical_bit: usize) -> Option<bool> {
        if qubit > self.current_state.num_qubits() {
            panic!("Qubit index out of range");
        }
        self.measurements.insert(qubit, classical_bit);
        // check if direct apply
        if self.direct_apply {
            let bit = self.current_state.measure_qubit(qubit);
            self.measured_state.insert(qubit, bit);
            Some(bit)
        } else {
            None
        }
    }

    /// Get the current quantum state
    pub fn state_vector(&self, filter_zero: bool) -> Vec<String> {
        println!("\nVector state: ");
        for state_str in self.current_state.get_amp_state(filter_zero) {
            println!("{}", state_str);
        }
        self.current_state.get_amp_state(filter_zero)
    }

    /// Get the quantum state from specified qubit
    pub fn state_qubit(&self, qubit: Qubit) -> Option<QubitState> {
        if let Some(&measured_bit) = self.measured_state.get(&qubit) {
            println!("\nQubit {} state:", qubit);
            println!(
                "Qubit {} already collapsed to classical bit: {}",
                qubit, measured_bit as u8
            );
            return None;
        }
        let state = self.current_state.get_qubit_state(qubit);
        Some(state)
    }

    /// Execute the circuit and return the measurement counts and probabilities.
    pub fn run(&self, shots: usize) -> RunResult {
        let start_time = Instant::now();
        println!("\nExecuting quantum circuit with {} shot(s)...", shots);

        // early return if no measurements
        if self.measurements.is_empty() {
            println!("\nWarning: No measurements defined in the circuit. Skipping measurements.");
            println!("\nExecution completed in: {:.3?}", start_time.elapsed());
            return RunResult {
                classical_bit: Vec::new(),
                shots: 0,
            };
        }

        // Run shots in parallel
        let results: Vec<Vec<bool>> = (0..shots)
            .into_par_iter()
            .map(|_shot| {
                let mut state = self.base_state.clone();

                for (_i, gate) in self.gates.iter().enumerate() {
                    gate.apply(&mut state, &self.qubit_indexing);
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
        RunResult {
            classical_bit: results,
            shots,
        }
    }

    /// Count all the statistical outcome after execution with enhanced formatting and analysis
    pub fn analyze(&self, rr: RunResult) {
        if rr.classical_bit.is_empty() {
            println!("\nNo results to count.");
            return;
        }

        let num_classical_bits = self.measurements.len();

        // Collect all outcomes using parallel processing
        let counts: HashMap<String, usize> = rr
            .classical_bit
            .par_iter()
            .map(|classical_register| {
                classical_register
                    .iter()
                    .rev()
                    .map(|b| if *b { '1' } else { '0' })
                    .collect::<String>()
            })
            .fold(
                HashMap::new,
                |mut acc: HashMap<String, usize>, bit_string| {
                    *acc.entry(bit_string).or_insert(0) += 1;
                    acc
                },
            )
            .reduce(HashMap::new, |mut acc, map| {
                for (key, value) in map {
                    *acc.entry(key).or_insert(0) += value;
                }
                acc
            });

        // Sort outcomes by binary value for better readability
        let mut sorted_outcomes: Vec<(String, usize)> = counts.into_iter().collect();
        sorted_outcomes.sort_by(|a, b| {
            // Sort by binary value (convert string to number)
            let a_val = usize::from_str_radix(&a.0, 2).unwrap_or(0);
            let b_val = usize::from_str_radix(&b.0, 2).unwrap_or(0);
            a_val.cmp(&b_val)
        });

        println!("\n{:=<65}", "");
        println!("MEASUREMENT RESULTS FROM {} SHOTS", rr.shots);
        println!("{:=<65}", "");

        // Individual qubit statistics (parallel processing)
        if num_classical_bits > 1 {
            println!("\nPer-Qubit Statistics:");
            println!("{:-<40}", "");

            let qubit_stats: Vec<(usize, usize, usize)> = (0..num_classical_bits)
                .into_par_iter()
                .map(|bit_idx| {
                    let (zeros, ones) = rr
                        .classical_bit
                        .par_iter()
                        .map(|result| if result[bit_idx] { (0, 1) } else { (1, 0) })
                        .reduce(|| (0, 0), |a, b| (a.0 + b.0, a.1 + b.1));
                    (bit_idx, zeros, ones)
                })
                .collect();

            for (bit_idx, zeros, ones) in qubit_stats {
                // Find which qubit this classical bit corresponds to
                let qubit_idx = self
                    .measurements
                    .iter()
                    .find(|(_, classical_bit)| **classical_bit == bit_idx)
                    .map(|(&qubit, _)| qubit)
                    .unwrap_or(bit_idx);

                println!("Qubit {} (classical bit {}):", qubit_idx, bit_idx);
                println!(
                    "  |0⟩: {:5} ({:5.2}%)",
                    zeros,
                    (zeros as f64 / rr.shots as f64) * 100.0
                );
                println!(
                    "  |1⟩: {:5} ({:5.2}%)",
                    ones,
                    (ones as f64 / rr.shots as f64) * 100.0
                );
                println!();
            }
        }

        // Overall outcome distribution
        println!("Complete State Measurements:");
        println!("{:-<65}", "");
        println!(
            "{:>12} | {:>8} | {:>8} | {:>8}",
            "State", "Count", "Percent", "Probability"
        );
        println!("{:-<65}", "");

        let total_shots = rr.shots as f64;
        let mut total_probability = 0.0;

        for (outcome, count) in &sorted_outcomes {
            let percentage = (*count as f64 / total_shots) * 100.0;
            let probability = *count as f64 / total_shots;
            total_probability += probability;

            println!(
                "{:>12} | {:>8} | {:>7.2}% | {:>8.4}",
                outcome, count, percentage, probability
            );
        }

        println!("{:-<65}", "");
        println!(
            "Total probability: {:.6} (should be ~1.0)",
            total_probability
        );

        // Statistical analysis
        println!("\nStatistical Analysis:");
        println!("{:-<40}", "");

        let expected_count = total_shots / (1 << num_classical_bits) as f64;
        let expected_probability = 1.0 / (1 << num_classical_bits) as f64;

        println!("Expected count per state: {:.1}", expected_count);
        println!(
            "Expected probability per state: {:.4}",
            expected_probability
        );

        // Calculate chi-squared goodness of fit
        let mut chi_squared = 0.0;
        let mut max_deviation = 0.0;
        let mut max_deviation_state = String::new();

        for (outcome, count) in &sorted_outcomes {
            let observed = *count as f64;
            let deviation = (observed - expected_count).abs();
            let chi_contribution = (observed - expected_count).powi(2) / expected_count;
            chi_squared += chi_contribution;

            if deviation > max_deviation {
                max_deviation = deviation;
                max_deviation_state = outcome.clone();
            }
        }

        println!("Chi-squared statistic: {:.4}", chi_squared);
        println!(
            "Largest deviation: {:.1} counts (state {})",
            max_deviation, max_deviation_state
        );

        // Entropy calculation (parallel)
        let entropy = -sorted_outcomes
            .par_iter()
            .map(|(_, count)| {
                let p = *count as f64 / total_shots;
                if p > 0.0 { p * p.log2() } else { 0.0 }
            })
            .sum::<f64>();

        let max_entropy = (num_classical_bits as f64).log2();
        println!("Shannon entropy: {:.4} bits", entropy);
        println!("Maximum entropy: {:.4} bits", max_entropy);
        println!(
            "Entropy ratio: {:.4} (1.0 = maximum randomness)",
            entropy / max_entropy
        );

        println!("{:=<60}", "");
    }
}
