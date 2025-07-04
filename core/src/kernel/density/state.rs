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
//! along with this program.  If not, see <http://www.gnu.org/licenses/>.\

use ndarray_linalg::Scalar;
use rand::Rng;

use crate::{
    Complex, Matrix, Qubit, Vector,
    constant::EPSILON,
    kernel::{BackendOperation, QuantumDebugger, QuantumState},
    ops,
};

use super::matrix::Density;

#[derive(Clone)]
pub struct DensityState {
    density: Density,
    num_qubits: usize,
}

impl DensityState {
    /// initialize empty state
    pub fn new() -> Self {
        Self {
            density: Density::new(Matrix::zeros((0, 0))),
            num_qubits: 0,
        }
    }

    /// replace existing density matrix
    pub fn set(&mut self, matrix: Matrix<Complex>) {
        let (rows, cols) = matrix.dim();
        assert_eq!(rows, cols, "Hamiltonian must be square");

        let num_qubits = (rows as f64).log2() as usize;
        assert_eq!(
            1 << num_qubits,
            rows,
            "Hamiltonian dimension must be 2^n for n qubits"
        );

        self.density = Density::new(matrix);
    }

    /// initialize from a pure state vector
    pub fn from_statevec(state: &Vector<Complex>) -> Self {
        let num_qubits = (state.len() as f64).log2() as usize;
        assert_eq!(
            1 << num_qubits,
            state.len(),
            "State vector length must be 2^n"
        );

        Self {
            density: Density::from_pure_statevec(state),
            num_qubits,
        }
    }

    /// Initialize from mixed states with probabilities
    pub fn from_mixed_states(states: &[Vector<Complex>], probs: &[f64]) -> Self {
        assert!(!states.is_empty(), "States cannot be empty");
        let num_qubits = (states[0].len() as f64).log2() as usize;

        Self {
            density: Density::from_mixed_statevec(states, probs),
            num_qubits,
        }
    }

    /// Initialize in maximally mixed state
    pub fn maximally_mixed(num_qubits: usize) -> Self {
        let dim = 1 << num_qubits;
        Self {
            density: Density::maximally_mixed(dim),
            num_qubits,
        }
    }

    /// Initialize from a thermal state given a Hamiltonian and inverse temperature
    /// The Hamiltonian should act on the full multi-qubit system
    pub fn from_thermal_state(hamiltonian: &Matrix<Complex>, beta: f64) -> Self {
        let (rows, cols) = hamiltonian.dim();
        assert_eq!(rows, cols, "Hamiltonian must be square");

        // Calculate number of qubits from Hamiltonian dimension
        let num_qubits = (rows as f64).log2() as usize;
        assert_eq!(
            1 << num_qubits,
            rows,
            "Hamiltonian dimension must be 2^n for n qubits"
        );

        Self {
            density: Density::thermal_state(hamiltonian, beta),
            num_qubits,
        }
    }

    /// check if density state is empty
    pub fn is_empty(&self) -> bool {
        self.density.matrix_as_ref().shape() == &[0, 0]
    }

    /// Get number of qubits
    pub fn num_qubits(&self) -> usize {
        self.num_qubits
    }

    /// get the reference of the density
    pub fn density_as_ref(&self) -> &Density {
        &self.density
    }

    /// Get all basis state probabilities
    fn density_probability(&self) -> Vec<f64> {
        let matrix = self.density.matrix_as_ref();
        let dim = 1 << self.num_qubits;
        (0..dim).map(|i| matrix[[i, i]].re()).collect()
    }

    /// Calculate measurement probability for a qubit
    fn calculate_measurement_probability(&self, qubit: usize, outcome: bool) -> f64 {
        let matrix = self.density.matrix_as_ref();
        let dim = 1 << self.num_qubits;
        let mut probability = 0.0;

        for i in 0..dim {
            let qubit_bit = (i >> qubit) & 1;
            let matches_outcome = (qubit_bit == 1) == outcome;
            if matches_outcome {
                probability += matrix[[i, i]].re();
            }
        }

        probability.max(0.0) // Ensure non-negative
    }

    /// Apply measurement projector and renormalize
    fn apply_measurement_projector(&mut self, qubit: usize, outcome: bool) {
        let dim = 1 << self.num_qubits;
        let matrix = self.density.matrix_as_ref().clone();
        let mut new_matrix = Matrix::zeros((dim, dim));

        // Calculate measurement probability for normalization
        let prob = self.calculate_measurement_probability(qubit, outcome);
        assert!(
            prob > EPSILON,
            "Measurement probability too small for outcome"
        );

        // Apply projector: P ρ P where P projects onto the measured outcome
        for i in 0..dim {
            for j in 0..dim {
                let i_bit = (i >> qubit) & 1;
                let j_bit = (j >> qubit) & 1;
                let i_matches = (i_bit == 1) == outcome;
                let j_matches = (j_bit == 1) == outcome;

                if i_matches && j_matches {
                    new_matrix[[i, j]] = matrix[[i, j]] / prob;
                }
            }
        }

        self.density = Density::new(new_matrix);
    }

    /// Collapse the state to a specific classical basis state
    /// Used internally after measuring all qubits
    fn collapse_to_classical_state(&mut self, basis_state: usize) {
        let dim = 1 << self.num_qubits;
        assert!(
            basis_state < dim,
            "Basis state {} out of range (0-{})",
            basis_state,
            dim - 1
        );

        // Create a new density matrix with only the measured state
        let mut new_matrix = Matrix::zeros((dim, dim));
        new_matrix[[basis_state, basis_state]] = Complex::new(1.0, 0.0);

        self.density = Density::new(new_matrix);
    }
}

impl QuantumState for DensityState {
    /// return number of qubits
    fn num_qubits(&self) -> usize {
        self.num_qubits()
    }

    /// Apply this gate to the given quantum state
    fn apply(&mut self, u: Matrix<Complex>) {
        let current_state = self.density_as_ref().matrix_as_ref();
        let u_dagger = ops::dagger(&u);
        let new_state = u.dot(current_state).dot(&u_dagger);
        self.set(new_state);
    }

    /// Dynamic boxed clone
    fn box_clone(&self) -> Box<dyn QuantumState> {
        Box::new(self.clone())
    }
}

impl BackendOperation for DensityState {
    /// reset the state vector
    fn reset_state(&mut self) {
        let matrix = Matrix::zeros((0, 0));
        self.density.replace_matrix(matrix);
    }

    /// Expand the system by adding a new qubit in a given state
    fn expand_state(&mut self, qubit_state: Vector<Complex>) {
        assert_eq!(qubit_state.len(), 2, "Qubit state must be 2-dimensional");

        let density = Density::from_pure_statevec(&qubit_state);

        if self.is_empty() {
            self.density = density;
            self.num_qubits = 1;
        } else {
            self.density = self.density.tensor_product(&density);
            self.num_qubits += 1;
        }
    }

    /// Measure a single qubit and collapse the state
    fn measure_qubit(&mut self, qubit: Qubit) -> bool {
        assert!(
            qubit < self.num_qubits,
            "Qubit {} out of range (0-{})",
            qubit,
            self.num_qubits - 1
        );

        // Calculate probabilities for |0⟩ and |1⟩
        let prob_0 = self.calculate_measurement_probability(qubit, false);
        let prob_1 = self.calculate_measurement_probability(qubit, true);

        // Ensure probabilities sum to 1 (within numerical precision)
        let total_prob = prob_0 + prob_1;
        assert!(
            (total_prob - 1.0).abs() < EPSILON,
            "Measurement probabilities don't sum to 1: {} + {} = {}",
            prob_0,
            prob_1,
            total_prob
        );

        // Generate random measurement outcome
        let mut rng = rand::rng();
        let rand_val: f64 = rng.random();
        let outcome = rand_val < prob_1;

        // Apply measurement projector and renormalize
        self.apply_measurement_projector(qubit, outcome);

        outcome
    }

    /// Measure all qubits in the computational basis and return the results
    /// Returns a vector of measurement outcomes where index i corresponds to qubit i
    /// After measurement, the state collapses to the measured classical state
    fn measure_all(&mut self) -> Vec<bool> {
        if self.num_qubits == 0 {
            return Vec::new();
        }

        // Get probabilities for all basis states
        let probabilities = self.density_probability();

        // Ensure probabilities sum to 1
        let total_prob: f64 = probabilities.iter().sum();
        assert!(
            (total_prob - 1.0).abs() < EPSILON,
            "Total probability doesn't sum to 1: {}",
            total_prob
        );

        // Generate random measurement outcome
        let mut rng = rand::rng();
        let rand_val: f64 = rng.random();

        // Find which basis state was measured using cumulative probability
        let mut cumulative_prob = 0.0;
        let mut measured_state = 0;

        for (i, &prob) in probabilities.iter().enumerate() {
            cumulative_prob += prob;
            if rand_val <= cumulative_prob {
                measured_state = i;
                break;
            }
        }

        // Convert the measured state index to individual qubit outcomes
        let mut outcomes = Vec::with_capacity(self.num_qubits);
        for qubit in 0..self.num_qubits {
            let bit = (measured_state >> qubit) & 1;
            outcomes.push(bit == 1);
        }

        // Collapse the state to the measured outcome
        self.collapse_to_classical_state(measured_state);

        outcomes
    }
}

impl QuantumDebugger for DensityState {
    type StateType = Matrix<Complex>;
    type BitStateType = Complex; // Diagonal element (probability as real part)
    type QubitStateType = (f64, f64); // (prob_0, prob_1) pair - real probabilities

    /// Returns the entire density matrix
    fn entire_state(&self, _filter_zero: bool) -> Self::StateType {
        self.density.matrix_as_ref().clone()
    }

    /// Returns the diagonal element (probability) for a specific computational basis state
    /// Following Qiskit: diagonal elements are the probabilities
    fn bit_state(&self, bits: &[bool]) -> Self::BitStateType {
        let num_qubits = self.num_qubits();
        if bits.len() != num_qubits {
            panic!(
                "Error: Expected {} bits for {}-qubit system, got {}",
                num_qubits,
                num_qubits,
                bits.len()
            );
        }

        // Convert boolean array to state index
        let mut state_index = 0;
        for (i, &bit) in bits.iter().enumerate() {
            if bit {
                state_index |= 1 << i;
            }
        }

        // Return the diagonal element (which contains the probability in the real part)
        let matrix = self.density.matrix_as_ref();
        matrix[[state_index, state_index]]
    }

    /// Returns (prob_0, prob_1) marginal probabilities for a specific qubit
    /// Following Qiskit: return actual probabilities, not complex amplitudes
    fn qubit_state(&self, qubit: Qubit) -> Self::QubitStateType {
        let num_qubits = self.num_qubits();
        if qubit >= num_qubits {
            panic!(
                "Error: Qubit index {} out of range. System has {} qubits.",
                qubit, num_qubits
            );
        }

        let prob_0 = self.calculate_measurement_probability(qubit, false);
        let prob_1 = self.calculate_measurement_probability(qubit, true);

        (prob_0, prob_1)
    }
}
