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

use ndarray_linalg::{Eigh, Scalar, Trace, UPLO};
use rand::Rng;

use crate::bit::QuantumBit;
use crate::constant::EPSILON;
use crate::io::{BackendOperation, QuantumDebugger, QuantumState};
use crate::operations::GateType;
use crate::types::{Complex, Matrix, Vector};
use crate::utils::{eigen, ops};

#[derive(Clone)]
pub struct Density {
    num_qubits: usize,
    matrix: Matrix<Complex>,
    dim: usize,
}

impl Density {
    /// new density matrix for mixed state
    pub fn new() -> Self {
        // set to |0> on 1 qubit by default
        let matrix: Matrix<Complex> = Matrix::<Complex>::zeros((0, 0));
        let (rows, _) = matrix.dim();
        Self {
            num_qubits: 0,
            matrix,
            dim: rows,
        }
    }

    /// new density matrix from a matrix
    pub fn from(matrix: Matrix<Complex>) -> Self {
        let (rows, cols) = matrix.dim();
        assert_eq!(rows, cols, "matrix must be square");
        let num_qubits = calculate_num_qubits(rows);
        // new density
        let density = Self {
            num_qubits,
            matrix,
            dim: rows,
        };
        assert!(density.is_valid(), "matrix is not a valid density matrix");
        density
    }

    /// new density matrix from pure state vector
    pub fn from_pure_statevec(state: &Vector<Complex>) -> Self {
        let dim = state.len();
        let num_qubits = calculate_num_qubits(dim);
        let mut matrix = Matrix::zeros((dim, dim));
        for i in 0..dim {
            for j in 0..dim {
                matrix[[i, j]] = state[i] * state[j].conj();
            }
        }
        Self {
            num_qubits,
            matrix,
            dim,
        }
    }

    /// new density matrix from state vector with mixed probabilities
    pub fn from_mixed_statevec(states: &[Vector<Complex>], probs: &[f64]) -> Self {
        assert_eq!(
            states.len(),
            probs.len(),
            "states and probabilities must have same length"
        );
        // normalize
        let prob_sum: f64 = probs.iter().sum();
        assert!(
            (prob_sum - 1.0).abs() < EPSILON,
            "probabilities must sum up to 1"
        );
        assert!(!states.is_empty(), "empty value length");
        let dim = states[0].len();
        let num_qubits = calculate_num_qubits(dim);
        let mut matrix = Matrix::zeros((dim, dim));
        for (state, &prob) in states.iter().zip(probs.iter()) {
            assert!(state.len() == dim, "all state must have equal dimension");
            let pure_density = Self::from_pure_statevec(state);
            matrix = matrix + pure_density.matrix * Complex::new(prob, 0.0);
        }
        Self {
            num_qubits,
            matrix,
            dim,
        }
    }

    /// new density matrix in a maximally mixed state
    pub fn maximally_mixed(dim: usize) -> Self {
        let num_qubits = calculate_num_qubits(dim);
        let factor = Complex::new(1.0 / dim as f64, 0.0);
        let matrix = Matrix::eye(dim) * factor;
        Self {
            num_qubits,
            matrix,
            dim,
        }
    }

    /// new density matrix in a thermal state
    pub fn thermal_state(hamiltonian: &Matrix<Complex>, beta: f64) -> Self {
        let (rows, cols) = hamiltonian.dim();
        assert_eq!(rows, cols, "hamiltonian must be square");
        let num_qubits = calculate_num_qubits(rows);
        // Compute hamiltonian eigen decomposition
        let (eigenvals, eigenvecs) = hamiltonian
            .eigh(UPLO::Lower)
            .expect("failed to compute eigenvalues");
        // Compute exp(-β * eigenvalues)
        let exp_eigenvals: Vector<f64> = eigenvals.mapv(|e| (-1.0 * beta * e).exp());
        // Compute partition function Z = Tr(exp(-βH))
        let partition_function: f64 = exp_eigenvals.sum();
        // Construct thermal state density matrix
        // ρ = (1/Z) * Σᵢ exp(-βEᵢ) |ψᵢ⟩⟨ψᵢ|
        // In matrix form: ρ = U * diag(exp(-βEᵢ)/Z) * U†
        let probs = &exp_eigenvals / partition_function;
        let diag_matrix = Matrix::from_diag(&probs.mapv(|p| Complex::new(p, 0.0)));
        let rho = eigenvecs
            .dot(&diag_matrix)
            .dot(&eigenvecs.t().mapv(|x| x.conj()));
        Self {
            num_qubits,
            matrix: rho,
            dim: rows,
        }
    }

    /// replace existing density matrix
    pub fn set(&mut self, matrix: Matrix<Complex>) {
        let (rows, cols) = matrix.dim();
        assert_eq!(rows, cols, "Matrix must be square");
        let num_qubits = calculate_num_qubits(rows);
        self.num_qubits = num_qubits;
        self.dim = rows;
        self.matrix = matrix;
    }

    /// check if density state is empty
    pub fn is_empty(&self) -> bool {
        self.matrix_as_ref().shape() == &[0, 0]
    }

    /// Get number of qubits
    pub fn num_qubits(&self) -> usize {
        self.num_qubits
    }

    /// Unitary evolution U ρ U†
    pub fn evolve_unitary(&self, unitary: &Matrix<Complex>) -> Self {
        let u_rho = unitary.dot(self.matrix_as_ref());
        let u_dagger = ops::dagger(unitary);
        let evolved = u_rho.dot(&u_dagger);
        Self {
            num_qubits: self.num_qubits,
            matrix: evolved,
            dim: self.dim(),
        }
    }

    /// Apply Kraus operators: ρ' = Σᵢ Kᵢ ρ Kᵢ†
    pub fn apply_kraus_ops(&self, kraus_ops: &[Matrix<Complex>]) -> Self {
        let mut result = Matrix::zeros((self.dim(), self.dim()));
        // apply operations
        for kraus_op in kraus_ops {
            let k_rho = kraus_op.dot(self.matrix_as_ref());
            let k_dagger = ops::dagger(&kraus_op);
            let k_rho_k_dag = k_rho.dot(&k_dagger);
            result = result + k_rho_k_dag;
        }
        Self {
            num_qubits: self.num_qubits,
            matrix: result,
            dim: self.dim(),
        }
    }

    /// Get reference to the density matrix
    pub fn matrix_as_ref(&self) -> &Matrix<Complex> {
        &self.matrix
    }

    /// Get dimension of the density matrix
    pub fn dim(&self) -> usize {
        self.dim
    }

    /// compute density trace
    pub fn trace(&self) -> Complex {
        self.matrix.trace().unwrap_or(Complex::new(0.0, 0.0))
    }

    /// check if matrix is a valid density matrix
    pub fn is_valid(&self) -> bool {
        self.is_hermitian()
            && self.is_positive_semidefinite()
            && (self.trace().re - 1.0).abs() < EPSILON
            && self.trace().im.abs() < EPSILON
    }

    /// Check if matrix is Hermitian
    pub(crate) fn is_hermitian(&self) -> bool {
        for i in 0..self.dim {
            for j in 0..self.dim {
                if (self.matrix[[i, j]] - self.matrix[[j, i]].conj()).abs() > EPSILON {
                    return false;
                }
            }
        }
        true
    }

    /// Check if matrix is positive semidefinite using eigenvalue decomposition
    pub(crate) fn is_positive_semidefinite(&self) -> bool {
        let eigenvals = eigen::eigen_values(&self.matrix);
        eigenvals.iter().all(|&ev| ev >= -1.0 * EPSILON)
    }

    /// Compute purity Tr(ρ²)
    pub fn purity(&self) -> f64 {
        let rho_squared = self.matrix.dot(&self.matrix);
        rho_squared.diag().sum().re()
    }

    /// Check if state is pure
    pub fn is_pure(&self) -> bool {
        let purity = self.purity();
        (purity - 1.0).abs() < EPSILON
    }

    /// Von Neumann entropy S(ρ) = -Tr(ρ log ρ)
    pub(crate) fn von_neumann_entropy(&self) -> f64 {
        let eigenvals = eigen::eigen_values(&self.matrix);
        let mut entropy = 0.0;
        for &eigenval in &eigenvals {
            if eigenval > EPSILON {
                entropy -= eigenval * eigenval.ln();
            }
        }
        entropy
    }

    /// Multiply density matrix by scalar
    pub fn scale(&self, scalar: f64) -> Self {
        let matrix = &self.matrix * Complex::new(scalar, 0.0);
        Self {
            num_qubits: self.num_qubits,
            matrix,
            dim: self.dim,
        }
    }

    /// Normalize density matrix to have trace 1
    pub fn normalize(&self) -> Self {
        let trace = self.trace().re().abs();
        assert!(trace > EPSILON, "cannot normalize zero trace matrix");
        self.scale(1.0 / trace)
    }

    /// Tensor product with other density
    pub fn tensor_product(&self, other: &Density) -> Self {
        let matrix = ops::kron(&self.matrix, &other.matrix);
        let dim = self.dim * other.dim;
        let num_qubits = self.num_qubits + other.num_qubits;
        Self {
            num_qubits,
            matrix,
            dim,
        }
    }

    /// Helper function to convert flat index to multi-index
    pub(crate) fn multi_index(&self, flat_idx: usize, dims: &[usize]) -> Vec<usize> {
        let mut result = Vec::with_capacity(dims.len());
        let mut remaining = flat_idx;
        for &dim in dims.iter().rev() {
            result.push(remaining % dim);
            remaining /= dim;
        }
        result.reverse();
        result
    }

    // Helper function to convert multi-index to flat index
    pub(crate) fn flat_index(&self, multi_idx: &[usize], dims: &[usize]) -> usize {
        let mut result = 0;
        let mut multiplier = 1;
        for (i, &idx) in multi_idx.iter().enumerate().rev() {
            result += idx * multiplier;
            multiplier *= dims[i];
        }
        result
    }

    /// Partial trace over specified subsystem indices
    /// subsystem_dims: dimensions of each subsystem
    /// trace_over: indices of subsystems to trace over
    pub fn partial_trace(&self, subsystem_dims: &[usize], trace_over: &[usize]) -> Density {
        let total_dim: usize = subsystem_dims.iter().product();
        assert_eq!(
            self.dim, total_dim,
            "subsystem dimensions don't match total dimension"
        );

        // Get kept subsystems
        let mut kept_subsystems = Vec::new();
        let mut kept_dims = Vec::new();
        for (i, &dim) in subsystem_dims.iter().enumerate() {
            if !trace_over.contains(&i) {
                kept_subsystems.push(i);
                kept_dims.push(dim);
            }
        }

        let new_dim: usize = kept_dims.iter().product();
        let mut result = Matrix::zeros((new_dim, new_dim));

        // Generate all possible index combinations
        for i in 0..self.dim {
            for j in 0..self.dim {
                let indices_i = self.multi_index(i, subsystem_dims);
                let indices_j = self.multi_index(j, subsystem_dims);

                // Check if traced indices match
                let mut trace_match = true;
                for &trace_idx in trace_over {
                    if indices_i[trace_idx] != indices_j[trace_idx] {
                        trace_match = false;
                        break;
                    }
                }

                if trace_match {
                    let kept_i: Vec<usize> =
                        kept_subsystems.iter().map(|&idx| indices_i[idx]).collect();
                    let kept_j: Vec<usize> =
                        kept_subsystems.iter().map(|&idx| indices_j[idx]).collect();

                    let flat_i = self.flat_index(&kept_i, &kept_dims);
                    let flat_j = self.flat_index(&kept_j, &kept_dims);

                    result[[flat_i, flat_j]] += self.matrix[[i, j]];
                }
            }
        }

        Density {
            num_qubits: 0,
            matrix: result,
            dim: new_dim,
        }
    }

    /// Get all basis state probabilities
    fn density_probability(&self) -> Vec<f64> {
        let matrix = self.matrix_as_ref();
        let dim = 1 << self.num_qubits;
        (0..dim).map(|i| matrix[[i, i]].re()).collect()
    }

    /// Calculate measurement probability for a qubit
    fn calculate_measurement_probability(&self, qubit: usize, outcome: bool) -> f64 {
        let matrix = self.matrix_as_ref();
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
        let matrix = self.matrix_as_ref().clone();
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

        self.matrix = new_matrix;
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

        self.matrix = new_matrix;
    }
}

impl QuantumState for Density {
    fn num_qubits(&self) -> usize {
        self.num_qubits()
    }

    fn apply(&mut self, u: Matrix<Complex>, _enumetared: GateType) {
        let current_state = self.matrix_as_ref();
        let u_dagger = ops::dagger(&u);
        let new_state = u.dot(current_state).dot(&u_dagger);
        self.set(new_state);
    }

    fn box_clone(&self) -> Box<dyn QuantumState> {
        Box::new(self.clone())
    }
}

impl BackendOperation for Density {
    fn expand_state(&mut self, qubit_state: Vector<Complex>) {
        assert_eq!(qubit_state.len(), 2, "Qubit state must be 2-dimensional");

        let density = Self::from_pure_statevec(&qubit_state);

        if self.is_empty() {
            *self = density;
        } else {
            *self = self.tensor_product(&density);
        }
    }

    fn measure_qubit(&mut self, qubit: usize) -> bool {
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

impl QuantumDebugger for Density {
    type StateType = Matrix<Complex>;
    type BitStateType = Complex;
    type QubitStateType = (f64, f64); 

    fn entire_state(&self, _filter_zero: bool) -> Self::StateType {
        self.matrix_as_ref().clone()
    }

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
        let matrix = self.matrix_as_ref();
        matrix[[state_index, state_index]]
    }

    fn qubit_state(&self, qubit: &QuantumBit) -> Self::QubitStateType {
        let num_qubits = self.num_qubits();
        let qubit = qubit.index();

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

/// Calculate number of qubits from matrix dimension
fn calculate_num_qubits(dim: usize) -> usize {
    if dim == 0 {
        return 0;
    }
    let num_qubits = (dim as f64).log2() as usize;
    assert_eq!(
        1 << num_qubits,
        dim,
        "Matrix dimension {} is not a power of 2",
        dim
    );
    num_qubits
}
