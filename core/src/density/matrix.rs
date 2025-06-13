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

use ndarray_linalg::{Eigh, SVD, Scalar, Trace, UPLO};

use crate::{
    constant::EPSILON,
    types::{Complex, Matrix, Vector},
    utils::{eigen, ops},
};

pub struct Density {
    matrix: Matrix<Complex>,
    dim: usize,
}

impl Density {
    /// new density matrix from a matrix
    pub fn new(matrix: Matrix<Complex>) -> Self {
        let (rows, cols) = matrix.dim();
        assert_eq!(rows, cols, "matrix must be square");
        // new density
        let density = Self { matrix, dim: rows };
        assert!(density.is_valid(), "matrix is not a valid density matrix");
        density
    }

    /// new density matrix from pure state vector
    pub fn from_pure_statevec(state: &Vector<Complex>) -> Self {
        let dim = state.len();
        let mut matrix = Matrix::zeros((dim, dim));
        for i in 0..dim {
            for j in 0..dim {
                matrix[[i, j]] = state[i] * state[j].conj();
            }
        }
        Self { matrix, dim }
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
        let mut matrix = Matrix::zeros((dim, dim));
        for (state, &prob) in states.iter().zip(probs.iter()) {
            assert!(state.len() == dim, "all state must have equal dimension");
            let pure_density = Self::from_pure_statevec(state);
            matrix = matrix + pure_density.matrix * Complex::new(prob, 0.0);
        }
        Self { matrix, dim }
    }

    /// new density matrix in a maximally mixed state
    pub fn maximally_mixed(dim: usize) -> Self {
        let factor = Complex::new(1.0 / dim as f64, 0.0);
        let matrix = Matrix::eye(dim) * factor;
        Self { matrix, dim }
    }

    /// new density matrix in a thermal state
    pub fn thermal_state(hamiltonian: &Matrix<Complex>, beta: f64) -> Self {
        let (rows, cols) = hamiltonian.dim();
        assert_eq!(rows, cols, "hamiltonian must be square");
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
            matrix: rho,
            dim: rows,
        }
    }

    /// Unitary evolution U ρ U†
    pub fn evolve_unitary(&self, unitary: &Matrix<Complex>) -> Self {
        let u_rho = unitary.dot(self.matrix_as_ref());
        let u_dagger = ops::dagger(unitary);
        let evolved = u_rho.dot(&u_dagger);
        Self {
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
        Self { matrix, dim }
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
            matrix: result,
            dim: new_dim,
        }
    }

    /// Trace distance (1/2) * ||ρ - σ||₁
    pub fn trace_distance(&self, other: &Density) -> f64 {
        assert_eq!(self.dim, other.dim, "dimensions must match");
        let diff = &self.matrix - &other.matrix;
        let svd = diff.svd(true, true).expect("SVD failed");
        let singular_values = svd.1;
        let trace_norm: f64 = singular_values.sum();
        0.5 * trace_norm
    }

    /// Fidelity F(ρ, σ) = [Tr(√(√ρ σ √ρ))]^2
    pub fn fidelity(&self, other: &Density) -> f64 {
        assert_eq!(self.dim, other.dim, "dimensions must match");

        if self.is_pure() && other.is_pure() {
            // For pure states: F(ρ, σ) = Tr(ρσ) = |⟨ψ|φ⟩|²
            let product = self.matrix.dot(&other.matrix);
            let trace_val = product.trace().unwrap();
            trace_val.norm_sqr()
        } else {
            // General case: Uhlmann fidelity
            let sqrt_rho = eigen::matrix_sqrt(self.matrix_as_ref());
            let product = sqrt_rho.dot(&other.matrix).dot(&sqrt_rho);
            let eigenvals = eigen::eigen_values(&product);
            // Fidelity = (∑ sqrt(λᵢ))²
            let trace_sqrt = eigenvals
                .iter()
                .map(|val| val.re().max(0.0).sqrt()) // only real, non-negative eigenvalues
                .sum::<f64>();
            trace_sqrt * trace_sqrt
        }
    }

    /// Bures distance d_B(ρ, σ) = √(2(1 - √F(ρ, σ)))
    pub fn bures_distance(&self, other: &Density) -> f64 {
        let fidelity = self.fidelity(other);
        (2.0 * (1.0 - fidelity.sqrt())).sqrt()
    }

    /// Hilbert-Schmidt distance ||ρ - σ||_HS = √Tr((ρ - σ)²)
    pub fn hilbert_schmidt_distance(&self, other: &Density) -> f64 {
        assert_eq!(self.dim, other.dim, "dimensions must match");

        let diff = &self.matrix - &other.matrix;
        let diff_squared = diff.dot(&diff);
        diff_squared.diag().sum().re.sqrt()
    }

    /// Operator norm (spectral norm) ||ρ - σ||_∞
    pub fn operator_norm_distance(&self, other: &Density) -> f64 {
        assert_eq!(self.dim, other.dim, "dimensions must match");

        let diff = &self.matrix - &other.matrix;
        let eigenvals = eigen::eigen_values(&diff);
        eigenvals.iter().map(|&ev| ev.abs()).fold(0.0, f64::max)
    }
}
