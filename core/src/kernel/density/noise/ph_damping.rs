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

use std::usize;

use ndarray::arr1;
use ndarray_linalg::{Scalar, Trace};
use rand::Rng;

use crate::{Complex, Matrix, constant::EPSILON, kernel::density::matrix::Density, ops};

/// Phase damping noise channel
/// Models pure dephasing without energy loss (T2* process)
/// Eliminates off-diagonal coherences while preserving populations
///
/// For single qubit: K0 = [[1, 0], [0, sqrt(1-γ)]], K1 = [[0, 0], [0, sqrt(γ)]]
/// For multi-qubit: Applied independently to each target qubit
pub struct PhaseDamping {
    prob: f64,
    num_qubits: usize,
    kraus_ops: Vec<Matrix<Complex>>,
}

impl PhaseDamping {
    /// Create new global phase damping channel (affects all qubits)
    pub fn new(num_qubits: usize, prob: f64) -> Self {
        assert!(num_qubits > 0, "number of qubits must be positive");
        assert!(
            (0.0..=1.0).contains(&prob),
            "phase damping probability must be 0.0 <= p <= 1.0"
        );

        let kraus_ops = if prob.abs() < EPSILON {
            vec![Matrix::eye(1 << num_qubits)]
        } else {
            let target_qubits = (0..num_qubits).collect::<Vec<usize>>().to_vec();
            Self::build_kraus_operators(&target_qubits, num_qubits, prob)
        };

        Self {
            prob,
            num_qubits,
            kraus_ops,
        }
    }

    /// Create new selective phase damping channel that applies to specific qubits
    pub fn new_selective(total_qubits: usize, target_qubits: Vec<usize>, prob: f64) -> Self {
        assert!(!target_qubits.is_empty(), "target qubits cannot be empty");
        assert!(
            target_qubits.iter().all(|&q| q < total_qubits),
            "target qubit indices must be less than total_qubits"
        );
        assert!(
            target_qubits.len() <= total_qubits,
            "cannot have more target qubits than total qubits",
        );

        let kraus_ops = if prob.abs() < EPSILON {
            vec![Matrix::eye(1 << total_qubits)]
        } else {
            Self::build_kraus_operators(&target_qubits, total_qubits, prob)
        };

        Self {
            prob,
            num_qubits: total_qubits,
            kraus_ops,
        }
    }

    /// Return the phase damping channel probability
    pub fn prob(&self) -> f64 {
        self.prob
    }

    /// Apply phase damping channel to a density matrix
    pub fn apply(&self, density: &Density) -> Density {
        let expected_dim = 1 << self.num_qubits;
        assert_eq!(
            density.dim(),
            expected_dim,
            "Expected density matrix of dim {}, got {}",
            expected_dim,
            density.dim()
        );
        density.apply_kraus_ops(&self.kraus_ops)
    }

    /// Sample which Kraus operator to apply based on probabilities
    pub fn sample_error<R: Rng>(&self, rng: &mut R) -> usize {
        let mut cumulative_prob = 0.0;
        let random_value: f64 = rng.random();

        for (idx, kraus_op) in self.kraus_ops.iter().enumerate() {
            // For phase damping, the probability is the trace of K†K
            let k_dagger = ops::dagger(kraus_op);
            let prob = k_dagger.dot(kraus_op).trace().unwrap().re();
            cumulative_prob += prob;

            if random_value <= cumulative_prob {
                return idx;
            }
        }
        self.kraus_ops.len() - 1
    }

    /// Apply a specific Kraus operator to the density matrix
    pub fn apply_sampled_error(&self, density: &Density, index: usize) -> Density {
        assert!(
            index < self.kraus_ops.len(),
            "Kraus operator index {} out of bounds",
            index
        );

        let k = &self.kraus_ops[index];
        let k_rho = k.dot(density.matrix_as_ref());
        let result = k_rho.dot(&ops::dagger(k));

        // Normalize the result
        let trace = result.trace().unwrap();
        assert!(
            trace.norm() > EPSILON,
            "Trace too close to zero; possible invalid Kraus op"
        );
        let normalized_result = result / trace;
        Density::new(normalized_result)
    }

    /// Sample a Kraus operator and apply it to the density matrix
    pub fn sample_and_apply<R: Rng>(&self, density: &Density, rng: &mut R) -> (usize, Density) {
        let index = self.sample_error(rng);
        let new_density = self.apply_sampled_error(density, index);
        (index, new_density)
    }

    /// Verify that Kraus operators satisfy completeness relation
    /// ∑ᵢ Kᵢ† Kᵢ = I
    pub fn verify_completeness(&self) -> bool {
        let dim = 1 << self.num_qubits;
        let mut sum = Matrix::<Complex>::zeros((dim, dim));

        for kraus_op in &self.kraus_ops {
            let k_dagger = ops::dagger(kraus_op);
            sum = sum + k_dagger.dot(kraus_op);
        }

        let identity = Matrix::<Complex>::eye(dim);
        let diff = sum - identity;

        diff.iter().all(|c| c.norm() < EPSILON)
    }

    /// Build Kraus operators for phase damping on specified qubits
    ///
    /// For n target qubits, this generates 2^n Kraus operators corresponding to
    /// all possible combinations of phase damping events on the target qubits.
    /// Each Kraus operator corresponds to a specific pattern of which qubits
    /// experience phase damping.
    fn build_kraus_operators(
        target_qubits: &[usize],
        total_qubits: usize,
        prob: f64,
    ) -> Vec<Matrix<Complex>> {
        let num_targets = target_qubits.len();
        let num_patterns = 1 << num_targets;
        let mut kraus_ops = Vec::new();

        for pattern in 0..num_patterns {
            let mut op = Matrix::<Complex>::from_elem((1, 1), Complex::new(1.0, 0.0)); // scalar identity

            for qubit in 0..total_qubits {
                let is_target = target_qubits.contains(&qubit);

                // If target, check the bit in the pattern
                let damping = if is_target {
                    let local_idx = target_qubits.iter().position(|&x| x == qubit).unwrap();
                    (pattern >> local_idx) & 1 == 1
                } else {
                    false
                };

                let local_op = if is_target {
                    if damping {
                        Matrix::from_diag(&arr1(&[
                            Complex::new(1.0, 0.0),
                            Complex::new(prob.sqrt(), 0.0),
                        ]))
                    } else {
                        Matrix::from_diag(&arr1(&[
                            Complex::new(1.0, 0.0),
                            Complex::new((1.0 - prob).sqrt(), 0.0),
                        ]))
                    }
                } else {
                    Matrix::eye(2)
                };

                op = ops::kron(&op, &local_op);
            }

            // Only add if not numerically zero
            let norm_sq: f64 = op.mapv(|c| c.norm_sqr()).sum();
            if norm_sq > EPSILON {
                kraus_ops.push(op);
            }
        }

        kraus_ops
    }
}
