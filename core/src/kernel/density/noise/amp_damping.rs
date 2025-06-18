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

use ndarray_linalg::{Scalar, Trace};
use rand::Rng;

use crate::{
    Complex, Matrix, QuantumGate, SingleQ::Identity, constant::EPSILON,
    kernel::density::matrix::Density, utils::ops,
};

/// Amplitude damping noise channel
/// Models energy loss in quantum systems (e.g., spontaneous emission)
/// Simulates the decay of |1⟩ to |0⟩ state with probability γ
pub struct AmplitudeDamping {
    gamma: f64,
    num_qubits: usize,
    kraus_ops: Vec<Matrix<Complex>>,
}

impl AmplitudeDamping {
    /// Create new global amplitude damping channel
    /// γ (gamma): damping parameter, probability of energy loss
    pub fn new(num_qubits: usize, gamma: f64) -> Self {
        assert!(num_qubits > 0, "number of qubits must be positive");
        assert!(
            (0.0..=1.0).contains(&gamma),
            "amplitude damping parameter must be 0.0 <= γ <= 1.0"
        );

        // Construct Kraus operators
        let kraus_ops = if gamma.abs() < EPSILON {
            // No damping, only identity
            let dim = 1 << num_qubits;
            let identity = Matrix::eye(dim).mapv(|x| Complex::new(x, 0.0));
            vec![identity]
        } else {
            Self::global_kraus_operators(num_qubits, gamma)
        };

        Self {
            gamma,
            num_qubits,
            kraus_ops,
        }
    }

    /// Create new selective amplitude damping channel that applies to specific qubits
    pub fn new_selective(total_qubits: usize, target_qubits: &[usize], gamma: f64) -> Self {
        assert!(!target_qubits.is_empty(), "target qubits cannot be empty");
        assert!(
            target_qubits.iter().all(|&q| q < total_qubits),
            "target qubit indices must be less than total_qubits"
        );
        assert!(
            target_qubits.len() <= total_qubits,
            "cannot have more target qubits than total qubits",
        );
        assert!(
            (0.0..=1.0).contains(&gamma),
            "amplitude damping parameter must be 0.0 <= γ <= 1.0"
        );

        // Construct selective channel
        let kraus_ops = if gamma.abs() < EPSILON {
            // No damping, only identity
            let dim = 1 << total_qubits;
            let identity = Matrix::eye(dim).mapv(|x| Complex::new(x, 0.0));
            vec![identity]
        } else {
            Self::selective_kraus_operators(total_qubits, target_qubits, gamma)
        };

        Self {
            gamma,
            num_qubits: total_qubits,
            kraus_ops,
        }
    }

    /// Return the amplitude damping parameter
    pub fn gamma(&self) -> f64 {
        self.gamma
    }

    /// Apply amplitude damping channel to density matrix
    pub fn apply(&self, density: &Density) -> Density {
        let expected_dim = 1 << self.num_qubits;
        assert_eq!(
            density.dim(),
            expected_dim,
            "density matrix dimension {} doesn't match expected dimension {} for {} qubit channel",
            density.dim(),
            expected_dim,
            self.num_qubits,
        );
        density.apply_kraus_ops(&self.kraus_ops)
    }

    /// Sample a Kraus operator index and return the collapsed state
    pub fn sample_kraus<R: Rng>(&self, density: &Density, rng: &mut R) -> (usize, Density) {
        let mut probabilities = Vec::with_capacity(self.kraus_ops.len());

        for kraus in &self.kraus_ops {
            let k_rho = kraus.dot(density.matrix_as_ref()).dot(&ops::dagger(kraus));
            let prob = k_rho.trace().unwrap_or(Complex::new(0.0, 0.0)).re(); // trace is real
            probabilities.push(prob);
        }

        // Normalize (in case of numerical drift)
        let sum: f64 = probabilities.iter().sum();
        let normalized: Vec<f64> = probabilities.iter().map(|p| p / sum).collect();

        // Sample index
        let mut acc = 0.0;
        let r = rng.random::<f64>();
        for (i, &p) in normalized.iter().enumerate() {
            acc += p;
            if r < acc {
                // Return the index and the resulting collapsed density
                let k = &self.kraus_ops[i];
                let new_rho = k.dot(density.matrix_as_ref()).dot(&ops::dagger(k));
                return (i, Density::new(new_rho));
            }
        }

        // Fallback to last index in case of float rounding
        let k = &self.kraus_ops[self.kraus_ops.len() - 1];
        let new_rho = k.dot(density.matrix_as_ref()).dot(&ops::dagger(k));
        (self.kraus_ops.len() - 1, Density::new(new_rho))
    }

    /// Check if Kraus operators satisfy the completeness relation
    /// Σᵢ Kᵢ† Kᵢ = I
    pub fn verify_completeness(&self) -> bool {
        let dim = 1 << self.num_qubits;
        let mut sum = Matrix::<Complex>::zeros((dim, dim));

        // Compute sum of Kᵢ† Kᵢ
        for kraus_op in &self.kraus_ops {
            let k_dagger = ops::dagger(kraus_op);
            sum = sum + k_dagger.dot(kraus_op);
        }

        // Check if sum equals identity matrix
        for i in 0..dim {
            for j in 0..dim {
                let expected = if i == j {
                    Complex::new(1.0, 0.0)
                } else {
                    Complex::new(0.0, 0.0)
                };
                if (sum[[i, j]] - expected).norm() > EPSILON {
                    return false;
                }
            }
        }
        true
    }

    /// Create single-qubit amplitude damping Kraus operators
    fn single_qubit_kraus_ops(gamma: f64) -> (Matrix<Complex>, Matrix<Complex>) {
        // K₀ = |0⟩⟨0| + √(1-γ)|1⟩⟨1|
        let mut k0 = Matrix::<Complex>::zeros((2, 2));
        k0[[0, 0]] = Complex::new(1.0, 0.0);
        k0[[1, 1]] = Complex::new((1.0 - gamma).sqrt(), 0.0);

        // K₁ = √γ|0⟩⟨1|
        let mut k1 = Matrix::<Complex>::zeros((2, 2));
        k1[[0, 1]] = Complex::new(gamma.sqrt(), 0.0);

        (k0, k1)
    }

    /// Construct Kraus operators for global amplitude damping channel
    fn global_kraus_operators(num_qubits: usize, gamma: f64) -> Vec<Matrix<Complex>> {
        let (k0, k1) = Self::single_qubit_kraus_ops(gamma);
        let num_ops = 2_usize.pow(num_qubits as u32);
        let mut kraus_ops = Vec::with_capacity(num_ops);

        // Generate all combinations using binary representation
        for i in 0..num_ops {
            let mut result = Matrix::from_elem((1, 1), Complex::new(1.0, 0.0));

            // Build tensor product based on binary representation
            for qubit in 0..num_qubits {
                let use_k1 = (i >> qubit) & 1 == 1;
                let kraus_op = if use_k1 { &k1 } else { &k0 };
                result = ops::kron(&result, kraus_op);
            }

            kraus_ops.push(result);
        }

        kraus_ops
    }

    /// Construct Kraus operators for selective amplitude damping channel
    fn selective_kraus_operators(
        total_qubits: usize,
        target_qubits: &[usize],
        gamma: f64,
    ) -> Vec<Matrix<Complex>> {
        let (k0, k1) = Self::single_qubit_kraus_ops(gamma);
        let identity = Identity::new().unitary_matrix();

        let num_target_qubits = target_qubits.len();
        let num_ops = 2_usize.pow(num_target_qubits as u32);
        let mut kraus_ops = Vec::with_capacity(num_ops);

        // Generate all combinations for target qubits
        for i in 0..num_ops {
            let mut result = Matrix::from_elem((1, 1), Complex::new(1.0, 0.0));
            let mut target_idx = 0;

            // Build tensor product for all qubits
            for qubit in 0..total_qubits {
                let kraus_for_qubit = if target_qubits.contains(&qubit) {
                    let use_k1 = (i >> target_idx) & 1 == 1;
                    target_idx += 1;
                    if use_k1 { &k1 } else { &k0 }
                } else {
                    &identity
                };

                result = ops::kron(&result, kraus_for_qubit);
            }

            kraus_ops.push(result);
        }

        kraus_ops
    }
}
