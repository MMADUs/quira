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

use rand::Rng;

use crate::{
    Complex, Matrix, QuantumGate,
    SingleQ::{Identity, PauliX, PauliY, PauliZ},
    constant::EPSILON,
    matrix::Density,
    utils::ops,
};

/// depolarizing noise channel  
/// general random noise by applying Pauli X, Y, or Z gates uniformly  
/// simulates complete loss of information due to uniform noise on all axes  
pub struct Depolarizing {
    prob: f64,
    num_qubits: usize,
    kraus_ops: Vec<Matrix<Complex>>,
}

impl Depolarizing {
    /// new global depolarizing channel
    pub fn new(num_qubits: usize, prob: f64) -> Self {
        assert!(num_qubits > 0, "number of qubits must be positive");
        assert!(
            (0.0..=1.0).contains(&prob),
            "depolarizing probability must be 0.0 <= p <= 1.0"
        );
        // construct channel
        let kraus_ops = if prob.abs() < EPSILON {
            // Only Identity operator
            vec![Identity::new().unitary_matrix()]
        } else {
            let pauli_ops = Self::global_pauli_operators(num_qubits);
            Self::construct_kraus_operators(num_qubits, prob, &pauli_ops)
        };

        Self {
            prob,
            num_qubits,
            kraus_ops,
        }
    }

    /// new mixed depolarizing channel that applies to specific qubits
    pub fn new_selective(total_qubits: usize, target_qubits: &[usize], prob: f64) -> Self {
        assert!(!target_qubits.is_empty(), "target qubits cannot be empty");
        assert!(
            target_qubits.iter().all(|&q| q < total_qubits),
            "target qubit indices must be less than total_qubits"
        );
        assert!(
            target_qubits.len() <= total_qubits,
            "cannot have more target qubits than total qubits",
        );
        // construct mixed channel
        let kraus_ops = if prob.abs() < EPSILON {
            // Only Identity
            vec![Identity::new().unitary_matrix()]
        } else {
            let num_target_qubits = target_qubits.len();
            let pauli_ops = Self::selective_pauli_operators(total_qubits, target_qubits);
            Self::construct_kraus_operators(num_target_qubits, prob, &pauli_ops)
        };

        Self {
            prob,
            num_qubits: total_qubits,
            kraus_ops,
        }
    }

    /// return the depolarizing channel probabilities
    pub fn prob(&self) -> f64 {
        self.prob
    }

    /// apply depolarizing channel to density
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

    /// samples a pauli operator index according to a depolarizing channel.
    pub fn sample_error<R: Rng>(&self, rng: &mut R) -> usize {
        let r: f64 = rng.random();
        let num_paulis = 4_usize.pow(self.num_qubits as u32);

        if r < 1.0 - self.prob {
            0
        } else {
            let span = num_paulis - 1;
            let offset = r - (1.0 - self.prob);
            let fraction = offset / self.prob;
            1 + (fraction * span as f64).floor() as usize
        }
    }

    /// apply a specific kraus operator to the density matrix
    pub fn apply_sampled_error(&self, density: &Density, index: usize) -> Density {
        assert!(
            index < self.kraus_ops.len(),
            "kraus operator index {} out of bounds",
            index
        );

        let k = &self.kraus_ops[index];
        let k_rho = k.dot(density.matrix_as_ref());
        let result = k_rho.dot(&ops::dagger(k));

        Density::new(result)
    }

    /// samples a kraus operator and applies it to the density matrix
    pub fn sample_and_apply<R: Rng>(&self, density: &Density, rng: &mut R) -> (usize, Density) {
        let index = self.sample_error(rng);
        let new_density = self.apply_sampled_error(density, index);
        (index, new_density)
    }

    /// check if kraus operators satisfy the completeness relation
    /// Σᵢ Kᵢ† Kᵢ = I
    pub fn verify_completeness(&self) -> bool {
        let dim = 1 << self.num_qubits;
        let mut sum = Matrix::<Complex>::zeros((dim, dim));
        // compute sum
        for kraus_op in &self.kraus_ops {
            let k_dagger = ops::dagger(&kraus_op);
            sum = sum + k_dagger.dot(kraus_op);
        }
        // check satisfication
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

    /// construct pauli operators for n-qubit depolarizing channel
    fn global_pauli_operators(num_qubits: usize) -> Vec<Matrix<Complex>> {
        let single_paulis = vec![
            Identity::new().unitary_matrix(),
            PauliX::new(0).unitary_matrix(),
            PauliY::new(0).unitary_matrix(),
            PauliZ::new(0).unitary_matrix(),
        ];
        let num_paulis = 4_usize.pow(num_qubits as u32);
        let mut n_paulis = Vec::with_capacity(num_paulis);
        // generate all combination using base-4 representation
        for i in 0..num_paulis {
            let mut pauli_indices = Vec::with_capacity(num_qubits);
            let mut temp_i = i;
            // convert index to base-4 representation
            for _ in 0..num_qubits {
                pauli_indices.push(temp_i % 4);
                temp_i /= 4;
            }
            pauli_indices.reverse();
            // compute tensor product of pauli operators
            let mut result = single_paulis[pauli_indices[0]].clone();
            for &pauli_idx in &pauli_indices[1..] {
                result = ops::kron(&result, &single_paulis[pauli_idx]);
            }
            n_paulis.push(result);
        }
        n_paulis
    }

    /// construct pauli operators for selective qubit depolarizing channel
    fn selective_pauli_operators(
        total_qubits: usize,
        target_qubits: &[usize],
    ) -> Vec<Matrix<Complex>> {
        let single_paulis = vec![
            Identity::new().unitary_matrix(),
            PauliX::new(0).unitary_matrix(),
            PauliY::new(0).unitary_matrix(),
            PauliZ::new(0).unitary_matrix(),
        ];
        let identity = &single_paulis[0];
        let num_target_qubits = target_qubits.len();
        let num_target_paulis = 4_usize.pow(num_target_qubits as u32);
        let mut selective_paulis = Vec::with_capacity(num_target_paulis);
        // generate all combination using base-4 representation
        for i in 0..num_target_paulis {
            let mut pauli_indices = Vec::with_capacity(num_target_qubits);
            let mut temp_i = i;
            // convert index to base-4 representation
            for _ in 0..num_target_qubits {
                pauli_indices.push(temp_i % 4);
                temp_i /= 4;
            }
            pauli_indices.reverse();
            // apply pauli operation for specific qubit
            let mut result = Matrix::from_elem((1, 1), Complex::new(1.0, 0.0));
            let mut target_idx = 0;
            for qubit in 0..total_qubits {
                let pauli_for_qubit = if target_qubits.contains(&qubit) {
                    let pauli_op = &single_paulis[pauli_indices[target_idx]];
                    target_idx += 1;
                    pauli_op
                } else {
                    identity
                };
                // compute tensor
                result = ops::kron(&result, pauli_for_qubit);
            }
            selective_paulis.push(result);
        }
        selective_paulis
    }

    /// construct kraus operators for n-qubit depolarizing channel
    fn construct_kraus_operators(
        num_qubits: usize,
        prob: f64,
        pauli_ops: &[Matrix<Complex>],
    ) -> Vec<Matrix<Complex>> {
        let num_paulis = 4_usize.pow(num_qubits as u32);
        let mut kraus_ops = Vec::with_capacity(num_paulis);
        // Identity term: sqrt(1-p) * I
        let coeff_identity = Complex::new((1.0 - prob).sqrt(), 0.0);
        kraus_ops.push(&pauli_ops[0] * coeff_identity);
        // Pauli terms: sqrt(p/(4^n - 1)) * Pᵢ for the (4^n - 1) non-identity Pauli operators
        let coeff_pauli = Complex::new((prob / (num_paulis - 1) as f64).sqrt(), 0.0);
        for i in 1..num_paulis {
            kraus_ops.push(&pauli_ops[i] * coeff_pauli);
        }
        kraus_ops
    }
}
