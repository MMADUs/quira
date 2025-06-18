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

use std::usize;

use rand::Rng;

use crate::{
    Complex, Matrix, QuantumGate,
    SingleQ::{Identity, PauliX, PauliY, PauliZ},
    constant::EPSILON,
    kernel::density::matrix::Density,
    utils::ops,
};

/// generic single pauli noise channel
/// applies a specific pauli gate (X, Y, or Z) with probability p
pub struct Flip {
    prob: f64,
    num_qubits: usize,
    kraus_ops: Vec<Matrix<Complex>>,
    channel_name: String,
}

impl Flip {
    /// create a new single pauli channel
    pub(crate) fn new_with_pauli(
        num_qubits: usize,
        prob: f64,
        pauli_op: Matrix<Complex>,
        name: &str,
    ) -> Self {
        assert!(num_qubits > 0, "number of qubits must be positive");
        assert!(
            (0.0..=1.0).contains(&prob),
            "{} probability must be 0.0 <= p <= 1.0",
            name
        );

        let kraus_ops = if prob.abs() < EPSILON {
            vec![Identity::new().unitary_matrix()]
        } else {
            Self::construct_kraus_operators(num_qubits, prob, pauli_op)
        };

        Self {
            prob,
            num_qubits,
            kraus_ops,
            channel_name: name.to_string(),
        }
    }

    /// create selective single pauli channel
    pub(crate) fn new_selective_with_pauli(
        total_qubits: usize,
        target_qubits: &[usize],
        prob: f64,
        pauli_op: Matrix<Complex>,
        name: &str,
    ) -> Self {
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
            vec![Identity::new().unitary_matrix()]
        } else {
            Self::construct_selective_kraus_operators(total_qubits, target_qubits, prob, pauli_op)
        };

        Self {
            prob,
            num_qubits: total_qubits,
            kraus_ops,
            channel_name: name.to_string(),
        }
    }

    /// return the error probability
    pub fn prob(&self) -> f64 {
        self.prob
    }

    /// return the channel name
    pub fn name(&self) -> &str {
        &self.channel_name
    }

    /// apply noise channel to density matrix
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

    /// sample error for each qubit independently
    pub fn sample_error<R: Rng>(&self, rng: &mut R) -> usize {
        let mut index = 0;
        for (i, _) in (0..self.num_qubits).enumerate() {
            if rng.random::<f64>() < self.prob {
                index |= 1 << i;
            }
        }
        index
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
    pub fn verify_completeness(&self) -> bool {
        let dim = 1 << self.num_qubits;
        let mut sum = Matrix::<Complex>::zeros((dim, dim));

        for kraus_op in &self.kraus_ops {
            let k_dagger = ops::dagger(&kraus_op);
            sum = sum + k_dagger.dot(kraus_op);
        }

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

    /// construct kraus operators for n-qubit single pauli channel
    fn construct_kraus_operators(
        num_qubits: usize,
        prob: f64,
        pauli_op: Matrix<Complex>,
    ) -> Vec<Matrix<Complex>> {
        let identity = Identity::new().unitary_matrix();
        let num_combinations = 1 << num_qubits; // 2^n combinations
        let mut kraus_ops = Vec::with_capacity(num_combinations);

        for pattern in 0..num_combinations {
            let mut result = Matrix::from_elem((1, 1), Complex::new(1.0, 0.0));
            let mut coeff = Complex::new(1.0, 0.0);

            for qubit in 0..num_qubits {
                let has_error = (pattern >> qubit) & 1 == 1;
                let gate = if has_error { &pauli_op } else { &identity };
                let prob_factor = if has_error { prob } else { 1.0 - prob };

                result = ops::kron(&result, gate);
                coeff *= Complex::new(prob_factor.sqrt(), 0.0);
            }

            kraus_ops.push(&result * coeff);
        }

        kraus_ops
    }

    /// construct selective kraus operators
    fn construct_selective_kraus_operators(
        total_qubits: usize,
        target_qubits: &[usize],
        prob: f64,
        pauli_op: Matrix<Complex>,
    ) -> Vec<Matrix<Complex>> {
        let identity = Identity::new().unitary_matrix();
        let num_combinations = 1 << target_qubits.len();
        let mut kraus_ops = Vec::with_capacity(num_combinations);

        for pattern in 0..num_combinations {
            let mut result = Matrix::from_elem((1, 1), Complex::new(1.0, 0.0));
            let mut coeff = Complex::new(1.0, 0.0);

            for qubit in 0..total_qubits {
                let gate = if let Some(target_idx) = target_qubits.iter().position(|&q| q == qubit)
                {
                    let has_error = (pattern >> target_idx) & 1 == 1;
                    let prob_factor = if has_error { prob } else { 1.0 - prob };
                    coeff *= Complex::new(prob_factor.sqrt(), 0.0);
                    if has_error { &pauli_op } else { &identity }
                } else {
                    &identity
                };

                result = ops::kron(&result, gate);
            }

            kraus_ops.push(&result * coeff);
        }

        kraus_ops
    }
}

/// bit flip noise channel - applies pauli-X gate with probability p
pub struct BitFlip;

impl BitFlip {
    /// new global bit flip channel
    pub fn new(num_qubits: usize, prob: f64) -> Flip {
        Flip::new_with_pauli(
            num_qubits,
            prob,
            PauliX::new(0).unitary_matrix(),
            "bit flip",
        )
    }

    /// new selective bit flip channel
    pub fn new_selective(total_qubits: usize, target_qubits: &[usize], prob: f64) -> Flip {
        Flip::new_selective_with_pauli(
            total_qubits,
            target_qubits,
            prob,
            PauliX::new(0).unitary_matrix(),
            "bit flip",
        )
    }
}

/// phase flip noise channel - applies pauli-Z gate with probability p
pub struct PhaseFlip;

impl PhaseFlip {
    /// new global phase flip channel
    pub fn new(num_qubits: usize, prob: f64) -> Flip {
        Flip::new_with_pauli(
            num_qubits,
            prob,
            PauliZ::new(0).unitary_matrix(),
            "phase flip",
        )
    }

    /// new selective phase flip channel
    pub fn new_selective(total_qubits: usize, target_qubits: &[usize], prob: f64) -> Flip {
        Flip::new_selective_with_pauli(
            total_qubits,
            target_qubits,
            prob,
            PauliZ::new(0).unitary_matrix(),
            "phase flip",
        )
    }
}

/// bit-phase flip noise channel - applies pauli-Y gate with probability p
pub struct BitPhaseFlip;

impl BitPhaseFlip {
    /// new global bit-phase flip channel
    pub fn new(num_qubits: usize, prob: f64) -> Flip {
        Flip::new_with_pauli(
            num_qubits,
            prob,
            PauliY::new(0).unitary_matrix(),
            "bit-phase flip",
        )
    }

    /// new selective bit-phase flip channel
    pub fn new_selective(total_qubits: usize, target_qubits: &[usize], prob: f64) -> Flip {
        Flip::new_selective_with_pauli(
            total_qubits,
            target_qubits,
            prob,
            PauliY::new(0).unitary_matrix(),
            "bit-phase flip",
        )
    }
}
