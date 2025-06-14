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

use crate::{Complex, Matrix, constant::EPSILON, matrix::Density, ops};

/// dephasing noise channel
/// general random noise by applying Z-phase flips randomly
/// simulates phase errors due to environmental interactions, leading to loss of quantum coherence
pub struct Dephasing {
    prob: f64,
    num_qubits: usize,
    kraus_ops: Vec<Matrix<Complex>>,
}

impl Dephasing {
    /// new global dephasing channel
    pub fn new(num_qubits: usize, prob: f64) -> Self {
        assert!(num_qubits > 0, "number of qubits must be positive");
        assert!(
            (0.0..=1.0).contains(&prob),
            "dephasing probability must be 0.0 <= p <= 1.0"
        );
        // construct channel
        let kraus_ops = if prob.abs() < EPSILON {
            vec![Matrix::eye(1 << num_qubits)]
        } else {
            Self::global_kraus_operators(num_qubits, prob)
        };

        Self {
            prob,
            num_qubits,
            kraus_ops,
        }
    }

    /// new mixed dephasing channel that applies to specific qubits
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
        // construct channel
        let kraus_ops = if prob.abs() < EPSILON {
            vec![Matrix::eye(1 << total_qubits)]
        } else {
            Self::selective_kraus_operators(target_qubits, total_qubits, prob)
        };

        Self {
            prob,
            num_qubits: total_qubits,
            kraus_ops,
        }
    }

    /// return the dephasing channel probabilities
    pub fn prob(&self) -> f64 {
        self.prob
    }

    /// apply dephasing channel to a density matrix
    pub fn apply(&self, density: &Density) -> Density {
        let expected_dim = 1 << self.num_qubits;
        assert_eq!(
            density.dim(),
            expected_dim,
            "Density matrix dimension {} doesn't match expected {} for {} qubits",
            density.dim(),
            expected_dim,
            self.num_qubits
        );

        density.apply_kraus_ops(&self.kraus_ops)
    }

    pub fn sample_error<R: Rng>(&self, rng: &mut R) -> usize {
        let mut cummulative_prob = 0.0;
        let random_value: f64 = rng.random();

        for (idx, kraus_op) in self.kraus_ops.iter().enumerate() {
            let prob = kraus_op.mapv(|x| x.norm_sqr()).sum();
            cummulative_prob += prob;

            if random_value <= cummulative_prob {
                return idx;
            }
        }
        self.kraus_ops.len() - 1
    }

    /// Verify that Kraus operators satisfy completeness relation
    /// ∑ᵢ Kᵢ† Kᵢ = I
    pub fn verify_completeness(&self) -> bool {
        let dim = 1 << self.num_qubits;
        let mut sum = Matrix::<Complex>::zeros((dim, dim));
        // compute sum
        for kraus_op in &self.kraus_ops {
            let k_dagger = ops::dagger(kraus_op);
            sum = sum + k_dagger.dot(kraus_op);
        }
        // check if sum is approximately identity
        let identity = Matrix::<Complex>::eye(dim);
        let diff = sum - identity;
        // check every element close to zero
        for i in 0..dim {
            for j in 0..dim {
                if diff[[i, j]].norm() > EPSILON {
                    return false;
                }
            }
        }
        true
    }

    /// apply singe qubit z operation to multi qubit operator
    fn apply_z_phase_to_qubit(
        num_qubits: usize,
        operator: &Matrix<Complex>,
        target_qubit: usize,
    ) -> Matrix<Complex> {
        let mut result = operator.clone();
        let total_dim = 1 << num_qubits;
        // for each computational basis state, apply phase flip if qubit is |1⟩
        for i in 0..total_dim {
            for j in 0..total_dim {
                // check if target qubit is |1⟩ in states i and j
                let bit_i = (i >> target_qubit) & 1;
                let bit_j = (j >> target_qubit) & 1;
                // apply phase: (-1)^bit for each state
                let sign_i = if bit_i == 1 { -1.0 } else { 1.0 };
                let sign_j = if bit_j == 1 { -1.0 } else { 1.0 };
                let phase = sign_i * sign_j;
                result[[i, j]] = operator[[i, j]] * Complex::new(phase, 0.0);
            }
        }
        result
    }

    /// construct the kraus operation for dephasing channel
    fn global_kraus_operators(num_qubits: usize, prob: f64) -> Vec<Matrix<Complex>> {
        let mut kraus_ops = Vec::new();
        let total_dim = 1 << num_qubits;

        for pattern in 0..total_dim as i32 {
            let num_dephasings = pattern.count_ones() as usize;

            let amplitude = if num_dephasings == 0 {
                ((1.0 - prob).powi(num_qubits as i32)).sqrt()
            } else {
                let p_power = prob.powi(num_dephasings as i32);
                let one_minus_p_power = (1.0 - prob).powi((num_qubits - num_dephasings) as i32);
                (p_power * one_minus_p_power).sqrt()
            };

            if amplitude.abs() < EPSILON {
                continue;
            }

            let mut operator = Matrix::eye(total_dim);

            for qubit in 0..num_qubits {
                if (pattern >> qubit) & 1 == 1 {
                    operator = Self::apply_z_phase_to_qubit(num_qubits, &operator, qubit);
                }
            }

            operator *= Complex::new(amplitude, 0.0);
            kraus_ops.push(operator);
        }
        kraus_ops
    }

    fn selective_kraus_operators(
        target_qubits: &[usize],
        total_qubits: usize,
        prob: f64,
    ) -> Vec<Matrix<Complex>> {
        let mut kraus_ops = Vec::new();
        let total_dim = 1 << total_qubits;
        let num_target_qubits = target_qubits.len();
        // construct all patterns for the target qubits
        for pattern in 0..(1 << num_target_qubits) as i32 {
            let num_dephasings = pattern.count_ones() as usize;

            let amplitude = if num_dephasings == 0 {
                ((1.0 - prob).powi(num_target_qubits as i32)).sqrt()
            } else {
                let p_power = prob.powi(num_dephasings as i32);
                let one_minus_p_power =
                    (1.0 - prob).powi((num_target_qubits - num_dephasings) as i32);
                (p_power * one_minus_p_power).sqrt()
            };

            if amplitude.abs() < EPSILON {
                continue;
            }

            let mut operator = Matrix::eye(total_dim);

            for (local_idx, &global_qubit) in target_qubits.iter().enumerate() {
                if (pattern >> local_idx) & 1 == 1 {
                    operator = Self::apply_z_phase_to_qubit(total_qubits, &operator, global_qubit);
                }
            }

            operator *= Complex::new(amplitude, 0.0);
            kraus_ops.push(operator);
        }
        kraus_ops
    }
}
