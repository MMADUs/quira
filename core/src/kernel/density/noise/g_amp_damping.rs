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
    Complex, Matrix, QuantumGate,
    SingleQ::Identity,
    constant::{BOLTZMANN, EPSILON, INF, PLANCK},
    kernel::density::matrix::Density,
    utils::ops,
};

/// Generalized amplitude damping noise channel
/// Models both energy loss (T₁ decay) and thermal effects
/// Includes both spontaneous emission and thermal excitation
///
/// This is equivalent to Qiskit's GeneralizedAmplitudeDampingError and
/// Cirq's GeneralizedAmplitudeDampingChannel
pub struct GeneralizedAmplitudeDamping {
    gamma: f64,
    p_excited: f64,
    num_qubits: usize,
    kraus_ops: Vec<Matrix<Complex>>,
}

impl GeneralizedAmplitudeDamping {
    /// Create new global generalized amplitude damping channel
    ///
    /// Parameters:
    /// - num_qubits: number of qubits in the system
    /// - gamma: damping parameter (0 ≤ γ ≤ 1), probability of energy decay
    /// - p_excited: thermal excitation probability (0 ≤ p ≤ 1)
    ///   - p = 0: pure amplitude damping (T=0)
    ///   - p > 0: finite temperature effects
    ///   - p = n̄/(n̄+1) where n̄ is mean thermal photon number
    pub fn new(num_qubits: usize, gamma: f64, p_excited: f64) -> Self {
        assert!(num_qubits > 0, "number of qubits must be positive");
        assert!(
            (0.0..=1.0).contains(&gamma),
            "damping parameter must be 0.0 <= γ <= 1.0"
        );
        assert!(
            (0.0..=1.0).contains(&p_excited),
            "excited state probability must be 0.0 <= p <= 1.0"
        );

        // Construct Kraus operators
        let kraus_ops = if gamma.abs() < EPSILON {
            // No damping, only identity
            let dim = 1 << num_qubits;
            let identity = Matrix::eye(dim).mapv(|x| Complex::new(x, 0.0));
            vec![identity]
        } else {
            Self::global_kraus_operators(num_qubits, gamma, p_excited)
        };

        Self {
            gamma,
            p_excited,
            num_qubits,
            kraus_ops,
        }
    }

    /// Create new selective generalized amplitude damping channel
    /// Applies damping only to specified target qubits
    pub fn new_selective(
        total_qubits: usize,
        target_qubits: &[usize],
        gamma: f64,
        p_excited: f64,
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
        assert!(
            (0.0..=1.0).contains(&gamma),
            "damping parameter must be 0.0 <= γ <= 1.0"
        );
        assert!(
            (0.0..=1.0).contains(&p_excited),
            "excited state probability must be 0.0 <= p <= 1.0"
        );

        // Construct selective channel
        let kraus_ops = if gamma.abs() < EPSILON {
            // No damping, only identity
            let dim = 1 << total_qubits;
            let identity = Matrix::eye(dim).mapv(|x| Complex::new(x, 0.0));
            vec![identity]
        } else {
            Self::selective_kraus_operators(total_qubits, target_qubits, gamma, p_excited)
        };

        Self {
            gamma,
            p_excited,
            num_qubits: total_qubits,
            kraus_ops,
        }
    }

    /// Create from temperature (Kelvin) and frequency
    /// Convenience constructor for thermal channels
    ///
    /// Parameters:
    /// - num_qubits: number of qubits
    /// - gamma: damping parameter
    /// - temperature: temperature in Kelvin
    /// - frequency: transition frequency in Hz
    pub fn from_thermal(num_qubits: usize, gamma: f64, temperature: f64, frequency: f64) -> Self {
        assert!(temperature > 0.0, "temperature must be positive");
        assert!(frequency > 0.0, "frequency must be positive");

        // Calculate thermal excitation probability
        // p = n̄/(n̄+1) where n̄ = 1/(exp(ℏω/kT) - 1)
        let hf_kt = (PLANCK * frequency) / (BOLTZMANN * temperature);
        let n_bar = 1.0 / (hf_kt.exp() - 1.0);
        let p_excited = n_bar / (n_bar + 1.0);

        Self::new(num_qubits, gamma, p_excited)
    }

    /// Return the damping parameter
    pub fn gamma(&self) -> f64 {
        self.gamma
    }

    /// Return the thermal excitation probability
    pub fn p_excited(&self) -> f64 {
        self.p_excited
    }

    /// Get effective temperature parameter
    /// Returns n̄ (mean thermal photon number)
    pub fn thermal_photon_number(&self) -> f64 {
        if self.p_excited >= 1.0 {
            INF
        } else {
            self.p_excited / (1.0 - self.p_excited)
        }
    }

    /// Apply generalized amplitude damping channel to density matrix
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
            let prob = k_rho.trace().unwrap_or(Complex::new(0.0, 0.0)).re();
            probabilities.push(prob.max(0.0)); // Ensure non-negative
        }

        // Normalize probabilities
        let sum: f64 = probabilities.iter().sum();
        if sum < EPSILON {
            // Fallback to uniform distribution if all probabilities are zero
            probabilities = vec![1.0 / self.kraus_ops.len() as f64; self.kraus_ops.len()];
        } else {
            probabilities = probabilities.iter().map(|p| p / sum).collect();
        }

        // Sample index
        let mut acc = 0.0;
        let r = rng.random::<f64>();
        for (i, &p) in probabilities.iter().enumerate() {
            acc += p;
            if r < acc {
                let k = &self.kraus_ops[i];
                let new_rho = k.dot(density.matrix_as_ref()).dot(&ops::dagger(k));
                let prob = probabilities[i];
                let normalized_rho = if prob > EPSILON {
                    new_rho.mapv(|x| x / Complex::new(prob, 0.0))
                } else {
                    new_rho
                };
                return (i, Density::new(normalized_rho));
            }
        }

        // Fallback to last index
        let k = &self.kraus_ops[self.kraus_ops.len() - 1];
        let new_rho = k.dot(density.matrix_as_ref()).dot(&ops::dagger(k));
        let prob = probabilities[probabilities.len() - 1];
        let normalized_rho = if prob > EPSILON {
            new_rho.mapv(|x| x / Complex::new(prob, 0.0))
        } else {
            new_rho
        };
        (self.kraus_ops.len() - 1, Density::new(normalized_rho))
    }

    /// Check if Kraus operators satisfy the completeness relation
    /// Σᵢ Kᵢ† Kᵢ = I
    pub fn verify_completeness(&self) -> bool {
        let dim = 1 << self.num_qubits;
        let mut sum = Matrix::<Complex>::zeros((dim, dim));

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

    /// Create single-qubit generalized amplitude damping Kraus operators
    ///
    /// The four Kraus operators are:
    /// K₀ = √(1-p) * (|0⟩⟨0| + √(1-γ)|1⟩⟨1|)  // No decay, ground state
    /// K₁ = √(1-p) * √γ|0⟩⟨1|                   // Decay to ground state  
    /// K₂ = √p * (√(1-γ)|0⟩⟨0| + |1⟩⟨1|)       // No decay, excited state
    /// K₃ = √p * √γ|1⟩⟨0|                       // Thermal excitation
    fn single_qubit_kraus_ops(gamma: f64, p_excited: f64) -> [Matrix<Complex>; 4] {
        let sqrt_gamma = gamma.sqrt();
        let sqrt_1_minus_gamma = (1.0 - gamma).sqrt();
        let sqrt_p = p_excited.sqrt();
        let sqrt_1_minus_p = (1.0 - p_excited).sqrt();

        // K₀: No decay from ground state distribution
        let mut k0 = Matrix::<Complex>::zeros((2, 2));
        k0[[0, 0]] = Complex::new(sqrt_1_minus_p, 0.0);
        k0[[1, 1]] = Complex::new(sqrt_1_minus_p * sqrt_1_minus_gamma, 0.0);

        // K₁: Decay to ground state
        let mut k1 = Matrix::<Complex>::zeros((2, 2));
        k1[[0, 1]] = Complex::new(sqrt_1_minus_p * sqrt_gamma, 0.0);

        // K₂: No decay from excited state distribution
        let mut k2 = Matrix::<Complex>::zeros((2, 2));
        k2[[0, 0]] = Complex::new(sqrt_p * sqrt_1_minus_gamma, 0.0);
        k2[[1, 1]] = Complex::new(sqrt_p, 0.0);

        // K₃: Thermal excitation
        let mut k3 = Matrix::<Complex>::zeros((2, 2));
        k3[[1, 0]] = Complex::new(sqrt_p * sqrt_gamma, 0.0);

        [k0, k1, k2, k3]
    }

    /// Construct Kraus operators for global generalized amplitude damping
    fn global_kraus_operators(
        num_qubits: usize,
        gamma: f64,
        p_excited: f64,
    ) -> Vec<Matrix<Complex>> {
        let single_ops = Self::single_qubit_kraus_ops(gamma, p_excited);
        let num_ops = 4_usize.pow(num_qubits as u32);
        let mut kraus_ops = Vec::with_capacity(num_ops);

        // Generate all combinations using base-4 representation
        for i in 0..num_ops {
            let mut result = Matrix::from_elem((1, 1), Complex::new(1.0, 0.0));
            let mut temp_i = i;

            // Build tensor product based on base-4 representation
            for _ in 0..num_qubits {
                let op_idx = temp_i % 4;
                temp_i /= 4;
                result = ops::kron(&result, &single_ops[op_idx]);
            }
            kraus_ops.push(result);
        }
        kraus_ops
    }

    /// Construct Kraus operators for selective generalized amplitude damping
    fn selective_kraus_operators(
        total_qubits: usize,
        target_qubits: &[usize],
        gamma: f64,
        p_excited: f64,
    ) -> Vec<Matrix<Complex>> {
        let single_ops = Self::single_qubit_kraus_ops(gamma, p_excited);
        let identity = Identity::new().unitary_matrix();

        let num_target_qubits = target_qubits.len();
        let num_ops = 4_usize.pow(num_target_qubits as u32);
        let mut kraus_ops = Vec::with_capacity(num_ops);

        // Generate all combinations for target qubits
        for i in 0..num_ops {
            let mut result = Matrix::from_elem((1, 1), Complex::new(1.0, 0.0));
            let mut temp_i = i;
            let mut target_idx = 0;

            // Build tensor product for all qubits
            for qubit in 0..total_qubits {
                let op_for_qubit = if target_qubits.contains(&qubit) {
                    let op_idx = temp_i % 4;
                    if target_idx == 0 {
                        temp_i /= 4;
                    }
                    target_idx += 1;
                    &single_ops[op_idx]
                } else {
                    &identity
                };
                result = ops::kron(&result, op_for_qubit);
            }
            kraus_ops.push(result);
        }
        kraus_ops
    }
}
