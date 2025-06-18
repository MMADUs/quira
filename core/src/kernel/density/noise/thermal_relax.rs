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

/// Thermal relaxation noise channel
/// Models realistic quantum decoherence including both T₁ (energy relaxation) and T₂ (dephasing)
///
/// This combines:
/// - Amplitude damping (T₁ decay): |1⟩ → |0⟩ transitions
/// - Pure dephasing (T₂* effects): phase randomization without energy loss
/// - Thermal excitation: |0⟩ → |1⟩ transitions at finite temperature
///
/// This is equivalent to Qiskit's ThermalRelaxationError and Cirq's ThermalNoiseModel
pub struct ThermalRelaxation {
    t1: f64,
    t2: f64,
    gate_time: f64,
    temperature: f64,
    frequency: f64,
    num_qubits: usize,
    kraus_ops: Vec<Matrix<Complex>>,
}

impl ThermalRelaxation {
    /// Create new global thermal relaxation channel
    ///
    /// Parameters:
    /// - num_qubits: number of qubits in the system
    /// - t1: T₁ time (energy relaxation time) in seconds
    /// - t2: T₂ time (dephasing time) in seconds, must satisfy T₂ ≤ 2T₁
    /// - gate_time: gate execution time in seconds
    /// - temperature: temperature in Kelvin (default: 0.0 for ground state)
    /// - frequency: qubit transition frequency in Hz (default: 5e9 Hz ≈ 5 GHz)
    pub fn new(
        num_qubits: usize,
        t1: f64,
        t2: f64,
        gate_time: f64,
        temperature: Option<f64>,
        frequency: Option<f64>,
    ) -> Self {
        assert!(num_qubits > 0, "number of qubits must be positive");
        assert!(t1 > 0.0, "T₁ time must be positive");
        assert!(t2 > 0.0, "T₂ time must be positive");
        assert!(
            t2 <= 2.0 * t1 + EPSILON,
            "T₂ must be ≤ 2T₁ (T₂ = {} > 2T₁ = {})",
            t2,
            2.0 * t1
        );
        assert!(gate_time >= 0.0, "gate time must be non-negative");

        let temperature = temperature.unwrap_or(0.0);
        let frequency = frequency.unwrap_or(5e9); // 5 GHz typical superconducting qubit frequency

        assert!(temperature >= 0.0, "temperature must be non-negative");
        assert!(frequency > 0.0, "frequency must be positive");

        // Calculate decay probabilities
        let (p1, p2, p_excited) =
            Self::calculate_probabilities(t1, t2, gate_time, temperature, frequency);

        // Construct Kraus operators
        let kraus_ops = if gate_time < EPSILON {
            // No time evolution, only identity
            let dim = 1 << num_qubits;
            let identity = Matrix::eye(dim).mapv(|x| Complex::new(x, 0.0));
            vec![identity]
        } else {
            Self::global_kraus_operators(num_qubits, p1, p2, p_excited)
        };

        Self {
            t1,
            t2,
            gate_time,
            temperature,
            frequency,
            num_qubits,
            kraus_ops,
        }
    }

    /// Create new selective thermal relaxation channel
    /// Applies thermal relaxation only to specified target qubits
    pub fn new_selective(
        total_qubits: usize,
        target_qubits: &[usize],
        t1: f64,
        t2: f64,
        gate_time: f64,
        temperature: Option<f64>,
        frequency: Option<f64>,
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
        assert!(total_qubits > 0, "total number of qubits must be positive");
        assert!(t1 > 0.0, "T₁ time must be positive");
        assert!(t2 > 0.0, "T₂ time must be positive");
        assert!(t2 <= 2.0 * t1 + EPSILON, "T₂ must be ≤ 2T₁");
        assert!(gate_time >= 0.0, "gate time must be non-negative");

        let temperature = temperature.unwrap_or(0.0);
        let frequency = frequency.unwrap_or(5e9);

        assert!(temperature >= 0.0, "temperature must be non-negative");
        assert!(frequency > 0.0, "frequency must be positive");

        // Calculate decay probabilities
        let (p1, p2, p_excited) =
            Self::calculate_probabilities(t1, t2, gate_time, temperature, frequency);

        // Construct selective channel
        let kraus_ops = if gate_time < EPSILON {
            // No time evolution, only identity
            let dim = 1 << total_qubits;
            let identity = Matrix::eye(dim).mapv(|x| Complex::new(x, 0.0));
            vec![identity]
        } else {
            Self::selective_kraus_operators(total_qubits, target_qubits, p1, p2, p_excited)
        };

        Self {
            t1,
            t2,
            gate_time,
            temperature,
            frequency,
            num_qubits: total_qubits,
            kraus_ops,
        }
    }

    /// Create thermal relaxation with explicit probabilities
    /// Direct specification of decay parameters
    pub fn from_probabilities(
        num_qubits: usize,
        p1: f64,        // Amplitude damping probability
        p2: f64,        // Pure dephasing probability
        p_excited: f64, // Thermal excitation probability
    ) -> Self {
        assert!(num_qubits > 0, "number of qubits must be positive");
        assert!((0.0..=1.0).contains(&p1), "p1 must be between 0 and 1");
        assert!((0.0..=1.0).contains(&p2), "p2 must be between 0 and 1");
        assert!(
            (0.0..=1.0).contains(&p_excited),
            "p_excited must be between 0 and 1"
        );

        let kraus_ops = Self::global_kraus_operators(num_qubits, p1, p2, p_excited);

        Self {
            t1: INF,
            t2: INF,
            gate_time: 0.0,
            temperature: 0.0,
            frequency: 5e9,
            num_qubits,
            kraus_ops,
        }
    }

    /// Return the T₁ time
    pub fn t1(&self) -> f64 {
        self.t1
    }

    /// Return the T₂ time
    pub fn t2(&self) -> f64 {
        self.t2
    }

    /// Return the gate time
    pub fn gate_time(&self) -> f64 {
        self.gate_time
    }

    /// Return the temperature
    pub fn temperature(&self) -> f64 {
        self.temperature
    }

    /// Return the frequency
    pub fn frequency(&self) -> f64 {
        self.frequency
    }

    /// Get the effective T₂* time (pure dephasing time)
    /// 1/T₂* = 1/T₂ - 1/(2T₁)
    pub fn t2_star(&self) -> f64 {
        if self.t1.is_infinite() || self.t2.is_infinite() {
            return INF;
        }
        let inv_t2_star = 1.0 / self.t2 - 1.0 / (2.0 * self.t1);
        if inv_t2_star <= 0.0 {
            INF
        } else {
            1.0 / inv_t2_star
        }
    }

    /// Get thermal occupation probability
    pub fn thermal_occupation(&self) -> f64 {
        if self.temperature <= 0.0 {
            return 0.0;
        }

        let hf_kt = (PLANCK * self.frequency) / (BOLTZMANN * self.temperature);
        1.0 / (hf_kt.exp() - 1.0)
    }

    /// Apply thermal relaxation channel to density matrix
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

    /// Calculate decay probabilities from T₁, T₂, gate time, and temperature
    fn calculate_probabilities(
        t1: f64,
        t2: f64,
        gate_time: f64,
        temperature: f64,
        frequency: f64,
    ) -> (f64, f64, f64) {
        // Amplitude damping probability: p₁ = 1 - exp(-t/T₁)
        let p1 = 1.0 - (-gate_time / t1).exp();

        // Pure dephasing probability from T₂*
        // 1/T₂* = 1/T₂ - 1/(2T₁), so T₂* = 1/(1/T₂ - 1/(2T₁))
        let inv_t2_star = 1.0 / t2 - 1.0 / (2.0 * t1);
        let p2 = if inv_t2_star <= 0.0 {
            0.0 // No pure dephasing if T₂ = 2T₁
        } else {
            let t2_star = 1.0 / inv_t2_star;
            1.0 - (-gate_time / t2_star).exp()
        };

        // Thermal excitation probability
        let p_excited = if temperature <= 0.0 {
            0.0
        } else {
            const BOLTZMANN: f64 = 1.380649e-23; // J/K
            const PLANCK: f64 = 6.62607015e-34; // J⋅s

            let hf_kt = (PLANCK * frequency) / (BOLTZMANN * temperature);
            let n_bar = 1.0 / (hf_kt.exp() - 1.0);
            n_bar / (n_bar + 1.0)
        };

        (p1, p2, p_excited)
    }

    /// Create single-qubit thermal relaxation Kraus operators (canonical standard form)
    ///
    /// This is the canonical generalized amplitude damping channel from Nielsen & Chuang
    /// (Eq. 8.79-8.82) which models thermal relaxation at finite temperature.
    ///
    /// The 4 canonical Kraus operators are:
    /// K₀ = (√(1-p)  0    )     K₁ = (0  √((1-p)γ))
    ///      (0       √(1-γ))          (0  0        )
    ///
    /// K₂ = (0  0)              K₃ = (√p      0    )
    ///      (√(pγ) 0)                 (0    √((1-γ)))
    ///
    /// where:
    /// - γ = 1 - exp(-Γt) is the amplitude damping probability
    /// - p = n̄/(n̄+1) is the thermal excitation probability  
    /// - n̄ is the mean thermal photon number
    ///
    /// For dephasing, we add a separate Z-rotation channel.
    fn single_qubit_kraus_ops(p1: f64, p2: f64, p_excited: f64) -> Vec<Matrix<Complex>> {
        let mut ops = Vec::new();

        let gamma = p1; // Amplitude damping probability γ
        let p = p_excited; // Thermal excitation probability

        // The 4 canonical operators for generalized amplitude damping:

        // K₀: Partial preservation
        let mut k0 = Matrix::<Complex>::zeros((2, 2));
        k0[[0, 0]] = Complex::new((1.0 - p).sqrt(), 0.0);
        k0[[1, 1]] = Complex::new((1.0 - gamma).sqrt(), 0.0);
        ops.push(k0);

        // K₁: Amplitude damping |1⟩ → |0⟩
        let mut k1 = Matrix::<Complex>::zeros((2, 2));
        k1[[0, 1]] = Complex::new(((1.0 - p) * gamma).sqrt(), 0.0);
        ops.push(k1);

        // K₂: Thermal excitation |0⟩ → |1⟩
        let mut k2 = Matrix::<Complex>::zeros((2, 2));
        k2[[1, 0]] = Complex::new((p * gamma).sqrt(), 0.0);
        ops.push(k2);

        // K₃: Thermal preservation
        let mut k3 = Matrix::<Complex>::zeros((2, 2));
        k3[[0, 0]] = Complex::new(p.sqrt(), 0.0);
        k3[[1, 1]] = Complex::new((1.0 - gamma).sqrt(), 0.0);
        ops.push(k3);

        // For pure dephasing, we use the standard phase damping channel
        // This is applied as a separate channel in sequence
        if p2 > 1e-12 {
            // Apply dephasing to all existing operators
            let lambda = p2; // Dephasing probability
            let mut dephased_ops = Vec::new();

            // For each existing Kraus operator, create dephased versions
            for op in &ops {
                // No dephasing version (multiply by √(1-λ))
                let mut no_dephase = op.clone();
                no_dephase *= Complex::new((1.0 - lambda).sqrt(), 0.0);
                dephased_ops.push(no_dephase);

                // Dephasing version (apply Z gate, multiply by √λ)
                let mut dephase = op.clone();
                dephase[[1, 1]] *= Complex::new(-1.0, 0.0); // Z gate effect
                dephase *= Complex::new(lambda.sqrt(), 0.0);
                dephased_ops.push(dephase);
            }

            ops = dephased_ops;
        }

        // Remove negligible operators
        ops.retain(|op| op.iter().any(|&x| x.norm() > 1e-12));

        ops
    }

    /// Construct Kraus operators for global thermal relaxation
    fn global_kraus_operators(
        num_qubits: usize,
        p1: f64,
        p2: f64,
        p_excited: f64,
    ) -> Vec<Matrix<Complex>> {
        let single_ops = Self::single_qubit_kraus_ops(p1, p2, p_excited);
        let num_single_ops = single_ops.len();
        let num_ops = num_single_ops.pow(num_qubits as u32);
        let mut kraus_ops = Vec::with_capacity(num_ops);

        // Generate all combinations using base-n representation
        for i in 0..num_ops {
            let mut result = Matrix::from_elem((1, 1), Complex::new(1.0, 0.0));
            let mut temp_i = i;

            // Build tensor product based on base-n representation
            for _ in 0..num_qubits {
                let op_idx = temp_i % num_single_ops;
                temp_i /= num_single_ops;
                result = ops::kron(&result, &single_ops[op_idx]);
            }

            // Only keep non-zero operators
            if result.iter().any(|&x| x.norm() > EPSILON) {
                kraus_ops.push(result);
            }
        }

        kraus_ops
    }

    /// Construct Kraus operators for selective thermal relaxation
    fn selective_kraus_operators(
        total_qubits: usize,
        target_qubits: &[usize],
        p1: f64,
        p2: f64,
        p_excited: f64,
    ) -> Vec<Matrix<Complex>> {
        let single_ops = Self::single_qubit_kraus_ops(p1, p2, p_excited);
        let identity = Identity::new().unitary_matrix();

        let num_target_qubits = target_qubits.len();
        let num_single_ops = single_ops.len();
        let num_ops = num_single_ops.pow(num_target_qubits as u32);
        let mut kraus_ops = Vec::with_capacity(num_ops);

        // Generate all combinations for target qubits
        for i in 0..num_ops {
            let mut result = Matrix::from_elem((1, 1), Complex::new(1.0, 0.0));
            let mut temp_i = i;
            let mut target_idx = 0;

            // Build tensor product for all qubits
            for qubit in 0..total_qubits {
                let op_for_qubit = if target_qubits.contains(&qubit) {
                    let op_idx = temp_i % num_single_ops;
                    if target_idx == 0 {
                        temp_i /= num_single_ops;
                    }
                    target_idx += 1;
                    &single_ops[op_idx]
                } else {
                    &identity
                };

                result = ops::kron(&result, op_for_qubit);
            }

            // Only keep non-zero operators
            if result.iter().any(|&x| x.norm() > EPSILON) {
                kraus_ops.push(result);
            }
        }

        kraus_ops
    }
}
