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

use crate::constant::{C_ONE, C_ZERO, EPSILON};
use crate::prelude::QubitState;
use crate::types::{Complex, Qubit, Vector};
use rand::Rng;

#[derive(Clone)]
/// the quantum state is the representation for the entire state of the qubit
/// the state is represented as a hilbert space.
pub struct QuantumStateVec {
    /// quantum amplitudes in complex vector
    amplitudes: Vector<Complex>,
}

impl QuantumStateVec {
    /// Initialize empty state
    pub fn new() -> Self {
        Self {
            amplitudes: Vector::from(Vec::new()),
        }
    }

    /// Check if quantum state is empty
    pub fn is_empty(&self) -> bool {
        self.amplitudes.len() == 0
    }

    /// Set existing state with new state
    pub fn set(&mut self, amplitudes: Vector<Complex>) {
        assert!(
            amplitudes.len().is_power_of_two(),
            "Amplitude of the vector must have length 2^n"
        );
        self.amplitudes = amplitudes;
        self.normalize();
    }

    /// Create new qubit from a set of vector amplitudes.
    pub fn from(amplitudes: Vector<Complex>) -> Self {
        let len = amplitudes.len();
        assert!(
            len.is_power_of_two(),
            "Amplitude of the vector must have length 2^n"
        );
        let mut qs = Self { amplitudes };
        qs.normalize();
        qs
    }
    
    /// Expanding the current state vector by kronecker product
    pub fn expand_state(&mut self, state: Vector<Complex>) {
        if self.is_empty() {
            self.amplitudes = Self::from(state).amplitudes;
            return;
        }
        let mut new_amplitudes = Vec::with_capacity(self.amplitudes.len() * state.len());
        for &a in self.amplitudes.iter() {
            for &b in state.iter() {
                new_amplitudes.push(a * b)
            }
        }
        self.amplitudes = Vector::from(new_amplitudes);
        self.normalize();
    }

    /// Returns the number of qubits (log2 of dimension).
    pub fn num_qubits(&self) -> usize {
        if self.is_empty() {
            0
        } else {
            (self.amplitudes.len() as f64).log2() as usize
        }
    }

    /// Get the reference of the entire amplitudes
    pub fn amplitudes_as_ref(&self) -> &Vector<Complex> {
        &self.amplitudes
    }

    /// Get the probability amplitude for a specific basis state
    pub fn amplitude(&self, basis_state: usize) -> Complex {
        self.amplitudes[basis_state]
    }

    /// Calculate the probability of measuring a specific basis state
    pub fn probability(&self, basis_state: usize) -> f64 {
        let amplitude = self.amplitudes[basis_state];
        // P(i) = |α[i]|^2 = α[i] * α[i]^*
        (amplitude * amplitude.conj()).re
    }

    /// Calculate the state vector norm (length)
    pub fn norm(&self) -> f64 {
        // ∥∣ψ⟩∥ = sqrt(sum(P(i)))
        self.amplitudes
            .iter()
            .map(|amp| (amp * amp.conj()).re)
            .sum::<f64>()
            .sqrt()
    }

    /// Normalize the state vector
    pub fn normalize(&mut self) {
        let norm = self.norm();
        // check near-zero to prevent division by zero or a very little number
        if norm > EPSILON {
            // replace every element by x/norm
            self.amplitudes.mapv_inplace(|x| x / norm);
        }
    }

    /// Return the inner product with another qubit.
    pub fn inner_product(&self, other: &QuantumStateVec) -> Complex {
        self.amplitudes.dot(&other.amplitudes)
    }

    /// Measure all qubits and collapse the state
    pub fn measure_all(&mut self) -> usize {
        let mut rng = rand::rng();
        // gen rand num between 0 and 1
        let rand_val: f64 = rng.random();
        // calculate cumlative probabilities
        let mut cumlative_prob = 0.0;
        for (idx, amplitude) in self.amplitudes.iter().enumerate() {
            cumlative_prob += (amplitude * amplitude.conj()).re;
            if cumlative_prob > rand_val {
                // collapse state to measured basis state
                let mut new_amplitudes = Vector::zeros(self.amplitudes.len());
                new_amplitudes[idx] = C_ONE;
                self.amplitudes = new_amplitudes;
                return idx;
            }
        }
        self.amplitudes.len() - 1
    }

    /// Measure a specific qubit and collapse the state accordingly
    pub fn measure_qubit(&mut self, qubit: usize) -> bool {
        if qubit >= self.num_qubits() {
            panic!("Qubit index out of range");
        }
        let mut rng = rand::rng();
        // calculate probability measuring |1⟩
        let mut prob_one = 0.0;
        for i in 0..self.amplitudes.len() {
            if (i & (1 << qubit)) != 0 {
                prob_one += (self.amplitudes[i] * self.amplitudes[i].conj()).re;
            }
        }
        // choose outcome based on probability
        let outcome = rng.random::<f64>() < prob_one;
        // collapse state based on measurement
        for i in 0..self.amplitudes.len() {
            let bit_is_set = (i & (1 << qubit)) != 0;
            if bit_is_set != outcome {
                self.amplitudes[i] = C_ZERO;
            }
        }
        // renormalize
        self.normalize();
        outcome
    }

    /// Return the quantum state as formatted strings.
    pub fn get_amp_state(&self, filter_zero: bool) -> Vec<String> {
        let num_qubits = self.num_qubits();
        self.amplitudes
            .iter()
            .enumerate()
            .filter(|(_, amp)| {
                if filter_zero {
                    amp.norm_sqr() > 1e-12
                } else {
                    true
                }
            })
            .map(|(i, amp)| {
                let prob = amp.norm_sqr();
                let bin = format!("{:0width$b}", i, width = num_qubits);
                format!(
                    "[{}] |{}⟩: (amp = {:.4} + {:.4}i) => (prob = {:.4}, {:.2}%)",
                    i,
                    bin,
                    amp.re,
                    amp.im,
                    prob,
                    prob * 100.0
                )
            })
            .collect()
    }

    /// Get specific state of a qubit
    pub fn get_qubit_state(&self, qubit: Qubit) -> QubitState {
        let num_qubits = self.num_qubits();

        if qubit >= num_qubits {
            panic!(
                "\nError: Qubit index {} out of range. System has {} qubits.",
                qubit, num_qubits
            );
        }

        println!("\nQubit {} state:", qubit);

        // Calculate marginal amplitudes for |0⟩ and |1⟩
        let mut amp_0 = Complex::new(0.0, 0.0);
        let mut amp_1 = Complex::new(0.0, 0.0);

        for (i, amp) in self.amplitudes.iter().enumerate() {
            let qubit_bit = (i >> qubit) & 1;
            if qubit_bit == 0 {
                amp_0 += amp;
            } else {
                amp_1 += amp;
            }
        }

        let prob_0 = amp_0.norm_sqr();
        let prob_1 = amp_1.norm_sqr();

        println!(
            "|0⟩: (amp = {:.4} + {:.4}i) => (prob = {:.4}, {:.2}%)",
            amp_0.re,
            amp_0.im,
            prob_0,
            prob_0 * 100.0
        );
        println!(
            "|1⟩: (amp = {:.4} + {:.4}i) => (prob = {:.4}, {:.2}%)",
            amp_1.re,
            amp_1.im,
            prob_1,
            prob_1 * 100.0
        );

        QubitState {
            amp_0,
            amp_1,
            prob_0,
            prob_1,
        }
    }
}
