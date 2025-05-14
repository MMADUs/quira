use crate::math::{Complex, Vector, ONE, ZERO};
use rand::Rng;

pub struct QuantumState {
    /// quantum amplitudes in complex vector
    amplitudes: Vector<Complex>,
    /// number of qubits involved in quantum states
    num_qubits: usize,
}

impl QuantumState {
    /// Create a new quantum state with n qubits in the |0...0⟩ states
    pub fn new(num_qubits: usize) -> Self {
        // left shift (<<) is equal to 2^n as n = number of shift
        // the size of qubit probability in vector space is mapped into 2^k
        // 1 * 2^k as k is the number of qubits
        let size = 1 << num_qubits;
        let mut amplitudes = Vector::zeros(size);
        // 1 + 0i = |0...0⟩ states
        amplitudes[0] = ONE;
        Self {
            amplitudes,
            num_qubits,
        }
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
        // 1e-10 = 1 * 10^-10 = 0.0000000001
        if norm > 1e-10 {
            // replace every element by x/norm
            self.amplitudes.mapv_inplace(|x| x / norm);
        }
    }

    /// Measure all qubits and collapse the state
    pub fn measure_all(&mut self) -> usize {
        let mut rng = rand::rng();
        // gen rand num between 0 and 1 
        let rand_val: f64 = rng.random();
        // calculate cumlative probabilities
        let mut cumlative_prob = 0.0;
        for(idx, amplitude) in self.amplitudes.iter().enumerate() {
            cumlative_prob += (amplitude * amplitude.conj()).re;
            if cumlative_prob > rand_val {
                // collapse state to measured basis state
                let mut new_amplitudes = Vector::zeros(self.amplitudes.len());
                new_amplitudes[idx] = ONE;
                self.amplitudes = new_amplitudes;
                return idx;
            }
        }
        self.amplitudes.len() - 1
    }

    /// Measure a specific qubit and collapse the state accordingly
    pub fn measure_qubit(&mut self, qubit: usize) -> bool {
        if qubit >= self.num_qubits {
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
                self.amplitudes[i] = ZERO;
            }
        }
        // renormalize
        self.normalize();
        outcome
    }
}
