use ndarray::{Array1, Array2};
use ndarray_linalg::*;
use rand::{rng, Rng};

/// Readout error probabilities for a single qubit
#[derive(Debug, Clone)]
pub struct SingleQubitReadoutError {
    /// Probability of measuring 1 when the true state is 0
    pub prob_0_to_1: f64,
    /// Probability of measuring 0 when the true state is 1
    pub prob_1_to_0: f64,
}

impl SingleQubitReadoutError {
    pub fn new(prob_0_to_1: f64, prob_1_to_0: f64) -> Self {
        assert!(prob_0_to_1 >= 0.0 && prob_0_to_1 <= 1.0);
        assert!(prob_1_to_0 >= 0.0 && prob_1_to_0 <= 1.0);
        Self { prob_0_to_1, prob_1_to_0 }
    }
    
    /// Create symmetric readout error (same error rate for both 0->1 and 1->0)
    pub fn symmetric(error_rate: f64) -> Self {
        Self::new(error_rate, error_rate)
    }
}

/// N-qubit readout error model
#[derive(Debug, Clone)]
pub struct ReadoutErrorModel {
    n_qubits: usize,
    /// Per-qubit readout errors (independent model)
    qubit_errors: Vec<SingleQubitReadoutError>,
    /// Optional: Full correlation matrix for correlated readout errors
    correlation_matrix: Option<Array2<f64>>,
}

impl ReadoutErrorModel {
    /// Create independent readout error model
    pub fn independent(qubit_errors: Vec<SingleQubitReadoutError>) -> Self {
        let n_qubits = qubit_errors.len();
        Self {
            n_qubits,
            qubit_errors,
            correlation_matrix: None,
        }
    }
    
    /// Create symmetric independent readout error model
    pub fn symmetric_independent(n_qubits: usize, error_rate: f64) -> Self {
        let qubit_errors = vec![SingleQubitReadoutError::symmetric(error_rate); n_qubits];
        Self::independent(qubit_errors)
    }
    
    /// Create correlated readout error model with full 2^n x 2^n confusion matrix
    pub fn correlated(n_qubits: usize, confusion_matrix: Array2<f64>) -> Self {
        let expected_size = 1 << n_qubits; // 2^n
        assert_eq!(confusion_matrix.shape(), [expected_size, expected_size]);
        
        // Verify it's a valid stochastic matrix (rows sum to 1)
        for row in confusion_matrix.rows() {
            let sum = row.sum();
            assert!((sum - 1.0).abs() < 1e-10, "Confusion matrix rows must sum to 1");
        }
        
        Self {
            n_qubits,
            qubit_errors: vec![], // Not used in correlated mode
            correlation_matrix: Some(confusion_matrix),
        }
    }
    
    /// Apply readout error to measurement outcomes
    pub fn apply_readout_error(&self, ideal_measurements: &[String]) -> Vec<String> {
        ideal_measurements.iter()
            .map(|measurement| self.apply_single_measurement(measurement))
            .collect()
    }
    
    /// Apply readout error to a single measurement string (e.g., "0101")
    pub fn apply_single_measurement(&self, ideal_measurement: &str) -> String {
        assert_eq!(ideal_measurement.len(), self.n_qubits);
        
        if let Some(ref confusion_matrix) = self.correlation_matrix {
            self.apply_correlated_error(ideal_measurement, confusion_matrix)
        } else {
            self.apply_independent_error(ideal_measurement)
        }
    }
    
    /// Apply independent readout errors to each qubit
    fn apply_independent_error(&self, ideal_measurement: &str) -> String {
        let mut rng = rng();
        let mut noisy_bits = String::new();
        
        for (i, bit_char) in ideal_measurement.chars().enumerate() {
            let ideal_bit = bit_char.to_digit(2).unwrap() as usize;
            let error_model = &self.qubit_errors[i];
            
            let flip_prob = match ideal_bit {
                0 => error_model.prob_0_to_1,
                1 => error_model.prob_1_to_0,
                _ => panic!("Invalid bit value: {}", ideal_bit),
            };
            
            let noisy_bit = if rng.random::<f64>() < flip_prob {
                1 - ideal_bit // Flip the bit
            } else {
                ideal_bit // Keep the bit
            };
            
            noisy_bits.push(char::from_digit(noisy_bit as u32, 2).unwrap());
        }
        
        noisy_bits
    }
    
    /// Apply correlated readout errors using confusion matrix
    fn apply_correlated_error(&self, ideal_measurement: &str, confusion_matrix: &Array2<f64>) -> String {
        let ideal_state = binary_string_to_index(ideal_measurement);
        let mut rng = rng();
        
        // Sample from the confusion matrix row corresponding to ideal state
        let probabilities = confusion_matrix.row(ideal_state);
        let noisy_state = sample_from_distribution(&probabilities.to_owned(), &mut rng);
        
        index_to_binary_string(noisy_state, self.n_qubits)
    }
    
    /// Get the confusion matrix (either from stored matrix or computed from independent errors)
    pub fn get_confusion_matrix(&self) -> Array2<f64> {
        if let Some(ref matrix) = self.correlation_matrix {
            matrix.clone()
        } else {
            self.compute_independent_confusion_matrix()
        }
    }
    
    /// Compute confusion matrix for independent readout errors
    fn compute_independent_confusion_matrix(&self) -> Array2<f64> {
        let n_states = 1 << self.n_qubits; // 2^n
        let mut confusion_matrix = Array2::<f64>::zeros((n_states, n_states));
        
        for ideal_state in 0..n_states {
            let ideal_bits = index_to_binary_string(ideal_state, self.n_qubits);
            
            // For each possible noisy outcome
            for noisy_state in 0..n_states {
                let noisy_bits = index_to_binary_string(noisy_state, self.n_qubits);
                
                // Calculate probability of this outcome given ideal state
                let mut prob = 1.0;
                for (i, (ideal_char, noisy_char)) in ideal_bits.chars().zip(noisy_bits.chars()).enumerate() {
                    let ideal_bit = ideal_char.to_digit(2).unwrap() as usize;
                    let noisy_bit = noisy_char.to_digit(2).unwrap() as usize;
                    let error_model = &self.qubit_errors[i];
                    
                    let bit_prob = match (ideal_bit, noisy_bit) {
                        (0, 0) => 1.0 - error_model.prob_0_to_1,
                        (0, 1) => error_model.prob_0_to_1,
                        (1, 0) => error_model.prob_1_to_0,
                        (1, 1) => 1.0 - error_model.prob_1_to_0,
                        _ => unreachable!(),
                    };
                    
                    prob *= bit_prob;
                }
                
                confusion_matrix[[ideal_state, noisy_state]] = prob;
            }
        }
        
        confusion_matrix
    }
    
    /// Apply readout error to quantum state measurement probabilities
    pub fn apply_to_probabilities(&self, ideal_probs: &Array1<f64>) -> Array1<f64> {
        let confusion_matrix = self.get_confusion_matrix();
        
        // Noisy probabilities = confusion_matrix^T * ideal_probs
        // (We transpose because confusion_matrix[i,j] = P(measure j | ideal i))
        confusion_matrix.t().dot(ideal_probs)
    }
    
    /// Create mitigation matrix (inverse of confusion matrix) for error correction
    pub fn create_mitigation_matrix(&self) -> Result<Array2<f64>, Box<dyn std::error::Error>> {
        let confusion_matrix = self.get_confusion_matrix();
        
        // For readout error mitigation, we need to invert the confusion matrix
        let mitigation_matrix = confusion_matrix.inv()?;
        
        Ok(mitigation_matrix)
    }
    
    /// Apply readout error mitigation to measured probabilities
    pub fn mitigate_readout_error(&self, noisy_probs: &Array1<f64>) -> Result<Array1<f64>, Box<dyn std::error::Error>> {
        let mitigation_matrix = self.create_mitigation_matrix()?;
        let mitigated_probs = mitigation_matrix.dot(noisy_probs);
        
        // Ensure probabilities are non-negative and normalized
        let mitigated_probs = mitigated_probs.mapv(|x| x.max(0.0));
        let sum = mitigated_probs.sum();
        let normalized_probs = mitigated_probs / sum;
        
        Ok(normalized_probs)
    }
}

/// Convert binary string to state index
fn binary_string_to_index(binary_str: &str) -> usize {
    usize::from_str_radix(binary_str, 2).unwrap()
}

/// Convert state index to binary string
fn index_to_binary_string(index: usize, n_bits: usize) -> String {
    format!("{:0width$b}", index, width = n_bits)
}

/// Sample from a discrete probability distribution
fn sample_from_distribution(probabilities: &Array1<f64>, rng: &mut impl Rng) -> usize {
    let random_value = rng.random::<f64>();
    let mut cumulative_prob = 0.0;
    
    for (i, &prob) in probabilities.iter().enumerate() {
        cumulative_prob += prob;
        if random_value <= cumulative_prob {
            return i;
        }
    }
    
    // Should not reach here if probabilities sum to 1
    probabilities.len() - 1
}
