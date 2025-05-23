use std::fmt;

use crate::types::Complex;

/// Result after long running a quantum circuit
pub struct RunResult {
    pub classical_bit: Vec<Vec<bool>>,
    pub shots: usize,
}

/// Result of a qubit state
pub struct QubitState {
    pub amp_0: Complex,
    pub amp_1: Complex,
    pub prob_0: f64,
    pub prob_1: f64,
}

/// Error types for quantum operations
#[derive(Debug)]
pub enum QuantumError {
    /// Error when trying to operate on incompatible qubits
    IncompatibleQubits,
    /// Error when matrix dimensions are invalid
    InvalidDimension,
    /// Error when qubit index is out of range
    QubitOutOfRange,
}

impl fmt::Display for QuantumError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            QuantumError::IncompatibleQubits => write!(f, "Incompatible qubits for operation"),
            QuantumError::InvalidDimension => write!(f, "Invalid matrix dimension"),
            QuantumError::QubitOutOfRange => write!(f, "Qubit index out of range"),
        }
    }
}
