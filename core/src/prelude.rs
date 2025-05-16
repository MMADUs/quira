use std::fmt;

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
