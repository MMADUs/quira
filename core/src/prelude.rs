use crate::{QuantumGate, Qubit, Vector, types::Complex};

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

/// A variant of qubit insertion
pub enum QubitToken {
    FROM(Vector<Complex>),
    ONES(usize),
    ZEROS(usize),
    RESET,
}

/// Quantum operation tokens
pub enum QuantumTokens {
    /// applying qubit to state
    Qubit(QubitToken),
    /// apply gate operation
    Operations(Box<dyn QuantumGate>),
    /// conditionally apply gate operation
    /// 'classical_bit', 'measured', 'operations'.
    Conditional((usize, bool, Box<dyn QuantumGate>)),
    /// qubit measurements
    /// 'qubit_index', 'classical_bit'.
    Measurements((Qubit, usize)),
    /// circuit barrier
    Barrier,
}
