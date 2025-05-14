mod single_qubit;
mod two_qubit;

use crate::math::{Complex, Matrix};
use crate::state::QuantumState;

pub trait QuantumGate: Send + Sync {
    /// Get the unitary matrix representation of the gate
    fn unitary_matrix(&self) -> Matrix<Complex>;

    /// Quantum gate name representation
    fn name(&self) -> String;

    /// Apply the gate to the quantum state
    fn apply(&self, state: &mut QuantumState) {}
}

pub trait SingleQubit: Send + Sync {
    /// Get the target qubit
    fn target_qubit(&self) -> usize;
}

pub trait TwoQubit: Send + Sync {
    /// Get the target qubit
    fn target_qubit(&self) -> usize;

    /// Get the control qubit
    fn control_qubit(&self) -> usize;

    /// Construct the matrix representation
    fn kak_decomposition(&self) {}
}
