pub mod single_qubit;
pub mod two_qubit;

use crate::{
    endian::{expand_unitary, QubitIndexing},
    kernel::QuantumState,
    types::{Complex, Matrix, Qubit},
};

pub trait QuantumGate: Send + Sync {
    /// Get the unitary matrix representation of the gate
    fn unitary_matrix(&self) -> Matrix<Complex>;

    /// Quantum gate name representation
    fn name(&self) -> String;

    /// Get the target qubits this gate operates on
    fn construct_targets(&self) -> Vec<Qubit>;
}

// Helper trait (optional)
pub trait GateApplyExt {
    fn apply<T: QuantumState>(&self, state: &mut T, indexing: &QubitIndexing);
}

// Implement for all references to types that implement QuantumGate
impl<T: QuantumGate + ?Sized> GateApplyExt for T {
    fn apply<U: QuantumState>(&self, state: &mut U, indexing: &QubitIndexing) {
        let n = state.num_qubits();
        let targets = self.construct_targets();
        let unitary = self.unitary_matrix();
        let u = expand_unitary(n, &targets, &unitary, indexing);
        state.apply(u);
    }
}

