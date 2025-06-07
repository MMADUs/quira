pub mod single_qubit;
pub mod two_qubit;

use crate::{
    endian::{QubitIndexing, expand_unitary},
    statevec::QuantumStateVec,
    types::{Complex, Matrix, Qubit},
};

pub trait QuantumGate: Send + Sync {
    /// Get the unitary matrix representation of the gate
    fn unitary_matrix(&self) -> Matrix<Complex>;

    /// Quantum gate name representation
    fn name(&self) -> String;

    /// Get the target qubits this gate operates on
    fn construct_targets(&self) -> Vec<Qubit>;

    /// Apply this gate to the given quantum state
    fn apply(&self, state: &mut QuantumStateVec, indexing: &QubitIndexing) {
        let n = state.num_qubits();
        let targets = self.construct_targets();
        let unitary = self.unitary_matrix();
        // apply gate to state
        let u = expand_unitary(n, &targets, &unitary, &indexing);
        let current_state = state.amplitudes_as_ref();
        let new_state = u.dot(current_state);
        state.set(new_state)
    }
}
