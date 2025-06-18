use crate::{Complex, Qubit, Vector};

pub mod density;
pub mod statevec;
pub mod tensornet;

/// every kernel (quantum state) implements backend operation
/// for executing circuit in the backend
pub trait BackendOperation {
    /// expand the kernel state
    fn expand_state(&mut self, state: Vector<Complex>);

    /// measure qubit on the kernel
    fn measure_qubit(&mut self, qubit: Qubit) -> bool;

    /// measure all qubit on the kernel
    fn measure_all(&mut self) -> Vec<bool>;
}

/// every kernel (quantum state) implements quantum state
/// for applying quantum gate to state
pub trait QuantumState {
    /// the state types
    type StateType;

    /// reference of the kernel state
    fn state_as_ref(&self) -> &Self::StateType;
}
