pub mod single_qubit;
pub mod two_qubit;

use crate::{
    endian::embed_gate,
    state::QuantumState,
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
    fn apply(&self, state: &mut QuantumState) {
        let n = state.num_qubits();
        let targets = self.construct_targets();
        let unitary = self.unitary_matrix();
        let u = embed_gate(n, &targets, &unitary);
        state.amplitudes = u.dot(&state.amplitudes)
    }

    ///// Multiply this gate with another gate (matrix multiplication)
    //fn multiply<G>(&self, other: &G) -> Result<Box<dyn G>, QuantumError>
    //where
    //    G: QuantumGate,
    //{
    //    if self.target_qubit() != other.target_qubit() {
    //        return Err(QuantumError::IncompatibleQubits);
    //    }
    //
    //    let self_matrix = self.unitary_matrix();
    //    let other_matrix = other.unitary_matrix();
    //
    //    // Perform matrix multiplication
    //    let product_matrix = self_matrix.dot(&other_matrix);
    //
    //    // Return a CompositeGate that represents the product
    //    Ok(Box::new(CompositeGate::new(
    //        self.target_qubit(),
    //        product_matrix,
    //        format!("{}*{}", self.name(), other.name()),
    //    )))
    //}
}

///// A composite gate created by multiplying two or more gates
//pub struct CompositeGate {
//    /// Target qubit index
//    target: Qubit,
//    /// Unitary matrix for this gate
//    matrix: Matrix<Complex>,
//    /// Name of the composite gate
//    name: String,
//}
//
//impl CompositeGate {
//    /// Create a new composite gate
//    pub fn new(target: Qubit, matrix: Matrix<Complex>, name: String) -> Self {
//        Self {
//            target,
//            matrix,
//            name,
//        }
//    }
//}
//
//impl QuantumGate for CompositeGate {
//    fn unitary_matrix(&self) -> Matrix<Complex> {
//        self.matrix.clone()
//    }
//
//    fn name(&self) -> String {
//        self.name.clone()
//    }
//}
