pub mod single;

use crate::prelude::QuantumError;
use crate::state::QuantumState;
use crate::types::{Complex, Matrix, Qubit};

use rayon::prelude::*;

pub trait QuantumGate: Send + Sync {
    /// Get the unitary matrix representation of the gate
    fn unitary_matrix(&self) -> Matrix<Complex>;

    /// Quantum gate name representation
    fn name(&self) -> String;

    /// Get the target qubit for this gate
    fn target_qubit(&self) -> Qubit;

    /// Apply the gate to a quantum state
    fn apply(&self, state: &mut QuantumState) {
        let n = state.num_qubits();
        let target = self.target_qubit();
        let matrix = self.unitary_matrix();

        // Apply the gate in parallel using chunks
        // Parallelization is effective for large state vectors
        if state.amplitudes().len() > 1024 {
            // Calculate chunk size for parallel processing
            let chunk_size = 2 << target;

            // Create a mutable slice for parallel processing
            let amplitudes = state.amplitudes_mut();

            // Process chunks in parallel
            amplitudes.par_chunks_mut(chunk_size).for_each(|chunk| {
                // For each chunk, process all pairs of amplitudes
                for i in 0..(chunk_size / 2) {
                    let i0 = i;
                    let i1 = i + (chunk_size / 2);

                    if i1 < chunk.len() {
                        let a0 = chunk[i0];
                        let a1 = chunk[i1];

                        // Apply matrix
                        chunk[i0] = matrix[[0, 0]] * a0 + matrix[[0, 1]] * a1;
                        chunk[i1] = matrix[[1, 0]] * a0 + matrix[[1, 1]] * a1;
                    }
                }
            });
        } else {
            // For smaller state vectors, process sequentially
            let size = 1 << n;
            let mask = 1 << target;

            // Create a copy of the current amplitudes
            let old_amplitudes = state.amplitudes().to_vec();

            for i in 0..size {
                // Skip if we've already processed this pair
                if i & mask != 0 {
                    continue;
                }

                let i0 = i; // Target qubit is 0
                let i1 = i | mask; // Target qubit is 1

                let a0 = old_amplitudes[i0];
                let a1 = old_amplitudes[i1];

                // Apply matrix - get a mutable reference to amplitudes
                let amplitudes = state.amplitudes_mut();
                amplitudes[i0] = matrix[[0, 0]] * a0 + matrix[[0, 1]] * a1;
                amplitudes[i1] = matrix[[1, 0]] * a0 + matrix[[1, 1]] * a1;
            }
        }
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

/// A composite gate created by multiplying two or more gates
pub struct CompositeGate {
    /// Target qubit index
    target: Qubit,
    /// Unitary matrix for this gate
    matrix: Matrix<Complex>,
    /// Name of the composite gate
    name: String,
}

impl CompositeGate {
    /// Create a new composite gate
    pub fn new(target: Qubit, matrix: Matrix<Complex>, name: String) -> Self {
        Self {
            target,
            matrix,
            name,
        }
    }
}

impl QuantumGate for CompositeGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        self.matrix.clone()
    }

    fn name(&self) -> String {
        self.name.clone()
    }

    fn target_qubit(&self) -> Qubit {
        self.target
    }
}
