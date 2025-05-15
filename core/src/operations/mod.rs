mod single;

use crate::math::{Complex, Matrix};
use crate::state::QuantumState;

pub trait QuantumGate: Send + Sync {
    /// Get the unitary matrix representation of the gate
    fn unitary_matrix(&self) -> Matrix<Complex>;

    /// Quantum gate name representation
    fn name(&self) -> String;

    /// Apply the gate to a quantum state
    fn apply(&self, state: &mut QuantumState) {
        let n = state.num_qubits;
        let target = self.target_qubit();
        let matrix = self.unitary_matrix();
        
        // Apply the gate in parallel using chunks
        // Parallelization is effective for large state vectors
        if state.amplitudes.len() > 1024 {
            // Get a mutable reference to the amplitudes data
            let amp_slice = state.amplitudes.as_slice_mut().unwrap();
            
            // Calculate chunk size for parallel processing
            let chunk_size = 2 << target;
            
            // Process chunks in parallel
            amp_slice.par_chunks_mut(chunk_size)
                .for_each(|chunk| {
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
            let old_amplitudes = state.amplitudes.clone();
            
            for i in 0..size {
                // Skip if we've already processed this pair
                if i & mask != 0 {
                    continue;
                }
                
                let i0 = i;         // Target qubit is 0
                let i1 = i | mask;  // Target qubit is 1
                
                let a0 = old_amplitudes[i0];
                let a1 = old_amplitudes[i1];
                
                // Apply matrix
                state.amplitudes[i0] = matrix[[0, 0]] * a0 + matrix[[0, 1]] * a1;
                state.amplitudes[i1] = matrix[[1, 0]] * a0 + matrix[[1, 1]] * a1;
            }
        }
    }
    
    /// Multiply this gate with another gate (matrix multiplication)
    fn multiply<G: Gate>(&self, other: &G) -> Result<Box<dyn Gate>, QuantumError> {
        if self.target_qubit() != other.target_qubit() {
            return Err(QuantumError::IncompatibleQubits);
        }
        
        let self_matrix = self.unitary_matrix();
        let other_matrix = other.unitary_matrix();
        
        // Perform matrix multiplication
        let product_matrix = self_matrix.dot(&other_matrix);
        
        // Return a CompositeGate that represents the product
        Ok(Box::new(CompositeGate::new(
            self.target_qubit(),
            product_matrix,
            format!("{}*{}", self.name(), other.name())
        )))
    }
}



pub trait TwoQubit: Send + Sync {
    /// Get the target qubit
    fn target_qubit(&self) -> usize;

    /// Get the control qubit
    fn control_qubit(&self) -> usize;

    /// Construct the matrix representation
    fn kak_decomposition(&self) {}
}
