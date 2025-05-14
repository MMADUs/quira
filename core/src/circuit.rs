use crate::operations::QuantumGate;

pub struct Circuit {
    num_qubits: usize,
    operations: Vec<Box<dyn QuantumGate>>,
}
