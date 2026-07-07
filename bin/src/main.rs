//! Simulation TEST

use quira::{
    include::{Density, QuantumCircuit, QuantumSimulator},
    operation::singleq::Hadamard,
    provider::{QuantumDebugger, QubitIndexing},
};

fn main() {
    // define circuit
    //
    let mut circuit = QuantumCircuit::new(2, 2);
    println!("num qubits: {}", circuit.num_qubits());
    let qreg = circuit.qreg_as_ref();
    circuit.add(Hadamard::new(&qreg[0]));
    let backend = Density::new();
    let mut qsim = QuantumSimulator::new(backend, QubitIndexing::LittleEndian);
    qsim.simulate(circuit);
    let state = qsim.backend().entire_state(false);
    let (rows, cols) = state.dim();
    for row in 0..rows {
        for col in 0..cols {
            print!("{}, ", state[[row, col]]);
        }
        println!()
    }
}

// todo list
// 1. refactor gate
// 2. stabilizer & integration
// 3. refactor noise channel
// 4. quantum fourier transform
// 5. back to the list on obsidian notes.
