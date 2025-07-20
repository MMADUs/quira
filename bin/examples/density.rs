//! Density Matrix Superposition example
//! density.rs
//!
//! run example with: cargo run --example density

use quira::{
    include::{ClassicalRegister, Density, QuantumCircuit, QuantumRegister, QuantumSimulator},
    operation::singleq::Hadamard,
    provider::{QuantumDebugger, QubitIndexing},
};

fn main() {
    // define registers
    //
    let qreg = QuantumRegister::new(1).label("qreg 1");
    let creg = ClassicalRegister::new(1).label("creg 1");
    // define circuit
    //
    let mut circuit = QuantumCircuit::with_reg(&qreg, &creg);
    // build circuit
    //
    circuit.add(Hadamard::new(&qreg[0]));
    circuit.measure_all();
    // define backend
    //
    let backend = Density::new();
    // simulate circuit
    //
    let mut qs = QuantumSimulator::new(backend, QubitIndexing::LittleEndian);
    qs.simulate(circuit);
    // see measured outcome in classical register
    //
    let creg_result = qs.creg();
    println!("classical bit: {}", creg_result.reg(&creg[0]));
    // inspect full state
    //
    let state = qs.backend().entire_state(false);
    let (rows, cols) = state.dim();
    for row in 0..rows {
        for col in 0..cols {
            print!("{}, ", state[[row, col]]);
        }
        println!()
    }
}
