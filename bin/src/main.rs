//! Simulation TEST

use quira::include::{
    ClassicalRegister, QuantumCircuit, QuantumRegister, QuantumSimulator, StateVec,
};
use quira::operation::singleq::Hadamard as H;
use quira::provider::{QuantumDebugger, QubitIndexing as QI};

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
    circuit.add(H::new(&qreg[0]));
    circuit.measure_all();
    // define backend 
    //
    let backend = StateVec::new();
    // simulate circuit 
    //
    let mut qs = QuantumSimulator::new(backend, QI::LittleEndian);
    qs.simulate(circuit);
    // see measured outcome in classical register 
    //
    let creg_result = qs.creg();
    println!("classical bit: {}", creg_result.reg(&creg[0]));
    // inspect full state 
    //
    let backend_state = qs.backend();
    for (bit, amplitude) in backend_state.entire_state(false) {
        println!("|{}>: {}", bit, amplitude);
    }
}

// todo list
// 1. refactor gate
// 2. stabilizer & integration
// 3. refactor noise channel
// 4. quantum fourier transform
// 5. back to the list on obsidian notes.
