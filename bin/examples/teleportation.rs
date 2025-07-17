//! Quantum Teleportation example
//! teleportation.rs
//!
//! run example with: cargo run --example teleportation

use quira::{
    common::{Complex, PI},
    include::{ClassicalRegister, QuantumCircuit, QuantumJob, QuantumRegister, QuantumSimulator, StateVec},
    operation::{
        singleq::{Hadamard as H, PauliX as X, PauliZ as Z},
        twoq::ControlledNot as CNOT,
    },
    provider::{QuantumDebugger, QubitIndexing as QI},
};

fn main() {
    // state we want to teleport
    //
    let theta = PI / 8.0;
    let alpha = Complex::new(theta.cos(), 0.0);
    let beta = Complex::new(0.0, theta.sin());
    println!("state to teleport");
    println!("{}|0> + {}|1> = 1", alpha, beta);
    // define register
    //
    let qreg = QuantumRegister::new(2);
    let creg = ClassicalRegister::new(2);
    // define circuit
    //
    let mut circuit = QuantumCircuit::with_reg(&qreg, &creg);
    let phi = circuit.qb_from_iter([alpha, beta]); // -> we teleport this qubit state
    // build circuit
    //
    circuit.add(H::new(&qreg[0]));
    circuit.add(CNOT::new(&qreg[0], &qreg[1]));
    circuit.add(CNOT::new(&phi, &qreg[0]));
    circuit.add(H::new(&phi));
    // define backend
    //
    let backend = StateVec::new();
    // simulate circuit
    //
    let mut qsim = QuantumSimulator::new(backend, QI::LittleEndian);
    qsim.simulate(circuit);
    // inspect state
    //
    let backend_state = qsim.backend();
    let (tp_alpha, tp_beta) = backend_state.qubit_state(&qreg[1]);
    println!("teleported state");
    println!("{}|0> + {}|1> = 1", tp_alpha, tp_beta);
    // quantum hardware execution
    //
    let mut circuit = QuantumCircuit::with_reg(&qreg, &creg);
    let phi = circuit.qb_from_iter([alpha, beta]);
    circuit.add(H::new(&qreg[0]));
    circuit.add(CNOT::new(&qreg[0], &qreg[1]));
    circuit.add(CNOT::new(&phi, &qreg[0]));
    circuit.add(H::new(&phi));
    circuit.measure(&phi, &creg[0]);
    circuit.measure(&qreg[0], &creg[1]);
    circuit.cond_add(&creg[0], 1, X::new(&qreg[1]));
    circuit.cond_add(&creg[1], 1, Z::new(&qreg[1]));
    // define backend 
    //
    let backend = StateVec::new();
    // execute job
    //
    let mut job = QuantumJob::new(backend, QI::LittleEndian);
    job.from_circuit(circuit);
    job.run(100000);
    let results = job.result();
    for result in results.unwrap() {
        // analyze each result
        // 
        for (bit, count) in result.get_counts() {
            println!("{}: {}", bit, count);
        }
        let (mb, mc) = result.get_most_frequent().unwrap();
        println!("|{}>: {}", mb, mc);
        result.print_results();
        println!("{}", result.histogram(20));
    }
}
