//! Quantum Teleportation example
//! teleportation.rs
//!
//! run example with: cargo run --example teleportation

use quira::QuantumCircuit as QC;
use quira::Qubit;
use quira::QubitIndexing as QI;
use quira::RunResult as RR;

use quira::SingleQ::{Hadamard as H, PauliX as X, PauliZ as Z};
use quira::TwoQ::ControlledNot as CNOT;

fn main() {
    // make a new quantum circuit
    let mut qc: QC = QC::new(true, QI::LittleEndian); // -> Qiskit uses Little Endian.

    // prepare qubit state to teleport
    let phi: Qubit = qc.qb_zero();

    // prepare Alice's and Bob's qubit
    let q0: Qubit = qc.qb_zero(); // -> Alice's qubit
    let q1: Qubit = qc.qb_zero(); // -> Bob's qubit

    // apply hadamard to phi
    qc.add(H::new(phi)); // -> superposition state we want to teleport

    qc.state_vector(false);
    // [0] |000⟩: (amp = 0.7071 + 0.0000i) => (prob = 0.5000, 50.00%)
    // [1] |001⟩: (amp = 0.7071 + 0.0000i) => (prob = 0.5000, 50.00%)
    // [2] |010⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [3] |011⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [4] |100⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [5] |101⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [6] |110⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [7] |111⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)

    // apply hadamard & controlled not
    qc.add(H::new(q0)); // -> superposition in Alice's qubit
    qc.add(CNOT::new(q0, q1)); // entangled Alice's and Bob's qubit

    qc.state_vector(false);
    // [0] |000⟩: (amp = 0.5000 + 0.0000i) => (prob = 0.2500, 25.00%)
    // [1] |001⟩: (amp = 0.5000 + 0.0000i) => (prob = 0.2500, 25.00%)
    // [2] |010⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [3] |011⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [4] |100⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [5] |101⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [6] |110⟩: (amp = 0.5000 + 0.0000i) => (prob = 0.2500, 25.00%)
    // [7] |111⟩: (amp = 0.5000 + 0.0000i) => (prob = 0.2500, 25.00%)

    // apply controlled not & hadamard
    qc.add(CNOT::new(phi, q0)); // -> entagled phi and Alice's qubit
    qc.add(H::new(phi)); // -> finishes off phi bell state

    qc.state_vector(true);
    // [0] |000⟩: (amp = 0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)
    // [1] |001⟩: (amp = 0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)
    // [2] |010⟩: (amp = 0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)
    // [3] |011⟩: (amp = -0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)
    // [4] |100⟩: (amp = 0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)
    // [5] |101⟩: (amp = -0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)
    // [6] |110⟩: (amp = 0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)
    // [7] |111⟩: (amp = 0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)

    let m0 = qc.measure(phi, 0); // measure phi to classical reg 0
    let m1 = qc.measure(q0, 1); // measure Alice's qubit to classical reg 1

    qc.state_vector(true);
    // [2] |010⟩: (amp = 0.7071 + 0.0000i) => (prob = 0.5000, 50.00%)
    // [6] |110⟩: (amp = 0.7071 + 0.0000i) => (prob = 0.5000, 50.00%)

    // apply conditional measurement
    // if any of the measurement is |1>, it applies the gate.
    qc.cond_add(m0, X::new(q1)); // -> if Alice measured |1>, Bob applies Pauli X
    qc.cond_add(m1, Z::new(q1)); // -> if Bob measured |1>, Alice applies Pauli Z

    println!("m0 = {}, m1 = {}", m0, m1); // -> superposition outcome

    qc.state_vector(true); // -> state is either Bob or Alice

    qc.state_qubit(q1); // -> verify, now phi state should be teleported to Bob's qubit state

    let shots: usize = 10000;
    let result: RR = qc.run(shots); // -> run cirucit for N shots

    qc.analyze(result); // -> analyze the statistic outcome from N shots result
}
