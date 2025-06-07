//! Quantum Teleportation example
//! teleportation.rs
//!
//! run example with: cargo run --example teleportation

use quira::QuantumSimulator as QS;
use quira::Qubit;
use quira::QubitIndexing as QI;

use quira::SingleQ::{Hadamard as H, PauliX as X, PauliZ as Z};
use quira::TwoQ::ControlledNot as CNOT;

fn main() {
    // make a new quantum circuit
    let mut qsim: QS = QS::new(QI::LittleEndian); // -> Qiskit uses Little Endian.

    // prepare qubit state to teleport
    let phi: Qubit = qsim.qb_zero();

    // prepare Alice's and Bob's qubit
    let q0: Qubit = qsim.qb_zero(); // -> Alice's qubit
    let q1: Qubit = qsim.qb_zero(); // -> Bob's qubit

    // apply hadamard to phi
    qsim.add(H::new(phi)); // -> superposition state we want to teleport

    qsim.state_vec(false);
    // [0] |000⟩: (amp = 0.7071 + 0.0000i) => (prob = 0.5000, 50.00%)
    // [1] |001⟩: (amp = 0.7071 + 0.0000i) => (prob = 0.5000, 50.00%)
    // [2] |010⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [3] |011⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [4] |100⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [5] |101⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [6] |110⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [7] |111⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)

    // apply hadamard & controlled not
    qsim.add(H::new(q0)); // -> superposition in Alice's qubit
    qsim.add(CNOT::new(q0, q1)); // entangled Alice's and Bob's qubit

    qsim.state_vec(false);
    // [0] |000⟩: (amp = 0.5000 + 0.0000i) => (prob = 0.2500, 25.00%)
    // [1] |001⟩: (amp = 0.5000 + 0.0000i) => (prob = 0.2500, 25.00%)
    // [2] |010⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [3] |011⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [4] |100⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [5] |101⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [6] |110⟩: (amp = 0.5000 + 0.0000i) => (prob = 0.2500, 25.00%)
    // [7] |111⟩: (amp = 0.5000 + 0.0000i) => (prob = 0.2500, 25.00%)

    // apply controlled not & hadamard
    qsim.add(CNOT::new(phi, q0)); // -> entagled phi and Alice's qubit
    qsim.add(H::new(phi)); // -> finishes off phi bell state

    qsim.state_vec(true);
    // [0] |000⟩: (amp = 0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)
    // [1] |001⟩: (amp = 0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)
    // [2] |010⟩: (amp = 0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)
    // [3] |011⟩: (amp = -0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)
    // [4] |100⟩: (amp = 0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)
    // [5] |101⟩: (amp = -0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)
    // [6] |110⟩: (amp = 0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)
    // [7] |111⟩: (amp = 0.3536 + 0.0000i) => (prob = 0.1250, 12.50%)

    let m0: bool = qsim.measure(phi, 0); // measure phi to classical reg 0
    let m1: bool = qsim.measure(q0, 1); // measure Alice's qubit to classical reg 1

    println!("m0 = {}, m1 = {}", m0 as u8, m1 as u8); // -> superposition outcome

    qsim.state_vec(true); // -> superposition state before teleportaion

    // apply conditional measurement for teleportation
    // if any of the measurement is |1>, it applies the gate.
    qsim.cond_add(m0, X::new(q1)); // -> if Phi measured |1>, Bob applies Pauli X
    qsim.cond_add(m1, Z::new(q1)); // -> if Alice measured |1>, Bob applies Pauli Z

    qsim.state_vec(true); // -> superposition state after teleportation

    qsim.state_qubit(q1); // -> verify, now phi state should be teleported to Bob's qubit state
    qsim.state_qubit(q0); // -> verify, now Alice's qubit has collapse
}
