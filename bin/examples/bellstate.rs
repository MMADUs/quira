//! Bell State example
//! bellstate.rs
//!
//! run example with: cargo run --example bellstate

use quira::QuantumSimulator as QS;
use quira::Qubit;
use quira::QubitIndexing as QI;

use quira::SingleQ::Hadamard as H;
use quira::TwoQ::ControlledNot as CNOT;

fn main() {
    // make a new quantum circuit
    let mut qsim: QS = QS::new(QI::LittleEndian); // -> Qiskit uses Little Endian.

    // initialize zero state qubit
    let q0: Qubit = qsim.qb_zero();
    let q1: Qubit = qsim.qb_zero();

    qsim.state_vec(false);
    // [0] |00⟩: (amp = 1.0000 + 0.0000i) => (prob = 1.0000, 100.00%)
    // [1] |01⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [2] |10⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [3] |11⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)

    // apply hadamard & controlled not
    qsim.add(H::new(q0)); // -> superposition
    qsim.add(CNOT::new(q0, q1)); // -> entanglement

    qsim.state_vec(true);
    // [0] |00⟩: (amp = 0.7071 + 0.0000i) => (prob = 0.5000, 50.00%)
    // [3] |11⟩: (amp = 0.7071 + 0.0000i) => (prob = 0.5000, 50.00%)

    let m0: bool = qsim.measure(q0, 0); // -> measure q0 to classical reg 0
    let m1: bool = qsim.measure(q1, 1); // -> measure q1 to classical reg 1

    println!("m0 = {}, m1 = {}", m0 as u8, m1 as u8); // -> superposition outcome
}
