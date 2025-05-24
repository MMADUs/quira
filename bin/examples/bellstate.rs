//! Bell State example
//! bellstate.rs
//!
//! run example with: cargo run --example bellstate

use quira::QuantumCircuit as QC;
use quira::Qubit;
use quira::QubitIndexing as QI;
use quira::RunResult as RR;

use quira::SingleQ::Hadamard as H;
use quira::TwoQ::ControlledNot as CNOT;

fn main() {
    // make a new quantum circuit
    let mut qc: QC = QC::new(true, QI::LittleEndian); // -> Qiskit uses Little Endian.

    // initialize zero state qubit
    let q0: Qubit = qc.qb_zero();
    let q1: Qubit = qc.qb_zero();

    qc.state_vector(false);
    // [0] |00⟩: (amp = 1.0000 + 0.0000i) => (prob = 1.0000, 100.00%)
    // [1] |01⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [2] |10⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)
    // [3] |11⟩: (amp = 0.0000 + 0.0000i) => (prob = 0.0000, 0.00%)

    // apply hadamard & controlled not
    qc.add(H::new(q0)); // -> superposition
    qc.add(CNOT::new(q0, q1)); // -> entanglement

    qc.state_vector(true);
    // [0] |00⟩: (amp = 0.7071 + 0.0000i) => (prob = 0.5000, 50.00%)
    // [3] |11⟩: (amp = 0.7071 + 0.0000i) => (prob = 0.5000, 50.00%)

    let m0: bool = qc.measure(q0, 0); // -> measure q0 to classical reg 0
    let m1: bool = qc.measure(q1, 1); // -> measure q1 to classical reg 1

    println!("m0 = {}, m1 = {}", m0, m1); // -> superposition outcome

    let shots: usize = 10000;
    let result: RR = qc.run(shots); // -> run circuit for N shots

    qc.analyze(result); // -> analyze the statistic outcome from N shots result
}
