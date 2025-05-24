//! Simulation TEST

use quira::QuantumCircuit as QC;
use quira::QubitIndexing as QI;
use quira::SingleQ::{Hadamard as H, PauliX as X, PauliZ as Z};
use quira::TwoQ::ControlledNot as CNOT;

fn main() {
    let mut qc = QC::new(true, QI::LittleEndian);

    let phi = qc.qb_zero();
    let q0 = qc.qb_zero();
    let q1 = qc.qb_zero();

    qc.add(H::new(phi));

    qc.state_vector(false);

    qc.add(H::new(q0));
    qc.add(CNOT::new(q0, q1));

    qc.state_vector(false);

    qc.add(CNOT::new(phi, q0));
    qc.add(H::new(phi));

    qc.state_vector(true);

    let m0 = qc.measure(phi, 0);
    let m1 = qc.measure(q0, 1);

    qc.state_vector(true);

    qc.cond_add(m0, X::new(q1));
    qc.cond_add(m1, Z::new(q1));

    println!("m0 = {}, m1 = {}", m0, m1);

    qc.state_vector(true);

    qc.state_qubit(q1);

    let shots = 10000;
    let result = qc.run(shots);

    qc.analyze(result);
}
