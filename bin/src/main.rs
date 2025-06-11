//! Simulation TEST

use quira::QuantumCircuit as QC;
use quira::QuantumSimulator as QS;
use quira::QubitIndexing as QI;

use quira::SingleQ::{Hadamard as H, PauliX as X, PauliZ as Z};
use quira::TwoQ::ControlledNot as CNOT;

use quira::density::*;

fn main() {
    test_linalg();
    let mut circuit = QC::new();

    let q = circuit.zeros(3);

    circuit.add(H::new(q[0]));
    circuit.add(H::new(q[1]));
    circuit.add(CNOT::new(q[1], q[2]));
    circuit.add(CNOT::new(q[0], q[1]));
    circuit.add(H::new(q[0]));
    circuit.measure(q[0], 0);
    circuit.measure(q[1], 1);
    circuit.cond_add(0, 1, X::new(q[2]));
    circuit.cond_add(1, 1, Z::new(q[2]));

    let mut simulator = QS::new(QI::LittleEndian);
    simulator.from_circuit(circuit);
    simulator.state_vec(true);
    simulator.state_qubit(q[2]);
    simulator.state_qubit(q[1]);
}
