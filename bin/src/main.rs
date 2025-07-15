//! Simulation TEST

use quira::kernel::QuantumDebugger;
use quira::*;
use quira::QuantumCircuit;
use quira::QuantumSimulator as QuiraSimulator;
use quira::QubitIndexing as QI;

use quira::SingleQ::{Hadamard as H, PauliX as X, PauliZ as Z};
use quira::TwoQ::ControlledNot as CNOT;

fn main() {
    let mut circuit_a = QuantumCircuit::new(1, 1);
    let q0 = circuit_a.qb_zero();
    circuit_a.add(H::new(q0));

    let mut circuit_b = QuantumCircuit::new(0, 0);
    circuit_b.add(H::new(q0));

    circuit_a.extend(circuit_b);

    let backend = StateVec::new();
    let mut sim = QuiraSimulator::new(backend, QI::LittleEndian);
    let state = sim.from_circuit(circuit_a).backend();

    for (bit, amplitude) in state.entire_state(false) {
        println!("|{}>: {}", bit, amplitude);
    }
}

// todo list
// 1. refactor gate
// 2. stabilizer & integration
// 3. refactor noise channel
// 4. quantum fourier transform
// 5. back to the list on obsidian notes.
