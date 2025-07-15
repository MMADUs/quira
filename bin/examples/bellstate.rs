//! Bell State example
//! bellstate.rs
//!
//! run example with: cargo run --example bellstate

use quira::QuantumCircuit;
use quira::QuantumJob;
use quira::QuantumSimulator as QuiraSimualtor;
use quira::QubitIndexing as QI;
use quira::kernel::QuantumDebugger;

use quira::SingleQ::Hadamard as H;
use quira::StateVec;
use quira::TwoQ::ControlledNot as CNOT;

fn main() {
    // build quantum circuit
    //
    let mut circuit = QuantumCircuit::new(2, 1);

    let q0 = circuit.qb_zero();
    let q1 = circuit.qb_zero();

    circuit.add(H::new(q0));
    circuit.add(CNOT::new(q0, q1));

    // simulating quantum state
    //
    let backend = StateVec::new();
    let mut sim = QuiraSimualtor::new(backend, QI::LittleEndian);
    let state = sim.from_circuit(circuit).backend();

    // debug state on simulation
    //
    for (bit, amplitude) in state.entire_state(false) {
        println!("|{}>: {}", bit, amplitude);
    }

    let qs0 = state.qubit_state(q0);
    let qs1 = state.qubit_state(q1);

    println!("qs0 = {}|0> + {}|1>", qs0.0, qs0.1);
    println!("qs1 = {}|0> + {}|1>", qs1.0, qs1.1);

    let bs = state.bit_state(&[false, false]); // |00>
    println!("|00>: {}", bs);

    // run measurement outcome
    //
    let mut circuit = QuantumCircuit::new(2, 1);

    let q0 = circuit.qb_zero();
    let q1 = circuit.qb_zero();

    circuit.add(H::new(q0));
    circuit.add(CNOT::new(q0, q1));
    circuit.measure_all();

    let backend = StateVec::new();
    let mut job = QuantumJob::new(backend, QI::LittleEndian);
    job.from_circuit(circuit);
    job.run(1000);
    let results = job.result();
    for result in results.unwrap() {
        for (bit, count) in result.get_counts() {
            println!("{}: {}", bit, count);
        }

        let (mb, mc) = result.get_most_frequent().unwrap();
        println!("|{}>: {}", mb, mc);

        result.print_results();
        println!("{}", result.histogram(20));
    }

    println!("hello quantum world!");
}
