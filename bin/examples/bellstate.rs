//! Bell State example
//! bellstate.rs
//!
//! run example with: cargo run --example bellstate

use quira::include::{ClassicalRegister, QuantumCircuit, QuantumJob, QuantumRegister, StateVec};
use quira::operation::singleq::Hadamard as H;
use quira::operation::twoq::ControlledNot as CNOT;
use quira::provider::QubitIndexing as QI;

fn main() {
    // define register
    //
    let qreg = QuantumRegister::new(2).label("qreg 1");
    let creg = ClassicalRegister::new(2).label("creg 1");
    // define circuit
    //
    let mut circuit = QuantumCircuit::with_reg(&qreg, &creg);
    // build circuit
    //
    circuit.add(H::new(&qreg[0]));
    circuit.add(CNOT::new(&qreg[0], &qreg[1]));
    circuit.measure_all();
    // define backend
    //
    let backend = StateVec::new();
    // define job
    //
    let mut job = QuantumJob::new(backend, QI::LittleEndian);
    job.from_circuit(circuit);
    // run job
    //
    job.run(1000000);
    // inspect job result
    //
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

    println!("hello quantum world!");
}
