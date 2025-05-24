//! Simulation TEST

use quira::QuantumCircuit as QC;
use quira::QubitIndexing as QI;

fn main() {
    let mut qc = QC::new(true, QI::LittleEndian);

    qc.zeros(10);
    qc.state_vector(false);
}
