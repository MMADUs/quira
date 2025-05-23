//! Simulation TEST

use quira::QuantumCircuit as QC;
use quira::QubitIndexing as QI;
use quira::SingleQ::Hadamard as H;
use quira::TwoQ::ControlledNot as CNOT;

fn main() {
    let mut qc = QC::new(true, QI::LittleEndian);

    let q0 = qc.qb_zero();
    let q1 = qc.qb_zero();
    //let q2 = qc.qb_zero();
    //let q3 = qc.qb_zero();

    qc.add_operation(H::new(q0));
    qc.add_operation(CNOT::new(q0, q1));
    //qc.add_operation(H::new(q1));
    //qc.add_operation(H::new(q2));
    //qc.add_operation(H::new(q3));

    qc.state_vector(true);

    qc.add_measurement(q0, 0);
    qc.add_measurement(q1, 1);
    //qc.add_measurement(q2, 2);
    //qc.add_measurement(q3, 3);
    
    let iterations = 10000;
    let results = qc.execute(iterations);
    qc.count_outcomes(results);
}
