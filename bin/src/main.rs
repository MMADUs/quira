use core::circuit::QuantumCircuit as QC;
use core::operations::single::clifford::Hadamard;

fn main() {
    let mut qc = QC::new(true);

    let q0 = qc.qb_zero();

    qc.add_operation(Hadamard::new(q0));

    qc.add_measurement(q0, 0);

    let iterations = 10000000;
    let _results = qc.execute(iterations);
}
