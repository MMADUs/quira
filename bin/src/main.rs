use core::circuit::QuantumCircuit as QC;

use core::operations::single_qubit::clifford::Hadamard as H;

fn main() {
    let mut qc = QC::new(true);

    let q0 = qc.qb_zero();
    let q1 = qc.qb_zero();
    let q2 = qc.qb_zero();
    let q3 = qc.qb_zero();

    //println!("{}, {}, {}, {}", q0, q1, q2, q3);

    qc.add_operation(H::new(q0));
    //qc.add_operation(H::new(q1));
    //qc.add_operation(H::new(q2));
    //qc.add_operation(H::new(q3));
    
    //qc.add_measurement(q1, 0);
    //qc.add_measurement(q1, 1);
    //qc.add_measurement(q2, 2);
    //qc.add_measurement(q3, 3);

    let iterations = 10000;
    let _results = qc.execute(iterations);
}
