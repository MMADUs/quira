use core::circuit::QuantumCircuit as QC;
use core::operations::single::clifford::Hadamard;
use core::operations::single::SingleQubit;
use core::types::Complex;

fn main() {
    let mut qc = QC::new(true);

    let h = Hadamard::new(0);
    let c1 = Complex::new(h.alpha_re(), h.alpha_im());
    let c2 = Complex::new(h.beta_re(), h.beta_im());
    println!("{} + {} = {}", c1, c2, c1 + c2);

    let q0 = qc.qb_zero();

    qc.add_operation(Hadamard::new(q0));

    qc.add_measurement(q0, 0);

    let iterations = 100;
    let _results = qc.execute(iterations);
}
