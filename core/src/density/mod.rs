use ndarray::array;
use ndarray_linalg::Inverse;

pub fn test_linalg() {
    let a = array![[1.0, 2.0], [3.0, 4.0],];

    let inv = a.inv().expect("matrix inversion failed");
    println!("Inverse:\n{:?}", inv);
}
