//! Copyright (c) 2024-2025 Quira, Inc.
//!
//! This file is part of Quira
//!
//! This program is free software: you can redistribute it and/or modify
//! it under the terms of the GNU Affero General Public License as published by
//! the Free Software Foundation, either version 3 of the License, or
//! (at your option) any later version.
//!
//! This program is distributed in the hope that it will be useful
//! but WITHOUT ANY WARRANTY; without even the implied warranty of
//! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//! GNU Affero General Public License for more details.
//!
//! You should have received a copy of the GNU Affero General Public License
//! along with this program.  If not, see <http://www.gnu.org/licenses/>.

use ndarray_linalg::{EigValsh, Eigh, Inverse, UPLO};

use crate::{
    constant::{EPSILON, INF},
    types::{Complex, Matrix, Vector},
};

pub fn eigen_rank(matrix: &Matrix<Complex>) -> usize {
    let eigenvals = eigen_values(matrix);
    eigenvals.iter().filter(|&&ev| ev > EPSILON).count()
}

pub fn eigen_values(matrix: &Matrix<Complex>) -> Vector<f64> {
    matrix
        .eigvalsh(UPLO::Lower)
        .expect("Failed to compute eigenvalues")
}

pub fn eigen_decomposition(matrix: &Matrix<Complex>) -> (Vector<f64>, Matrix<Complex>) {
    matrix
        .eigh(UPLO::Lower)
        .expect("Failed to compute eigen decomposition")
}

pub fn matrix_sqrt(matrix: &Matrix<Complex>) -> Matrix<Complex> {
    let (eigenvals, eigenvecs) = eigen_decomposition(matrix);

    let sqrt_eigenvals: Vector<Complex> =
        eigenvals.mapv(|val| Complex::new(val.clamp(0.0, INF).sqrt(), 0.0));

    let sqrt_diag = Matrix::from_diag(&sqrt_eigenvals);
    eigenvecs
        .dot(&sqrt_diag)
        .dot(&eigenvecs.t().mapv(|x| x.conj()))
}

pub fn matrix_log(matrix: &Matrix<Complex>) -> Option<Matrix<Complex>> {
    let (eigenvalues, eigenvectors) = eigen_decomposition(matrix);
    let dim = matrix.nrows();

    // Diagonal matrix of log(eigenvalues)
    let mut log_diag = Matrix::<Complex>::zeros((dim, dim));
    for i in 0..dim {
        if eigenvalues[i] <= 0.0 {
            return None; // log undefined or divergent
        }
        log_diag[[i, i]] = Complex::new(eigenvalues[i].ln(), 0.0);
    }

    let vinv = eigenvectors.inv().ok()?;
    let result = eigenvectors.dot(&log_diag).dot(&vinv);
    Some(result)
}
