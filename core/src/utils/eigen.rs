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

use ndarray_linalg::{EigValsh, Eigh, UPLO};

use crate::{
    constant::EPSILON,
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
