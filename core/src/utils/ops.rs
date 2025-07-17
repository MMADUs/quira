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

use ndarray::LinalgScalar;

use crate::types::{Complex, Matrix};

/// Computes the Kronecker product of two square matrices a and b.
/// Both matrices must be square of any dimension.
/// Returns a matrix of size (a_dim * b_dim) x (a_dim * b_dim).
pub fn kron<T>(a: &Matrix<T>, b: &Matrix<T>) -> Matrix<T>
where
    T: LinalgScalar,
{
    let dim_a = a.shape()[0];
    let dim_b = b.shape()[0];
    assert_eq!(a.shape()[1], dim_a, "Matrix a must be square");
    assert_eq!(b.shape()[1], dim_b, "Matrix b must be square");
    let dim = dim_a * dim_b;
    let mut result = Matrix::zeros((dim, dim));
    for (mut chunk, elem) in result
        .exact_chunks_mut((dim_b, dim_b))
        .into_iter()
        .zip(a.iter())
    {
        let v: Matrix<T> = Matrix::from_elem((dim_b, dim_b), *(elem)) * b;
        chunk.assign(&v);
    }
    result
}

/// Matrix dagger (conjugate tranpose) operation.
pub fn dagger(u: &Matrix<Complex>) -> Matrix<Complex> {
    u.mapv(|x| x.conj()).t().to_owned()
}

// Matrix partial transpose operation.
pub fn partial_transpose<T>(matrix: &Matrix<T>, dim_a: usize, dim_b: usize) -> Matrix<T>
where
    T: LinalgScalar,
{
    let dim = dim_a * dim_b;
    let (rows, cols) = matrix.dim();
    assert_eq!((rows, cols), (dim, dim), "dimensions must match");
    let mut result = Matrix::<T>::zeros((dim, dim));
    for a1 in 0..dim_a {
        for b1 in 0..dim_b {
            for a2 in 0..dim_a {
                for b2 in 0..dim_b {
                    let from_row = a1 * dim_b + b1;
                    let from_col = a2 * dim_b + b2;
                    let to_row = a1 * dim_b + b2;
                    let to_col = a2 * dim_b + b1;

                    result[[to_row, to_col]] = matrix[[from_row, from_col]];
                }
            }
        }
    }
    result
}
