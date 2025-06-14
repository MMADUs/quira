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

use crate::{
    ops,
    types::{Complex, Matrix},
};

/// Qubit index ordering mode when applying unitary to quantum state.
#[derive(Debug, Clone)]
pub enum QubitIndexing {
    LittleEndian,
    BigEndian,
}

/// Builds a permutation matrix P of size 2^n x 2^n
/// that moves the qubits in targets to the front (lowest bits)
/// while preserving the order of other qubits.
pub(crate) fn permutation_matrix(n: usize, targets: &[usize]) -> Matrix<Complex> {
    let dim = 1 << n;
    let mut p = Matrix::<Complex>::zeros((dim, dim));
    for i in 0..dim {
        // Extract bits of i as bool vector
        let bits: Vec<bool> = (0..n).map(|bit| (i >> bit) & 1 == 1).collect();
        // Rearrange bits: targets first, then others
        let mut permuted_bits = Vec::with_capacity(n);
        // Add target bits
        for &t in targets {
            permuted_bits.push(bits[t]);
        }
        // Add remaining bits in ascending order
        for pos in 0..n {
            if !targets.contains(&pos) {
                permuted_bits.push(bits[pos]);
            }
        }
        // Convert permuted bits back to integer index
        let mut j = 0;
        for (idx, &b) in permuted_bits.iter().enumerate() {
            if b {
                j |= 1 << idx;
            }
        }
        p[(j, i)] = Complex::new(1.0, 0.0);
    }
    p
}

/// Expanding a k-qubit unitary into an n-qubit quantum state (Hilbert space),
/// acting on the qubits specified in targets.
///
/// Let:
/// - n = total number of qubits
/// - k = number of qubits that the gate U acts on
/// - U = 2^k x 2^k unitary matrix
/// - I = identity matrix of size 2^(n-k) x 2^(n-k)
/// - P = permutation matrix that moves target qubits to the front (lowest bits)
/// - P† = the conjugate transpose (inverse) of P.
///
/// Then the expanded operator is:
///
/// For BigEndian indexing:
///     U_expanded = P† * ( U ⊗ I ⊗ ... ⊗ I ) * P
///
/// For LittleEndian indexing:
///     U_expanded = P† * ( I ⊗ ... ⊗ I ⊗ U ) * P
///
pub(crate) fn expand_unitary(
    n: usize,
    targets: &[usize],
    u: &Matrix<Complex>,
    indexing: &QubitIndexing,
) -> Matrix<Complex> {
    if targets.is_empty() {
        return Matrix::eye(1 << n); // do nothing
    }

    let k = targets.len();
    assert_eq!(u.shape(), &[1 << k, 1 << k], "Unitary must be 2^k x 2^k");

    let p = permutation_matrix(n, &targets);
    let p_dagger = ops::dagger(&p);

    let i = Matrix::<Complex>::eye(1 << (n - k));

    let u_kron = match indexing {
        QubitIndexing::BigEndian => ops::kron(&u, &i),
        QubitIndexing::LittleEndian => ops::kron(&i, &u),
    };

    p_dagger.dot(&u_kron).dot(&p)
}
