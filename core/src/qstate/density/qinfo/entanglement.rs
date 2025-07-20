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

use ndarray_linalg::SVD;

use crate::bit::QuantumBit;
use crate::constant::EPSILON;
use crate::kernel::density::matrix::Density;
use crate::operations::singleq::PauliY;
use crate::ops::QuantumGate;
use crate::types::{Complex, Matrix, Vector};
use crate::utils::{eigen, ops};

/// Concurrence for two-qubit states
pub fn concurrence(density: &Density) -> f64 {
    assert_eq!(density.dim(), 2, "concurrence is only for 2 qubit system");
    // Y^⊗n where n is n-qubit system
    let pauli_y = PauliY::new(&QuantumBit::new(0)).unitary_matrix();
    let y_tensor_y = ops::kron(&pauli_y, &pauli_y);
    // ρ̃ = (Y^⊗n) ρ* (Y^⊗n)
    let rho_conj = density.matrix_as_ref().mapv(|x| x.conj());
    let rho_tilde = y_tensor_y.dot(&rho_conj).dot(&y_tensor_y);
    // R = ρ ρ̃
    let r_matrix = density.matrix_as_ref().dot(&rho_tilde);
    let eigenvals: Vector<f64> = eigen::eigen_values(&r_matrix);
    // Convert Vector<f64> to Vec<f64> for sorting
    let mut eigenvals_vec: Vec<f64> = eigenvals.to_vec();
    // Sort eigenvalues in descending order
    eigenvals_vec.sort_by(|a, b| b.partial_cmp(a).unwrap());
    // Take square roots
    let sqrt_eigenvals: Vec<f64> = eigenvals_vec
        .iter()
        .map(|&x| if x >= 0.0 { x.sqrt() } else { 0.0 })
        .collect();
    // Concurrence = max(0, λ₁ - λ₂ - λ₃ - λ₄)
    if sqrt_eigenvals.len() >= 4 {
        let concurrence =
            sqrt_eigenvals[0] - sqrt_eigenvals[1] - sqrt_eigenvals[2] - sqrt_eigenvals[3];
        concurrence.max(0.0)
    } else {
        0.0
    }
}

/// Entanglement of formation (for two-qubit states)
pub fn entanglement_of_formation(density: &Density) -> f64 {
    let c = concurrence(density);
    if c == 0.0 {
        return 0.0;
    }
    // make as internal function
    let h = |x: f64| {
        if x == 0.0 || x == 1.0 {
            0.0
        } else {
            -1.0 * x * x.log2() - (1.0 - x) * (1.0 - x).log2()
        }
    };
    // compute x
    let x = (1.0 + (1.0 - c * c).sqrt()) / 2.0;
    h(x)
}

/// Logarithmic negativity (for bipartite systems)
pub fn logarithmic_negativity(density: &Density, subsystem_dims: &[usize]) -> f64 {
    assert_eq!(
        subsystem_dims.len(),
        2,
        "logarithmic negativity requires exactly 2 subsystems"
    );
    // Partial transpose with respect to the second subsystem
    let rho_pt = ops::partial_transpose(
        density.matrix_as_ref(),
        subsystem_dims[0],
        subsystem_dims[1],
    );
    // Compute eigenvalues of the partial transpose
    let eigenvals = eigen::eigen_values(&rho_pt);
    // Trace norm = sum of absolute values of all eigenvalues
    let trace_norm: f64 = eigenvals.iter().map(|&ev| ev.abs()).sum();
    trace_norm.log2()
}

/// Negativity (for bipartite systems)
pub fn negativity(density: &Density, subsystem_dims: &[usize]) -> f64 {
    assert_eq!(
        subsystem_dims.len(),
        2,
        "negativity requires exactly 2 subsystems"
    );
    // Partial transpose of the density matrix
    let rho_pt = ops::partial_transpose(
        density.matrix_as_ref(),
        subsystem_dims[0],
        subsystem_dims[1],
    );
    // Eigenvalues of the partially transposed matrix
    let eigenvals = eigen::eigen_values(&rho_pt);
    // Sum of the absolute values of the negative eigenvalues
    eigenvals
        .iter()
        .filter(|&&ev| ev < 0.0)
        .map(|&ev| ev.abs())
        .sum()
}

// Schmidt decomposition for bipartite pure states
// Returns (Schmidt coefficients, left vectors, right vectors)
pub fn schmidt_decomposition(
    state_vector: &Vector<Complex>,
    dim_a: usize,
    dim_b: usize,
) -> Option<(Vec<f64>, Vec<Vector<Complex>>, Vec<Vector<Complex>>)> {
    assert_eq!(
        state_vector.len(),
        dim_a * dim_b,
        "State vector size must match dim_a * dim_b"
    );

    // Reshape state vector into a dim_a x dim_b matrix
    let mut state_matrix = Matrix::<Complex>::zeros((dim_a, dim_b));
    for i in 0..dim_a {
        for j in 0..dim_b {
            state_matrix[[i, j]] = state_vector[i * dim_b + j];
        }
    }

    // Perform SVD: A = U Σ V†
    let svd_result = state_matrix.svd(true, true).ok()?;
    let (u_opt, sigma, vt_opt) = svd_result;

    let u = u_opt?;
    let vt = vt_opt?;

    let mut schmidt_coeffs = Vec::new();
    let mut left_vecs = Vec::new();
    let mut right_vecs = Vec::new();

    for i in 0..sigma.len() {
        let s = sigma[i];
        if s > EPSILON {
            schmidt_coeffs.push(s);

            // Left singular vector (column of U)
            let u_vec = u.column(i).to_owned();
            left_vecs.push(u_vec);

            // Right singular vector (row of V†, take conjugate transpose)
            let v_row = vt.row(i);
            let v_col = v_row.mapv(|x| x.conj());
            right_vecs.push(v_col);
        }
    }

    Some((schmidt_coeffs, left_vecs, right_vecs))
}

// Helper function to extract state vector from pure density matrix
fn extract_pure_state_vector(density: &Density) -> Option<Vector<Complex>> {
    if !density.is_pure() {
        return None;
    }

    // Find the eigenvalue closest to 1 and its eigenvector
    let (eigenvals, eigenvecs) = eigen::eigen_decomposition(density.matrix_as_ref());

    let mut max_idx = 0;
    let mut max_val = eigenvals[0];

    for (i, &val) in eigenvals.iter().enumerate() {
        if val > max_val {
            max_val = val;
            max_idx = i;
        }
    }

    if (max_val - 1.0).abs() < EPSILON {
        Some(eigenvecs.column(max_idx).to_owned())
    } else {
        None
    }
}

/// Schmidt rank (number of non-zero Schmidt coefficients)
pub fn schmidt_rank(density: &Density, dim_a: usize, dim_b: usize) -> Option<usize> {
    let state_vector = extract_pure_state_vector(density)?;
    schmidt_decomposition(&state_vector, dim_a, dim_b).map(|(coeffs, _, _)| coeffs.len())
}

/// Schmidt number (effective number of terms in Schmidt decomposition)
pub fn schmidt_number(density: &Density, dim_a: usize, dim_b: usize) -> Option<f64> {
    let state_vector = extract_pure_state_vector(density)?;
    schmidt_decomposition(&state_vector, dim_a, dim_b).map(|(coeffs, _, _)| {
        let sum_squares: f64 = coeffs.iter().map(|&c| c.powi(4)).sum();
        if sum_squares > EPSILON {
            1.0 / sum_squares
        } else {
            0.0
        }
    })
}
