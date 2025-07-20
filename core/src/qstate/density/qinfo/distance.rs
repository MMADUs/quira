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

use ndarray_linalg::{SVD, Scalar, Trace};

use crate::kernel::density::matrix::Density;

/// Trace distance (1/2) * ||ρ - σ||₁
pub fn trace_distance(rho: &Density, alpha: &Density) -> f64 {
    assert_eq!(rho.dim(), alpha.dim(), "dimensions must match");
    let diff = rho.matrix_as_ref() - alpha.matrix_as_ref();
    let svd = diff.svd(true, true).expect("SVD failed");
    let singular_values = svd.1;
    let trace_norm: f64 = singular_values.sum();
    0.5 * trace_norm
}

/// Fidelity F(ρ, σ) = [Tr(√(√ρ σ √ρ))]^2
pub fn fidelity(rho: &Density, alpha: &Density) -> f64 {
    assert_eq!(rho.dim(), alpha.dim(), "dimensions must match");

    if rho.is_pure() && alpha.is_pure() {
        // For pure states: F(ρ, σ) = Tr(ρσ) = |⟨ψ|φ⟩|²
        let product = rho.matrix_as_ref().dot(alpha.matrix_as_ref());
        let trace_val = product.trace().unwrap();
        trace_val.norm_sqr()
    } else {
        // General case: Uhlmann fidelity
        let sqrt_rho = crate::utils::eigen::matrix_sqrt(rho.matrix_as_ref());
        let product = sqrt_rho.dot(alpha.matrix_as_ref()).dot(&sqrt_rho);
        let eigenvals = crate::utils::eigen::eigen_values(&product);
        // Fidelity = (∑ sqrt(λᵢ))²
        let trace_sqrt = eigenvals
            .iter()
            .map(|val| val.re().max(0.0).sqrt()) // only real, non-negative eigenvalues
            .sum::<f64>();
        trace_sqrt * trace_sqrt
    }
}

/// Bures distance d_B(ρ, σ) = √(2(1 - √F(ρ, σ)))
pub fn bures_distance(rho: &Density, alpha: &Density) -> f64 {
    let fidelity = fidelity(rho, alpha);
    (2.0 * (1.0 - fidelity.sqrt())).sqrt()
}

/// Hilbert-Schmidt distance ||ρ - σ||_HS = √Tr((ρ - σ)²)
pub fn hilbert_schmidt_distance(rho: &Density, alpha: &Density) -> f64 {
    assert_eq!(rho.dim(), alpha.dim(), "dimensions must match");

    let diff = rho.matrix_as_ref() - alpha.matrix_as_ref();
    let diff_squared = diff.dot(&diff);
    diff_squared.diag().sum().re.sqrt()
}

/// Operator norm (spectral norm) ||ρ - σ||_∞
pub fn operator_norm_distance(rho: &Density, alpha: &Density) -> f64 {
    assert_eq!(rho.dim(), alpha.dim(), "dimensions must match");

    let diff = rho.matrix_as_ref() - alpha.matrix_as_ref();
    let eigenvals = crate::utils::eigen::eigen_values(&diff);
    eigenvals.iter().map(|&ev| ev.abs()).fold(0.0, f64::max)
}
