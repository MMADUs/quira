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

use ndarray_linalg::Scalar;

use crate::{
    constant::{EPSILON, INF},
    eigen,
};

use super::matrix::Density;

pub struct Entropy {
    density: Density,
}

impl Entropy {
    pub fn from(density: Density) -> Self {
        Self { density }
    }

    /// Shannon entropy in the computational basis (diagonal elements)
    pub fn shannon_entropy(&self) -> f64 {
        let diag = self.density.matrix_as_ref().diag();
        diag.iter()
            .map(|c| c.re()) // Extract real part
            .filter(|&p| p > EPSILON) // Avoid log(0)
            .map(|p| -p * p.ln())
            .sum()
    }

    /// Compute linear entropy 1 - Tr(ρ²)
    pub fn linear_entropy(&self) -> f64 {
        1.0 - self.density.purity()
    }

    /// Von Neumann entropy S(ρ) = -Tr(ρ log ρ)
    pub fn von_neumann_entropy(&self) -> f64 {
        self.density.von_neumann_entropy()
    }

    /// Relative entropy S(ρ||σ) = Tr(ρ(log ρ - log σ))
    /// Note: This is a simplified implementation using diagonal approximation
    pub fn relative_entropy(&self, other: &Density) -> f64 {
        assert_eq!(self.density.dim(), other.dim(), "dimensions must match");

        let eigenvals_a = eigen::eigen_values(self.density.matrix_as_ref());
        let eigenvals_b = eigen::eigen_values(other.matrix_as_ref());

        // For a proper implementation, both matrices should be diagonalized
        // in the same basis. This is a simplified version.
        let mut relative_entropy = 0.0;
        for i in 0..eigenvals_a.len() {
            let rho_i = eigenvals_a[i];
            let sigma_i = eigenvals_b[i];
            if rho_i > EPSILON && sigma_i > EPSILON {
                relative_entropy += rho_i * (rho_i.ln() - sigma_i.ln());
            } else if rho_i > EPSILON && sigma_i <= EPSILON {
                return f64::INFINITY; // Divergence
            }
        }
        relative_entropy
    }

    /// Mutual information I(A:B) = S(ρ_A) + S(ρ_B) - S(ρ_AB)
    pub fn mutual_information(&self, subsystem_dims: &[usize]) -> f64 {
        assert_eq!(
            subsystem_dims.len(),
            2,
            "mutual information requires exactly 2 subsystems"
        );

        let rho_a = self.density.partial_trace(subsystem_dims, &[1]);
        let rho_b = self.density.partial_trace(subsystem_dims, &[0]);

        let s_a = rho_a.von_neumann_entropy();
        let s_b = rho_b.von_neumann_entropy();
        let s_ab = self.von_neumann_entropy();

        s_a + s_b - s_ab
    }

    /// Conditional entropy S(A|B) = S(ρ_AB) - S(ρ_B)
    pub fn conditional_entropy(&self, subsystem_dims: &[usize], conditioned_on: usize) -> f64 {
        assert_eq!(
            subsystem_dims.len(),
            2,
            "conditional entropy requires exactly 2 subsystems"
        );
        assert!(conditioned_on < 2, "conditioned_on must be 0 or 1");

        let trace_out = if conditioned_on == 0 {
            vec![1]
        } else {
            vec![0]
        };

        let rho_conditioned = self.density.partial_trace(subsystem_dims, &trace_out);
        let s_joint = self.von_neumann_entropy();
        let s_conditioned = rho_conditioned.von_neumann_entropy();

        s_joint - s_conditioned
    }

    /// Entanglement entropy (same as von Neumann entropy for bipartite systems)
    pub fn entanglement_entropy(&self, subsystem_dims: &[usize], subsystem_idx: usize) -> f64 {
        assert!(
            subsystem_idx < subsystem_dims.len(),
            "subsystem index out of bounds"
        );

        let mut trace_over = Vec::new();
        for i in 0..subsystem_dims.len() {
            if i != subsystem_idx {
                trace_over.push(i);
            }
        }

        let reduced_density = self.density.partial_trace(subsystem_dims, &trace_over);
        reduced_density.von_neumann_entropy()
    }

    /// Renyi entropy of order α: S_α(ρ) = (1/(1-α)) * log(Tr(ρ^α))
    pub fn renyi_entropy(&self, alpha: f64) -> f64 {
        assert!(
            alpha > 0.0 && alpha != 1.0,
            "alpha must be positive and not equal to 1"
        );

        if (alpha - 1.0).abs() < EPSILON {
            return self.von_neumann_entropy();
        }

        let eigenvals = eigen::eigen_values(self.density.matrix_as_ref());
        let mut trace_rho_alpha = 0.0;

        for &eigenval in &eigenvals {
            if eigenval > EPSILON {
                trace_rho_alpha += eigenval.powf(alpha);
            }
        }

        if trace_rho_alpha <= EPSILON {
            return INF;
        }

        (1.0 / (1.0 - alpha)) * trace_rho_alpha.ln()
    }

    pub fn tsallis_entropy(&self, q: f64) -> f64 {
        assert!(
            q > 0.0 && (q - 1.0).abs() > EPSILON,
            "q must be > 0 and ≠ 1"
        );

        let eigenvals = eigen::eigen_values(self.density.matrix_as_ref());
        let trace_rho_q: f64 = eigenvals
            .iter()
            .filter(|&&λ| λ > EPSILON)
            .map(|&λ| λ.powf(q))
            .sum();

        (1.0 - trace_rho_q) / (q - 1.0)
    }

    /// Min entropy (Renyi entropy of order infinity)
    pub fn min_entropy(&self) -> f64 {
        let eigenvals = eigen::eigen_values(self.density.matrix_as_ref());
        let max_eigenval = eigenvals.iter().fold(0.0f64, |acc, &x| acc.max(x));

        if max_eigenval <= EPSILON {
            INF
        } else {
            -1.0 * max_eigenval.ln()
        }
    }

    /// Max entropy (Renyi entropy of order 0)
    pub fn max_entropy(&self) -> f64 {
        let rank = eigen::eigen_rank(self.density.matrix_as_ref());
        (rank as f64).ln()
    }
}
