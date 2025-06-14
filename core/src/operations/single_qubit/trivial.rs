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

use super::SingleQubitGate;
use crate::operations::QuantumGate;
use crate::types::{Complex, Matrix, Qubit};
use ndarray::array;

#[derive(Debug, Clone)]
/// Represents the Identity gate, which leaves the qubit state unchanged.
///
/// The Identity gate is a fundamental quantum gate that performs no operation on the qubit state.
/// It is equivalent to a "no operation" instruction and is often used in quantum circuit construction
/// and as a placeholder when no action is needed on a particular qubit.
///
/// The matrix form is:
///
/// I = [ [ 1, 0 ],
///       [ 0, 1 ] ]
///
/// This gate leaves the quantum state completely unchanged. It's useful in quantum algorithms,
/// error correction, and circuit design for padding or alignment purposes.
pub struct Identity;

impl Identity {
    pub fn new() -> Self {
        Self
    }
}

impl QuantumGate for Identity {
    /// constructs the 2x2 unitary matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("I")
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![]
    }
}

impl SingleQubitGate for Identity {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        1.0
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        0.0
    }

    /// returns the real part of beta
    fn beta_re(&self) -> f64 {
        0.0
    }

    /// returns the imaginary part of beta
    fn beta_im(&self) -> f64 {
        0.0
    }

    /// returns the global phase
    fn global_phase(&self) -> f64 {
        0.0
    }
}
