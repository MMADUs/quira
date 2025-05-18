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

use super::{SingleQubit, SingleQubitGate};
use crate::operations::QuantumGate;
use crate::types::{Complex, Matrix, Qubit};
use ndarray::array;

#[derive(Debug, Clone)]
/// Represents an arbitrary single-qubit quantum gate.
///
/// This structure allows the creation of any unitary 2×2 matrix for single-qubit operations
/// by specifying the complex values alpha and beta, along with a global phase factor.
///
/// The matrix form is:
///
/// U = e^(i*global_phase) * [ [ α, -β*],
///                            [ β,  α*] ]
///
/// where α* and β* are complex conjugates of α and β respectively.
/// The parameters must satisfy |α|² + |β|² = 1 to ensure unitarity.
pub struct Arbitrary {
    qubit: Qubit,
    alpha_re: f64,
    alpha_im: f64,
    beta_re: f64,
    beta_im: f64,
    global_phase: f64,
}

impl Arbitrary {
    pub fn new(
        qubit: Qubit,
        alpha_re: f64,
        alpha_im: f64,
        beta_re: f64,
        beta_im: f64,
        global_phase: f64,
    ) -> Self {
        // check if all parameters are zero
        if alpha_re == 0.0 && alpha_im == 0.0 && beta_re == 0.0 && beta_im == 0.0 {
            panic!("Invalid gate parameters: all alpha and beta components are zero");
        }

        // calculate the squared norm to verify unitarity
        let norm: f64 =
            alpha_re.powf(2.0) + alpha_im.powf(2.0) + beta_re.powf(2.0) + beta_im.powf(2.0);

        // verify that the norm is approximately 1 (unitarity condition)
        if (norm - 1.0).abs() > 1e-6 {
            panic!(
                "Unitarity violation: |α|² + |β|² must equal 1.0, got {}",
                norm
            );
        }

        Self {
            qubit,
            alpha_re,
            alpha_im,
            beta_re,
            beta_im,
            global_phase,
        }
    }
}

impl QuantumGate for Arbitrary {
    /// construct the 2x2 unitary matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        // calculate e^(i*global_phase) to get the global phase prefix
        let pref = Complex::new(0.0, self.global_phase).exp();
        array![
            [
                pref * Complex::new(self.alpha_re, self.alpha_im),
                pref * Complex::new(-1.0 * self.beta_re, self.beta_im)
            ],
            [
                pref * Complex::new(self.beta_re, self.beta_im),
                pref * Complex::new(self.alpha_re, -1.0 * self.alpha_im)
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        String::from("AR")
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target_qubit()]
    }
}

impl SingleQubit for Arbitrary {
    /// returns the index of the qubit this gate operates on.
    fn target_qubit(&self) -> Qubit {
        self.qubit
    }
}

impl SingleQubitGate for Arbitrary {
    /// returns the arbitrary real part of alpha
    fn alpha_re(&self) -> f64 {
        self.alpha_re
    }

    /// returns the arbitrary imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        self.alpha_im
    }

    /// returns the arbitrary real part of beta
    fn beta_re(&self) -> f64 {
        self.beta_re
    }

    /// returns the arbitrary imaginary part of beta
    fn beta_im(&self) -> f64 {
        self.beta_im
    }

    /// returns the arbitrary global phase
    fn global_phase(&self) -> f64 {
        self.global_phase
    }
}

#[derive(Debug, Clone)]
/// Represents the universal single-qubit gate (U gate) that can be used to implement any single-qubit operation.
///
/// The U gate is parameterized by three angles (theta, phi, lambda) and can represent any single-qubit unitary
/// transformation. It is considered a universal single-qubit gate because any other single-qubit gate can be
/// expressed in terms of the U gate with appropriate parameter choices.
///
/// The matrix form is:
///
/// U(θ, φ, λ) = [ [ cos(θ/2)                   , -e^(iλ)*sin(θ/2)    ],
///                [ e^(iφ)*sin(θ/2)            , e^(i(φ+λ))*cos(θ/2) ] ]
///
/// Which expands to:
///
/// U(θ, φ, λ) = [ [ cos(θ/2)                   , -sin(θ/2)(cos(λ) + i*sin(λ))    ],
///                [ sin(θ/2)(cos(φ) + i*sin(φ)), cos(θ/2)(cos(φ+λ) + i*sin(φ+λ)) ] ]
///
/// Special cases:
/// - U(π,0,π) = X (Pauli-X)
/// - U(π/2,0,π) = H (Hadamard)
/// - U(0,0,π/2) = S (Phase gate)
/// - U(0,0,π/4) = T (π/8 gate)
pub struct UGate {
    qubit: Qubit,
    theta: f64,
    phi: f64,
    lambda: f64,
}

impl UGate {
    pub fn new(qubit: Qubit, theta: f64, phi: f64, lambda: f64) -> Self {
        Self {
            qubit,
            theta,
            phi,
            lambda,
        }
    }
}

impl QuantumGate for UGate {
    /// construct the 2x2 unitary matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        array![
            [
                Complex::new(c, 0.0),
                Complex::new(-1.0 * s * self.lambda.cos(), -1.0 * s * self.lambda.sin())
            ],
            [
                Complex::new(s * self.phi.cos(), s * self.phi.sin()),
                Complex::new(
                    c * (self.phi.cos() * self.lambda.cos() - self.phi.sin() * self.lambda.sin()),
                    c * (self.phi.cos() * self.lambda.sin() + self.phi.sin() * self.lambda.cos())
                )
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("U({:.4}, {:.4}, {:.4})", self.theta, self.phi, self.lambda)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target_qubit()]
    }
}

impl SingleQubit for UGate {
    /// returns the index of the qubit this gate operates on.
    fn target_qubit(&self) -> Qubit {
        self.qubit
    }
}

impl SingleQubitGate for UGate {
    /// returns the arbitrary real part of alpha
    fn alpha_re(&self) -> f64 {
        (self.theta / 2.0).cos() * self.global_phase().cos()
    }

    /// returns the arbitrary imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        (self.theta / 2.0).cos() * self.global_phase().sin()
    }

    /// returns the arbitrary real part of beta
    fn beta_re(&self) -> f64 {
        (-1.0) * (self.theta / 2.0).sin() * (self.lambda - self.phi).cos()
    }

    /// returns the arbitrary imaginary part of beta
    fn beta_im(&self) -> f64 {
        (-1.0) * (self.theta / 2.0).sin() * (self.lambda - self.phi).sin()
    }

    /// returns the arbitrary global phase
    fn global_phase(&self) -> f64 {
        (self.phi + self.lambda) / 2.0
    }
}
