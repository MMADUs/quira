/*
Copyright (c) 2024-2025 Quira, Inc.

This file is part of Quira

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

use ndarray::array;

use crate::bit::QuantumBit;
use crate::operations::GateType;
use crate::operations::single_qubit::SingleQubitType;
use crate::ops::QuantumGate;
use crate::types::{Complex, Matrix};

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
    target: QuantumBit,
    theta: f64,
    phi: f64,
    lambda: f64,
}

impl UGate {
    pub fn new(target: &QuantumBit, theta: f64, phi: f64, lambda: f64) -> Self {
        Self {
            target: target.clone(),
            theta,
            phi,
            lambda,
        }
    }
}

impl QuantumGate for UGate {
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

    fn name(&self) -> String {
        format!(
            "U(target={}, theta={:.4}, phi={:.4}, lambda={:.4})",
            self.target, self.theta, self.phi, self.lambda
        )
    }

    fn construct_targets(&self) -> Vec<usize> {
        vec![self.target.index()]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::UGate(Self::new(
            &self.target,
            self.theta,
            self.phi,
            self.lambda,
        )))
    }
}
