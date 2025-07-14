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

use crate::operations::QuantumGate;
use crate::operations::single_qubit::{SingleQubitGate, SingleQubitType};
use crate::{Complex, GateType, Matrix, Qubit};

#[derive(Debug, Clone)]
/// Represents a phase shift operation that only affects the |0⟩ state of a qubit.
///
/// This gate applies a phase shift of theta to the |0⟩ state while leaving the |1⟩ state unchanged.
///
/// The matrix form is:
///
/// PhaseShiftState0(θ) = [ [ cos(θ) + i*sin(θ), 0 ],
///                         [ 0                , 1 ] ]
///
///                     = [ [ e^(iθ), 0 ],
///                         [ 0     , 1 ] ]
///
/// This is the complement to PhaseShiftState1 and allows for precise control of qubit phases.
pub struct U0Gate {
    target: Qubit,
    theta: f64,
}

impl U0Gate {
    pub fn new(target: Qubit, theta: f64) -> Self {
        Self { target, theta }
    }
}

impl QuantumGate for U0Gate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = self.theta;
        array![
            [
                Complex::new(theta.cos(), theta.sin()),
                Complex::new(0.0, 0.0)
            ],
            [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)]
        ]
    }

    fn name(&self) -> String {
        format!(
            "Phase-Shift-0(target={}, theta={:.4})",
            self.target, self.theta
        )
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::U0Gate(Self::new(self.target, self.theta)))
    }
}

impl SingleQubitGate for U0Gate {
    fn alpha_re(&self) -> f64 {
        (self.theta / 2.0).cos()
    }

    fn alpha_im(&self) -> f64 {
        (self.theta / 2.0).sin()
    }

    fn beta_re(&self) -> f64 {
        0.0
    }

    fn beta_im(&self) -> f64 {
        0.0
    }

    fn global_phase(&self) -> f64 {
        self.theta / 2.0
    }
}

#[derive(Debug, Clone)]
/// Represents a phase shift operation that only affects the |1⟩ state of a qubit.
///
/// This gate applies a phase shift of theta to the |1⟩ state while leaving the |0⟩ state unchanged.
///
/// The matrix form is:
///
/// PhaseShiftState1(θ) = [ [ 1, 0                 ],
///                         [ 0, cos(θ) + i*sin(θ) ] ]
///
///                     = [ [ 1, 0      ],
///                         [ 0, e^(iθ) ] ]
///
/// This is useful for adjusting relative phases between the computational basis states.
pub struct U1Gate {
    target: Qubit,
    theta: f64,
}

impl U1Gate {
    pub fn new(target: Qubit, theta: f64) -> Self {
        Self { target, theta }
    }
}

impl QuantumGate for U1Gate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = self.theta;
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [
                Complex::new(0.0, 0.0),
                Complex::new(theta.cos(), theta.sin())
            ]
        ]
    }

    fn name(&self) -> String {
        format!(
            "Phase-Shift-1(target={}, theta={:.4})",
            self.target, self.theta
        )
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::U1Gate(Self::new(self.target, self.theta)))
    }
}

impl SingleQubitGate for U1Gate {
    fn alpha_re(&self) -> f64 {
        (self.theta / 2.0).cos()
    }

    fn alpha_im(&self) -> f64 {
        (-1.0) * (self.theta / 2.0).sin()
    }

    fn beta_re(&self) -> f64 {
        0.0
    }

    fn beta_im(&self) -> f64 {
        0.0
    }

    fn global_phase(&self) -> f64 {
        self.theta / 2.0
    }
}

#[derive(Debug, Clone)]
/// Represents the U2 gate, a single-qubit rotation gate.
///
/// U2(φ, λ) = (1/√2) * [ [      1,    -e^(iλ) ],
///                       [ e^(iφ), e^(i(φ+λ)) ] ]
///
/// This gate performs a rotation equivalent to RZ(φ) RY(π/2) RZ(λ).
pub struct U2Gate {
    target: Qubit,
    phi: f64,    // φ parameter
    lambda: f64, // λ parameter
}

impl U2Gate {
    pub fn new(target: Qubit, phi: f64, lambda: f64) -> Self {
        Self {
            target,
            phi,
            lambda,
        }
    }
}

impl QuantumGate for U2Gate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let inv_sqrt2 = 1.0 / (2.0_f64).sqrt();
        let exp_i_phi = Complex::new(0.0, self.phi).exp();
        let exp_i_lambda = Complex::new(0.0, self.lambda).exp();
        let exp_i_phi_lambda = Complex::new(0.0, self.phi + self.lambda).exp();

        array![
            [Complex::new(inv_sqrt2, 0.0), -exp_i_lambda * inv_sqrt2],
            [exp_i_phi * inv_sqrt2, exp_i_phi_lambda * inv_sqrt2]
        ]
    }

    fn name(&self) -> String {
        format!(
            "U2(target={}, phi={:.3}, lambda={:.3})",
            self.target, self.phi, self.lambda
        )
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::U2Gate(Self::new(
            self.target,
            self.phi,
            self.lambda,
        )))
    }
}

impl SingleQubitGate for U2Gate {
    fn alpha_re(&self) -> f64 {
        1.0 / (2.0_f64).sqrt()
    }

    fn alpha_im(&self) -> f64 {
        0.0
    }

    fn beta_re(&self) -> f64 {
        self.phi.cos() / (2.0_f64).sqrt()
    }

    fn beta_im(&self) -> f64 {
        self.phi.sin() / (2.0_f64).sqrt()
    }

    fn global_phase(&self) -> f64 {
        (self.phi + self.lambda) / 2.0
    }
}

#[derive(Debug, Clone)]
/// Represents the U3 gate, the most general single-qubit unitary gate.
///
/// U3(θ, φ, λ) = [ [ cos(θ/2),           -e^(iλ)sin(θ/2)    ],
///                [ e^(iφ)sin(θ/2),     e^(i(φ+λ))cos(θ/2) ] ]
///
/// This gate can represent any single-qubit rotation.
pub struct U3Gate {
    target: Qubit,
    theta: f64,  // θ parameter
    phi: f64,    // φ parameter
    lambda: f64, // λ parameter
}

impl U3Gate {
    pub fn new(target: Qubit, theta: f64, phi: f64, lambda: f64) -> Self {
        Self {
            target,
            theta,
            phi,
            lambda,
        }
    }
}

impl QuantumGate for U3Gate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let cos_half_theta = (self.theta / 2.0).cos();
        let sin_half_theta = (self.theta / 2.0).sin();
        let exp_i_phi = Complex::new(0.0, self.phi).exp();
        let exp_i_lambda = Complex::new(0.0, self.lambda).exp();
        let exp_i_phi_lambda = Complex::new(0.0, self.phi + self.lambda).exp();

        array![
            [
                Complex::new(cos_half_theta, 0.0),
                -exp_i_lambda * sin_half_theta
            ],
            [
                exp_i_phi * sin_half_theta,
                exp_i_phi_lambda * cos_half_theta
            ]
        ]
    }

    fn name(&self) -> String {
        format!(
            "U3(target={}, θ={:.3}, φ={:.3}, λ={:.3})",
            self.target, self.theta, self.phi, self.lambda
        )
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::U3Gate(Self::new(
            self.target,
            self.theta,
            self.phi,
            self.lambda,
        )))
    }
}

impl SingleQubitGate for U3Gate {
    fn alpha_re(&self) -> f64 {
        (self.theta / 2.0).cos()
    }

    fn alpha_im(&self) -> f64 {
        0.0
    }

    fn beta_re(&self) -> f64 {
        self.phi.cos() * (self.theta / 2.0).sin()
    }

    fn beta_im(&self) -> f64 {
        self.phi.sin() * (self.theta / 2.0).sin()
    }

    fn global_phase(&self) -> f64 {
        (self.phi + self.lambda) / 2.0
    }
}
