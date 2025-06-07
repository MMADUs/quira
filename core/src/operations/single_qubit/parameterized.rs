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
use crate::constant::PI;
use ndarray::array;

#[derive(Debug, Clone)]
/// Represents a rotation around the X-axis of the Bloch sphere.
///
/// This gate rotates the qubit state around the X-axis by an angle theta.
///
/// The matrix form is:
///
/// RX(θ) = [ [ cos(θ/2)   , -i*sin(θ/2) ],
///           [ -i*sin(θ/2), cos(θ/2)    ] ]
///
/// This gate is equivalent to e^(-i*θ*X/2) where X is the Pauli-X matrix.
pub struct RotateX {
    target: Qubit,
    theta: f64,
}

impl RotateX {
    pub fn new(target: Qubit, theta: f64) -> Self {
        Self { target, theta }
    }
}

impl QuantumGate for RotateX {
    /// construct the 2x2 unitary matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(0.0, -1.0 * s)],
            [Complex::new(0.0, -1.0 * s), Complex::new(c, 0.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("RX(target={}, theta={:.4})", self.target, self.theta)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for RotateX {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        (self.theta / 2.0).cos()
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
        (-1.0) * (self.theta / 2.0).sin()
    }

    /// returns the global phase
    fn global_phase(&self) -> f64 {
        0.0
    }
}

#[derive(Debug, Clone)]
/// Represents a rotation around the Y-axis of the Bloch sphere.
///
/// This gate rotates the qubit state around the Y-axis by an angle theta.
///
/// The matrix form is:
///
/// RY(θ) = [ [ cos(θ/2), -sin(θ/2) ],
///           [ sin(θ/2), cos(θ/2)  ] ]
///
/// This gate is equivalent to e^(-i*θ*Y/2) where Y is the Pauli-Y matrix.
pub struct RotateY {
    target: Qubit,
    theta: f64,
}

impl RotateY {
    pub fn new(target: Qubit, theta: f64) -> Self {
        Self { target, theta }
    }
}

impl QuantumGate for RotateY {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(-1.0 * s, 0.0)],
            [Complex::new(s, 0.0), Complex::new(c, 0.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("RY(target={}, theta={:.4})", self.target, self.theta)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for RotateY {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        (self.theta / 2.0).cos()
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        0.0
    }

    /// returns the real part of beta
    fn beta_re(&self) -> f64 {
        (self.theta / 2.0).sin()
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

#[derive(Debug, Clone)]
/// Represents a rotation around the Z-axis of the Bloch sphere.
///
/// This gate rotates the qubit state around the Z-axis by an angle theta.
///
/// The matrix form is:
///
/// RZ(θ) = [ [ e^(-iθ/2), 0        ],
///           [ 0        , e^(iθ/2) ] ]
///
/// This gate is equivalent to e^(-i*θ*Z/2) where Z is the Pauli-Z matrix.
pub struct RotateZ {
    target: Qubit,
    theta: f64,
}

impl RotateZ {
    pub fn new(target: Qubit, theta: f64) -> Self {
        RotateZ { target, theta }
    }
}

impl QuantumGate for RotateZ {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        array![
            [Complex::new(c, -1.0 * s), Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), Complex::new(c, s)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("RZ(target={}, theta={:.4})", self.target, self.theta)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for RotateZ {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        (self.theta / 2.0).cos()
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        (-1.0) * (self.theta / 2.0).sin()
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

#[derive(Debug, Clone)]
/// Represents a rotation in the XY-plane of the Bloch sphere.
///
/// This gate rotates the qubit state by angle theta in the direction specified by angle phi in the XY-plane.
///
/// The matrix form is:
///
/// RXY(θ, φ) = [ [ cos(θ/2)                   , -sin(θ/2)(sin(φ) + i*cos(φ)) ],
///               [ sin(θ/2)(sin(φ) - i*cos(φ)), cos(θ/2)                     ] ]
///
/// This gate combines aspects of both RX and RY rotations and allows for arbitrary rotations in the XY-plane.
pub struct RotateXY {
    target: Qubit,
    theta: f64,
    phi: f64,
}

impl RotateXY {
    pub fn new(target: Qubit, theta: f64, phi: f64) -> Self {
        Self { target, theta, phi }
    }
}

impl QuantumGate for RotateXY {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        let vx: f64 = (self.phi).cos();
        let vy: f64 = (self.phi).sin();
        array![
            [
                Complex::new(c, 0.0),
                Complex::new(-1.0 * s * vy, -1.0 * s * vx)
            ],
            [Complex::new(s * vy, -1.0 * s * vx), Complex::new(c, 0.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "RXY(target={}, theta={:.4}, phi={:.4})",
            self.target, self.theta, self.phi
        )
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for RotateXY {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        (self.theta / 2.0).cos()
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        0.0
    }

    /// returns the real part of beta
    fn beta_re(&self) -> f64 {
        let s = (self.theta / 2.0).sin();
        let vy = (self.phi).sin();
        s * vy
    }

    /// returns the imaginary part of beta
    fn beta_im(&self) -> f64 {
        let s = (self.theta / 2.0).sin();
        let vx = (self.phi).cos();
        (-1.0) * s * vx
    }

    /// returns the global phase
    fn global_phase(&self) -> f64 {
        0.0
    }
}

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
    /// construct the 2x2 matrix representing the gate.
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

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "Phase-Shift-0(target={}, theta={:.4})",
            self.target, self.theta
        )
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for U0Gate {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        (self.theta / 2.0).cos()
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        (self.theta / 2.0).sin()
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
    /// construct the 2x2 matrix representing the gate.
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

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "Phase-Shift-1(target={}, theta={:.4})",
            self.target, self.theta
        )
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for U1Gate {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        (self.theta / 2.0).cos()
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        (-1.0) * (self.theta / 2.0).sin()
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
        Self { target, phi, lambda }
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
        format!("U2(target={}, phi={:.3}, lambda={:.3})", self.target, self.phi, self.lambda)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
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
        Self { target, theta, phi, lambda }
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
            [Complex::new(cos_half_theta, 0.0), -exp_i_lambda * sin_half_theta],
            [exp_i_phi * sin_half_theta, exp_i_phi_lambda * cos_half_theta]
        ]
    }

    fn name(&self) -> String {
        format!("U3(target={}, θ={:.3}, φ={:.3}, λ={:.3})", 
                self.target, self.theta, self.phi, self.lambda)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
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

#[derive(Debug, Clone)]
/// Represents the X^t gate (X to the power of t).
///
/// X^t = [ [ cos(πt/2) - i*sin(πt/2), -i*sin(πt/2) ],
///         [ -i*sin(πt/2),           cos(πt/2) + i*sin(πt/2) ] ]
///
/// When t=1, this reduces to the standard Pauli-X gate.
pub struct XPowGate {
    target: Qubit,
    exponent: f64, // t parameter
}

impl XPowGate {
    pub fn new(target: Qubit, exponent: f64) -> Self {
        Self { target, exponent }
    }
}

impl QuantumGate for XPowGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let angle = PI * self.exponent / 2.0;
        let cos_angle = angle.cos();
        let sin_angle = angle.sin();
        
        array![
            [Complex::new(cos_angle, -sin_angle), Complex::new(0.0, -sin_angle)],
            [Complex::new(0.0, -sin_angle), Complex::new(cos_angle, sin_angle)]
        ]
    }

    fn name(&self) -> String {
        format!("X^{:.3}(target={})", self.exponent, self.target)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for XPowGate {
    fn alpha_re(&self) -> f64 {
        (PI * self.exponent / 2.0).cos()
    }

    fn alpha_im(&self) -> f64 {
        -(PI * self.exponent / 2.0).sin()
    }

    fn beta_re(&self) -> f64 {
        0.0
    }

    fn beta_im(&self) -> f64 {
        -(PI * self.exponent / 2.0).sin()
    }

    fn global_phase(&self) -> f64 {
        0.0
    }
}

#[derive(Debug, Clone)]
/// Represents the Y^t gate (Y to the power of t).
///
/// Y^t = [ [ cos(πt/2) - i*sin(πt/2), -sin(πt/2) - i*cos(πt/2) ],
///         [ sin(πt/2) - i*cos(πt/2),  cos(πt/2) + i*sin(πt/2) ] ]
///
/// When t=1, this reduces to the standard Pauli-Y gate.
pub struct YPowGate {
    target: Qubit,
    exponent: f64, // t parameter
}

impl YPowGate {
    pub fn new(target: Qubit, exponent: f64) -> Self {
        Self { target, exponent }
    }
}

impl QuantumGate for YPowGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let angle = PI * self.exponent / 2.0;
        let cos_angle = angle.cos();
        let sin_angle = angle.sin();
        
        array![
            [Complex::new(cos_angle, -sin_angle), Complex::new(-sin_angle, -cos_angle)],
            [Complex::new(sin_angle, -cos_angle), Complex::new(cos_angle, sin_angle)]
        ]
    }

    fn name(&self) -> String {
        format!("Y^{:.3}(target={})", self.exponent, self.target)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for YPowGate {
    fn alpha_re(&self) -> f64 {
        (PI * self.exponent / 2.0).cos()
    }

    fn alpha_im(&self) -> f64 {
        -(PI * self.exponent / 2.0).sin()
    }

    fn beta_re(&self) -> f64 {
        -(PI * self.exponent / 2.0).sin()
    }

    fn beta_im(&self) -> f64 {
        -(PI * self.exponent / 2.0).cos()
    }

    fn global_phase(&self) -> f64 {
        0.0
    }
}

#[derive(Debug, Clone)]
/// Represents the Z^t gate (Z to the power of t).
///
/// Z^t = [ [ e^(-iπt/2), 0         ],
///         [ 0,          e^(iπt/2) ] ]
///
/// When t=1, this reduces to the standard Pauli-Z gate.
pub struct ZPowGate {
    target: Qubit,
    exponent: f64, // t parameter
}

impl ZPowGate {
    pub fn new(target: Qubit, exponent: f64) -> Self {
        Self { target, exponent }
    }
}

impl QuantumGate for ZPowGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let angle = PI * self.exponent / 2.0;
        let exp_neg_i_angle = Complex::new(0.0, -angle).exp();
        let exp_pos_i_angle = Complex::new(0.0, angle).exp();
        
        array![
            [exp_neg_i_angle, Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), exp_pos_i_angle]
        ]
    }

    fn name(&self) -> String {
        format!("Z^{:.3}(target={})", self.exponent, self.target)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for ZPowGate {
    fn alpha_re(&self) -> f64 {
        (PI * self.exponent / 2.0).cos()
    }

    fn alpha_im(&self) -> f64 {
        -(PI * self.exponent / 2.0).sin()
    }

    fn beta_re(&self) -> f64 {
        0.0
    }

    fn beta_im(&self) -> f64 {
        0.0
    }

    fn global_phase(&self) -> f64 {
        0.0
    }
}

#[derive(Debug, Clone)]
/// Represents the PhasedXPow gate.
///
/// This gate applies an X^t rotation with an additional phase offset.
/// PhasedXPow(t, φ) = Z^(-φ/π) * X^t * Z^(φ/π)
///
/// The matrix form is:
/// [ [ cos(πt/2) - i*sin(πt/2)*cos(φ), -i*sin(πt/2)*sin(φ) - sin(πt/2)*cos(φ) ],
///   [ -i*sin(πt/2)*sin(φ) + sin(πt/2)*cos(φ), cos(πt/2) + i*sin(πt/2)*cos(φ) ] ]
pub struct PhasedXPowGate {
    target: Qubit,
    exponent: f64, // t parameter
    phase_exponent: f64, // φ parameter
}

impl PhasedXPowGate {
    pub fn new(target: Qubit, exponent: f64, phase_exponent: f64) -> Self {
        Self { target, exponent, phase_exponent }
    }
}

impl QuantumGate for PhasedXPowGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let t_angle = PI * self.exponent / 2.0;
        let cos_t = t_angle.cos();
        let sin_t = t_angle.sin();
        let cos_phi = self.phase_exponent.cos();
        let sin_phi = self.phase_exponent.sin();
        
        array![
            [
                Complex::new(cos_t, -sin_t * cos_phi),
                Complex::new(-sin_t * cos_phi, -sin_t * sin_phi)
            ],
            [
                Complex::new(sin_t * cos_phi, -sin_t * sin_phi),
                Complex::new(cos_t, sin_t * cos_phi)
            ]
        ]
    }

    fn name(&self) -> String {
        format!("PhasedX^{:.3}(target={}, φ={:.3})", 
                self.exponent, self.target, self.phase_exponent)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for PhasedXPowGate {
    fn alpha_re(&self) -> f64 {
        (PI * self.exponent / 2.0).cos()
    }

    fn alpha_im(&self) -> f64 {
        -(PI * self.exponent / 2.0).sin() * self.phase_exponent.cos()
    }

    fn beta_re(&self) -> f64 {
        (PI * self.exponent / 2.0).sin() * self.phase_exponent.cos()
    }

    fn beta_im(&self) -> f64 {
        -(PI * self.exponent / 2.0).sin() * self.phase_exponent.sin()
    }

    fn global_phase(&self) -> f64 {
        0.0
    }
}

#[derive(Debug, Clone)]
/// Represents the H^t gate (Hadamard to the power of t).
///
/// H^t gate interpolates between the identity (t=0) and Hadamard (t=1) gates.
/// 
/// H^t = [ [ cos(πt/4) + sin(πt/4), cos(πt/4) - sin(πt/4) ],
///         [ cos(πt/4) - sin(πt/4), cos(πt/4) + sin(πt/4) ] ] / √2
///
/// When t=1, this reduces to the standard Hadamard gate.
pub struct HPowGate {
    target: Qubit,
    exponent: f64, // t parameter
}

impl HPowGate {
    pub fn new(target: Qubit, exponent: f64) -> Self {
        Self { target, exponent }
    }
}

impl QuantumGate for HPowGate {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let angle = PI * self.exponent / 4.0;
        let cos_angle = angle.cos();
        let sin_angle = angle.sin();
        let inv_sqrt2 = 1.0 / (2.0_f64).sqrt();
        
        let a = (cos_angle + sin_angle) * inv_sqrt2;
        let b = (cos_angle - sin_angle) * inv_sqrt2;
        
        array![
            [Complex::new(a, 0.0), Complex::new(b, 0.0)],
            [Complex::new(b, 0.0), Complex::new(a, 0.0)]
        ]
    }

    fn name(&self) -> String {
        format!("H^{:.3}(target={})", self.exponent, self.target)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for HPowGate {
    fn alpha_re(&self) -> f64 {
        let angle = PI * self.exponent / 4.0;
        (angle.cos() + angle.sin()) / (2.0_f64).sqrt()
    }

    fn alpha_im(&self) -> f64 {
        0.0
    }

    fn beta_re(&self) -> f64 {
        let angle = PI * self.exponent / 4.0;
        (angle.cos() - angle.sin()) / (2.0_f64).sqrt()
    }

    fn beta_im(&self) -> f64 {
        0.0
    }

    fn global_phase(&self) -> f64 {
        0.0
    }
}

#[derive(Debug, Clone)]
/// Represents a rotation around an arbitrary axis defined by spherical coordinates on the Bloch sphere.
///
/// This gate rotates the qubit state by angle theta around an axis defined by the spherical coordinates
/// (spherical_theta, spherical_phi).
///
/// The matrix form is:
///
/// RotateAroundSphericalAxis(θ, θ_s, φ_s) =
///
/// [ [ cos(θ/2) - i*sin(θ/2)*cos(θ_s)                      , -i*sin(θ/2)*(sin(θ_s)*sin(φ_s) + i*sin(θ_s)*cos(φ_s)) ],
///   [ i*sin(θ/2)*(sin(θ_s)*sin(φ_s) - i*sin(θ_s)*cos(φ_s)), cos(θ/2) + i*sin(θ/2)*cos(θ_s)                        ] ]
///
/// where θ is the rotation angle, θ_s is the polar angle, and φ_s is the azimuthal angle
/// in spherical coordinates defining the rotation axis.
///
/// This is a generalized rotation gate that can represent any single-qubit unitary operation.
pub struct RotateAroundSphericalAxis {
    target: Qubit,
    theta: f64,
    spherical_theta: f64,
    spherical_phi: f64,
}

impl RotateAroundSphericalAxis {
    pub fn new(target: Qubit, theta: f64, spherical_theta: f64, spherical_phi: f64) -> Self {
        Self {
            target,
            theta,
            spherical_theta,
            spherical_phi,
        }
    }
}

impl QuantumGate for RotateAroundSphericalAxis {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        let vx: f64 = (self.spherical_theta).sin() * (self.spherical_phi).cos();
        let vy: f64 = (self.spherical_theta).sin() * (self.spherical_phi).sin();
        let vz: f64 = (self.spherical_theta).cos();
        array![
            [
                Complex::new(c, -1.0 * s * vz),
                Complex::new(-1.0 * s * vy, -1.0 * s * vx)
            ],
            [Complex::new(s * vy, -1.0 * s * vx), Complex::new(c, s * vz)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "RAS(target={}, theta={:.4}, sphere_theta={:.4}, sphere_phi{:.4})",
            self.target, self.theta, self.spherical_theta, self.spherical_phi
        )
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for RotateAroundSphericalAxis {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        (self.theta / 2.0).cos()
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        let s = (self.theta / 2.0).sin();
        let vz = (self.spherical_theta).cos();
        (-1.0) * s * vz
    }

    /// returns the real part of beta
    fn beta_re(&self) -> f64 {
        let s = (self.theta / 2.0).sin();
        let vy = (self.spherical_phi).sin();
        let st = (self.spherical_theta).sin();
        s * vy * st
    }

    /// returns the imaginary part of beta
    fn beta_im(&self) -> f64 {
        let s = (self.theta / 2.0).sin();
        let vx = (self.spherical_phi).cos();
        let st = (self.spherical_theta).sin();
        (-1.0) * s * vx * st
    }

    /// returns the global phase
    fn global_phase(&self) -> f64 {
        0.0
    }
}
