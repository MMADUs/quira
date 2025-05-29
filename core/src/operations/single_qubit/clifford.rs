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
use crate::constant::PI;
use crate::operations::QuantumGate;
use crate::types::{Complex, Matrix, Qubit};
use ndarray::array;

#[derive(Debug, Clone)]
/// Represents the Pauli-X gate (also known as the bit-flip or NOT gate).
///
/// This gate flips the state of a qubit from |0⟩ to |1⟩ and vice versa.
/// It is the quantum equivalent of the classical NOT gate.
///
/// The matrix form is:
///
/// X = [ [ 0, 1 ],
///       [ 1, 0 ] ]
pub struct PauliX {
    target: Qubit,
}

impl PauliX {
    pub fn new(target: Qubit) -> Self {
        Self { target }
    }
}

impl QuantumGate for PauliX {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(0.0, 0.0), Complex::new(1.0, 0.0)],
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("X(target={})", self.target)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for PauliX {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        0.0
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
        -1.0
    }

    /// returns the global phase
    fn global_phase(&self) -> f64 {
        PI / 2.0
    }
}

#[derive(Debug, Clone)]
/// Represents the Pauli-Y gate.
///
/// This gate performs a rotation around the Y-axis by π radians,
/// with an additional phase factor.
///
/// The matrix form is:
///
/// Y = [ [ 0, -i ],
///       [ i,  0 ] ]
pub struct PauliY {
    target: Qubit,
}

impl PauliY {
    pub fn new(target: Qubit) -> Self {
        Self { target }
    }
}

impl QuantumGate for PauliY {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(0.0, 0.0), Complex::new(0.0, -1.0)],
            [Complex::new(0.0, 1.0), Complex::new(0.0, 0.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("Y(target={})", self.target)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for PauliY {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        0.0
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        0.0
    }

    /// returns the real part of beta
    fn beta_re(&self) -> f64 {
        1.0
    }

    /// returns the imaginary part of beta
    fn beta_im(&self) -> f64 {
        0.0
    }

    /// returns the global phase
    fn global_phase(&self) -> f64 {
        PI / 2.0
    }
}

#[derive(Debug, Clone)]
/// Represents the Pauli-Z gate (also known as the phase-flip gate).
///
/// This gate leaves the |0⟩ state unchanged and flips the sign of the |1⟩ state.
/// It represents a rotation around the Z-axis by π radians.
///
/// The matrix form is:
///
/// Z = [ [ 1,  0 ],
///       [ 0, -1 ] ]
pub struct PauliZ {
    target: Qubit,
}

impl PauliZ {
    pub fn new(target: Qubit) -> Self {
        Self { target }
    }
}

impl QuantumGate for PauliZ {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), Complex::new(-1.0, 0.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("Z(target={})", self.target)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for PauliZ {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        0.0
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        -1.0
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
        PI / 2.0
    }
}

#[derive(Debug, Clone)]
/// Represents the square root of the Pauli-X gate.
///
/// This gate performs half the action of the Pauli-X gate, equivalent to
/// a 90-degree rotation around the X-axis.
///
/// The matrix form is:
///
/// √X = [ [ (1+i)/2, (1-i)/2 ],
///        [ (1-i)/2, (1+i)/2 ] ]
///
/// Applying this gate twice is equivalent to applying the X gate once.
pub struct SqrtPauliX {
    target: Qubit,
}

impl SqrtPauliX {
    pub fn new(target: Qubit) -> Self {
        Self { target }
    }
}

impl QuantumGate for SqrtPauliX {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = PI / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(0.0, -1.0 * s)],
            [Complex::new(0.0, -1.0 * s), Complex::new(c, 0.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("Sqrt-X(target={})", self.target)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for SqrtPauliX {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        (PI / 4.0).cos()
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
        (PI / 4.0).sin() * (-1.0)
    }

    /// returns the global phase
    fn global_phase(&self) -> f64 {
        0.0
    }
}

#[derive(Debug, Clone)]
/// Represents the inverse of the square root of the Pauli-X gate.
///
/// This gate is the adjoint (inverse) of the √X gate, equivalent to
/// a -90-degree rotation around the X-axis.
///
/// The matrix form is the complex conjugate transpose of the √X matrix.
///
/// Applying this gate twice is equivalent to applying the X gate once,
/// and applying it after the √X gate results in the identity gate.
pub struct InvSqrtPauliX {
    target: Qubit,
}

impl InvSqrtPauliX {
    pub fn new(target: Qubit) -> Self {
        Self { target }
    }
}

impl QuantumGate for InvSqrtPauliX {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = PI / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(0.0, 1.0 * s)],
            [Complex::new(0.0, 1.0 * s), Complex::new(c, 0.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("Inv-Sqrt-X(target={})", self.target)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for InvSqrtPauliX {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        (PI / 4.0).cos()
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
        (PI / 4.0).sin()
    }

    /// returns the global phase
    fn global_phase(&self) -> f64 {
        0.0
    }
}

#[derive(Debug, Clone)]
/// Represents the square root of the Pauli-Y gate.
///
/// This gate performs half the action of the Pauli-Y gate, equivalent to
/// a 90-degree rotation around the Y-axis.
///
/// Applying this gate twice is equivalent to applying the Y gate once.
pub struct SqrtPauliY {
    target: Qubit,
}

impl SqrtPauliY {
    pub fn new(target: Qubit) -> Self {
        Self { target }
    }
}

impl QuantumGate for SqrtPauliY {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = PI / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(-1.0 * s, 0.0)],
            [Complex::new(s, 0.0), Complex::new(c, 0.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("Sqrt-Y(target={})", self.target)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for SqrtPauliY {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        (PI / 4.0).cos()
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        0.0
    }

    /// returns the real part of beta
    fn beta_re(&self) -> f64 {
        (PI / 4.0).sin()
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
/// Represents the inverse of the square root of the Pauli-Y gate.
///
/// This gate is the adjoint (inverse) of the √Y gate, equivalent to
/// a -90-degree rotation around the Y-axis.
///
/// The matrix form is the complex conjugate transpose of the √Y matrix.
///
/// Applying this gate twice is equivalent to applying the Y gate once,
/// and applying it after the √Y gate results in the identity gate.
pub struct InvSqrtPauliY {
    target: Qubit,
}

impl InvSqrtPauliY {
    pub fn new(target: Qubit) -> Self {
        Self { target }
    }
}

impl QuantumGate for InvSqrtPauliY {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = (-1.0 * PI) / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        array![
            [Complex::new(c, 0.0), Complex::new(-1.0 * s, 0.0)],
            [Complex::new(s, 0.0), Complex::new(c, 0.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("Inv-Sqrt-Y(target={})", self.target)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for InvSqrtPauliY {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        (PI / 4.0).cos()
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        0.0
    }

    /// returns the real part of beta
    fn beta_re(&self) -> f64 {
        (-1.0) * (PI / 4.0).sin()
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
/// Represents the Hadamard gate, a fundamental quantum gate.
///
/// This gate creates a superposition of the |0⟩ and |1⟩ states.
/// It maps |0⟩ to (|0⟩ + |1⟩)/√2 and |1⟩ to (|0⟩ - |1⟩)/√2.
///
/// The matrix form is:
///
/// H = (1/√2) * [ [ 1,  1 ],
///                [ 1, -1 ] ]
///
/// The Hadamard gate is self-inverse: applying it twice results in the identity.
pub struct Hadamard {
    target: Qubit,
}

impl Hadamard {
    pub fn new(target: Qubit) -> Self {
        Self { target }
    }
}

impl QuantumGate for Hadamard {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let f = 1.0 / ((2.0_f64).sqrt());
        array![
            [Complex::new(f, 0.0), Complex::new(f, 0.0)],
            [Complex::new(f, 0.0), Complex::new(-1.0 * f, 0.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("H(target={})", self.target)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for Hadamard {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        0.0
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        -1.0 / ((2.0_f64).sqrt())
    }

    /// returns the real part of beta
    fn beta_re(&self) -> f64 {
        0.0
    }

    /// returns the imaginary part of beta
    fn beta_im(&self) -> f64 {
        -1.0 / ((2.0_f64).sqrt())
    }

    /// returns the global phase
    fn global_phase(&self) -> f64 {
        PI / 2.0
    }
}

#[derive(Debug, Clone)]
/// Represents the S gate (also known as the phase gate or Z90 gate).
///
/// This gate leaves the |0⟩ state unchanged and maps |1⟩ to i|1⟩,
/// introducing a π/2 phase shift to the |1⟩ state.
///
/// The matrix form is:
///
/// S = [ [ 1, 0 ],
///       [ 0, i ] ]
///
/// Applying this gate twice is equivalent to applying the Z gate once.
pub struct SGate {
    target: Qubit,
}

impl SGate {
    pub fn new(target: Qubit) -> Self {
        Self { target }
    }
}

impl QuantumGate for SGate {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), Complex::new(0.0, 1.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("S(target={})", self.target)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for SGate {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        1.0 / ((2.0_f64).sqrt())
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        -1.0 / ((2.0_f64).sqrt())
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
        PI / 4.0
    }
}

#[derive(Debug, Clone)]
/// Represents the inverse of the S gate.
///
/// This gate leaves the |0⟩ state unchanged and maps |1⟩ to -i|1⟩,
/// introducing a -π/2 phase shift to the |1⟩ state.
///
/// The matrix form is:
///
/// S† = [ [ 1,  0 ],
///        [ 0, -i ] ]
///
/// Applying this gate after the S gate results in the identity gate.
pub struct InvSGate {
    target: Qubit,
}

impl InvSGate {
    pub fn new(target: Qubit) -> Self {
        Self { target }
    }
}

impl QuantumGate for InvSGate {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [Complex::new(1.0, 0.0), Complex::new(0.0, 0.0)],
            [Complex::new(0.0, 0.0), Complex::new(0.0, -1.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("Inv-S(target={})", self.target)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for InvSGate {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        1.0 / ((2.0_f64).sqrt())
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        1.0 / ((2.0_f64).sqrt())
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
        (-1.0 * PI) / 4.0
    }
}

#[derive(Debug, Clone)]
/// Represents the √X (square root of X) gate, also known as the √NOT gate.
///
/// This gate performs half of a NOT operation. Applying it twice is equivalent to a full X (NOT) gate.
///
/// The matrix form is:
///
/// SX = [ [ (1+i)/2, (1-i)/2 ],
///        [ (1-i)/2, (1+i)/2 ] ]
///
/// This gate is equivalent to e^(i*π/4) * RX(π/2) where RX is the rotation around X-axis.
/// It maps |0⟩ to (|0⟩ + i|1⟩)/√2 and |1⟩ to (i|0⟩ + |1⟩)/√2.
pub struct SXGate {
    target: Qubit,
}

impl SXGate {
    pub fn new(target: Qubit) -> Self {
        Self { target }
    }
}

impl QuantumGate for SXGate {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = PI / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        let gp = Complex::new(0.0, PI / 4.0).exp();
        array![
            [gp * Complex::new(c, 0.0), gp * Complex::new(0.0, -1.0 * s)],
            [gp * Complex::new(0.0, -1.0 * s), gp * Complex::new(c, 0.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("SX(target={})", self.target)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for SXGate {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        (PI / 4.0).cos()
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
        (-1.0) * (PI / 4.0).sin()
    }

    /// returns the global phase
    fn global_phase(&self) -> f64 {
        PI / 4.0
    }
}

#[derive(Debug, Clone)]
/// Represents the inverse of the √X gate (SX†).
///
/// This gate is the Hermitian conjugate of the SX gate. Applying SX followed by SX† returns to the original state.
///
/// The matrix form is:
///
/// SX† = [ [ (1-i)/2, (1-i)/2 ],
///         [ (1-i)/2, (1+i)/2 ] ]
///
/// This gate is equivalent to e^(i*π/4) * RX(-π/2) where RX is the rotation around X-axis.
pub struct InvSXGate {
    target: Qubit,
}

impl InvSXGate {
    pub fn new(target: Qubit) -> Self {
        Self { target }
    }
}

impl QuantumGate for InvSXGate {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let theta: f64 = PI / 2.0;
        let c: f64 = (theta / 2.0).cos();
        let s: f64 = (theta / 2.0).sin();
        let gp = Complex::new(0.0, PI / 4.0).exp();
        array![
            [gp * Complex::new(c, 0.0), gp * Complex::new(0.0, 1.0 * s)],
            [gp * Complex::new(0.0, 1.0 * s), gp * Complex::new(c, 0.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!("Inv-SX(target={})", self.target)
    }

    /// construct targets for quantum state
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }
}

impl SingleQubitGate for InvSXGate {
    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        (PI / 4.0).cos()
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
        (PI / 4.0).sin()
    }

    /// returns the global phase
    fn global_phase(&self) -> f64 {
        PI / 4.0
    }
}
