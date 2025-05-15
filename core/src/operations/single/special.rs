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

use super::SingleQubit;
use crate::math::{Complex, Matrix, PI};
use crate::operations::QuantumGate;
use ndarray::array;

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
    qubit: usize,
}

impl SXGate {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
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
        String::from("SX")
    }
}

impl SingleQubit for SXGate {
    /// returns the index of the qubit this gate operates on.
    fn target_qubit(&self) -> usize {
        self.qubit
    }

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
    qubit: usize,
}

impl InvSXGate {
    pub fn new(qubit: usize) -> Self {
        Self { qubit }
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
        String::from("Inv-SX")
    }
}

impl SingleQubit for InvSXGate {
    /// returns the index of the qubit this gate operates on.
    fn target_qubit(&self) -> usize {
        self.qubit
    }

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

#[derive(Debug, Clone)]
/// Represents a generalized π (Pi) pulse quantum gate.
///
/// This gate flips the state of a qubit around an axis in the XY plane of the Bloch sphere.
/// The angle theta determines the orientation of the rotation axis in the XY plane.
///
/// The matrix form is:
///
/// GPi(θ) = [ [ 0              , cos(θ)-i*sin(θ) ],
///            [ cos(θ)+i*sin(θ), 0               ] ]
///
/// When θ=0, this becomes equivalent to the Pauli-X gate.
pub struct GPi {
    qubit: usize,
    theta: f64,
}

impl GPi {
    pub fn new(qubit: usize, theta: f64) -> Self {
        Self { qubit, theta }
    }
}

impl QuantumGate for GPi {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta).cos();
        let s: f64 = (self.theta).sin();
        array![
            [Complex::new(0.0, 0.0), Complex::new(c, -1.0 * s)],
            [Complex::new(c, s), Complex::new(0.0, 0.0)]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        String::from("GPi")
    }
}

impl SingleQubit for GPi {
    /// returns the index of the qubit this gate operates on.
    fn target_qubit(&self) -> usize {
        self.qubit
    }

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
        (self.theta).sin()
    }

    /// returns the imaginary part of beta
    fn beta_im(&self) -> f64 {
        (self.theta).cos() * -1.0
    }

    /// returns the global phase
    fn global_phase(&self) -> f64 {
        PI / 2.0
    }
}

#[derive(Debug, Clone)]
/// Represents a generalized π/2 (Pi/2) pulse quantum gate.
///
/// This gate performs a 90-degree rotation around an axis in the XY plane of the Bloch sphere.
/// The angle theta determines the orientation of the rotation axis in the XY plane.
///
/// The matrix form is:
///
/// GPi2(θ) = (1/√2) * [ [ 1              , -sin(θ)-i*cos(θ) ],
///                      [ sin(θ)-i*cos(θ), 1                ] ]
///
/// When θ=0, this is equivalent to a 90-degree rotation around the X-axis (√X gate without the global phase).
pub struct GPi2 {
    qubit: usize,
    theta: f64,
}

impl GPi2 {
    pub fn new(qubit: usize, theta: f64) -> Self {
        Self { qubit, theta }
    }
}

impl QuantumGate for GPi2 {
    /// construct the 2x2 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta).cos();
        let s: f64 = (self.theta).sin();
        array![
            [Complex::new(1.0, 0.0), Complex::new(-1.0 * s, -1.0 * c)],
            [Complex::new(s, -1.0 * c), Complex::new(1.0, 0.0)]
        ] / 2.0_f64.sqrt()
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        String::from("GPi2")
    }
}

impl SingleQubit for GPi2 {
    /// returns the index of the qubit this gate operates on.
    fn target_qubit(&self) -> usize {
        self.qubit
    }

    /// returns the real part of alpha
    fn alpha_re(&self) -> f64 {
        1.0 / 2.0_f64.sqrt()
    }

    /// returns the imaginary part of alpha
    fn alpha_im(&self) -> f64 {
        0.0
    }

    /// returns the real part of beta
    fn beta_re(&self) -> f64 {
        (self.theta).sin() / 2.0_f64.sqrt()
    }

    /// returns the imaginary part of beta
    fn beta_im(&self) -> f64 {
        (self.theta).cos() * -1.0 / 2.0_f64.sqrt()
    }

    /// returns the global phase
    fn global_phase(&self) -> f64 {
        0.0
    }
}
