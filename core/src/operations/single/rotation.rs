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
use crate::operations::QuantumGate;
use crate::types::{Complex, Matrix, Qubit};
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
    qubit: Qubit,
    theta: f64,
}

impl RotateX {
    pub fn new(qubit: Qubit, theta: f64) -> Self {
        Self { qubit, theta }
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
        format!("RX({:.4})", self.theta)
    }

    /// returns the index of the qubit this gate operates on.
    fn target_qubit(&self) -> Qubit {
        self.qubit
    }
}

impl SingleQubit for RotateX {
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
    qubit: Qubit,
    theta: f64,
}

impl RotateY {
    pub fn new(qubit: Qubit, theta: f64) -> Self {
        Self { qubit, theta }
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
        format!("RY({:.4})", self.theta)
    }

    /// returns the index of the qubit this gate operates on.
    fn target_qubit(&self) -> Qubit {
        self.qubit
    }
}

impl SingleQubit for RotateY {
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
    qubit: Qubit,
    theta: f64,
}

impl RotateZ {
    pub fn new(qubit: Qubit, theta: f64) -> Self {
        RotateZ { qubit, theta }
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
        format!("RZ({:.4})", self.theta)
    }

    /// returns the index of the qubit this gate operates on.
    fn target_qubit(&self) -> Qubit {
        self.qubit
    }
}

impl SingleQubit for RotateZ {
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
    qubit: Qubit,
    theta: f64,
    phi: f64,
}

impl RotateXY {
    pub fn new(qubit: Qubit, theta: f64, phi: f64) -> Self {
        Self { qubit, theta, phi }
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
        format!("RXY({:.4}, {:.4})", self.theta, self.phi)
    }

    /// returns the index of the qubit this gate operates on.
    fn target_qubit(&self) -> Qubit {
        self.qubit
    }
}

impl SingleQubit for RotateXY {
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
pub struct PhaseShiftState1 {
    qubit: Qubit,
    theta: f64,
}

impl PhaseShiftState1 {
    pub fn new(qubit: Qubit, theta: f64) -> Self {
        Self { qubit, theta }
    }
}

impl QuantumGate for PhaseShiftState1 {
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
        format!("Phase-Shift-1({:.4})", self.theta)
    }

    /// returns the index of the qubit this gate operates on.
    fn target_qubit(&self) -> Qubit {
        self.qubit
    }
}

impl SingleQubit for PhaseShiftState1 {
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
pub struct PhaseShiftState0 {
    qubit: Qubit,
    theta: f64,
}

impl PhaseShiftState0 {
    pub fn new(qubit: Qubit, theta: f64) -> Self {
        Self { qubit, theta }
    }
}

impl QuantumGate for PhaseShiftState0 {
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
        format!("Phase-Shift-0({:.4})", self.theta)
    }

    /// returns the index of the qubit this gate operates on.
    fn target_qubit(&self) -> Qubit {
        self.qubit
    }
}

impl SingleQubit for PhaseShiftState0 {
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
    qubit: Qubit,
    theta: f64,
    spherical_theta: f64,
    spherical_phi: f64,
}

impl RotateAroundSphericalAxis {
    pub fn new(qubit: Qubit, theta: f64, spherical_theta: f64, spherical_phi: f64) -> Self {
        Self {
            qubit,
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
        format!("RAS({:.4}, {:.4}, {:.4})", self.theta, self.spherical_theta, self.spherical_phi)
    }

    /// returns the index of the qubit this gate operates on.
    fn target_qubit(&self) -> Qubit {
        self.qubit
    }
}

impl SingleQubit for RotateAroundSphericalAxis {
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
