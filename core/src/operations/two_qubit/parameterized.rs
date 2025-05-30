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

use crate::{
    constant::PI,
    operations::QuantumGate,
    types::{Complex, Matrix, Qubit},
};

use ndarray::array;

#[derive(Debug, Clone)]
/// Represents a phase-shifted controlled-Z gate.
///
/// This gate applies a phase shift to the control qubit and then performs
/// a controlled-Z operation. The phase shift is applied as e^(iφ) when the
/// control qubit is in the |1⟩ state.
///
/// The matrix form is:
/// PSCZ = [ [ 1,      0,      0,           0 ],
///          [ 0, e^(iφ),      0,           0 ],
///          [ 0,      0, e^(iφ),           0 ],
///          [ 0,      0,      0, e^(i(2φ+π)) ] ]
///
/// # Parameters
/// - `control`: The control qubit index
/// - `target`: The target qubit index  
/// - `phi`: The phase shift parameter φ
pub struct PhaseShiftedControlledZ {
    control: Qubit,
    target: Qubit,
    phi: f64,
}

impl PhaseShiftedControlledZ {
    pub fn new(control: Qubit, target: Qubit, phi: f64) -> Self {
        Self {
            control,
            target,
            phi,
        }
    }
}

impl QuantumGate for PhaseShiftedControlledZ {
    /// construct the 4x4 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let phi: f64 = self.phi;
        let c: f64 = phi.cos();
        let s: f64 = phi.sin();
        let cos: f64 = (2.0 * phi + PI).cos();
        let sin: f64 = (2.0 * phi + PI).sin();
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(c, s),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(c, s),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(cos, sin),
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "PSCZ(control={}, target={}, phi={:.4})",
            self.control, self.target, self.phi
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
/// Represents a phase-shifted controlled phase gate.
///
/// This gate combines a phase shift on intermediate states with a controlled
/// phase operation. It applies phase shifts e^(iφ) on the |01⟩ and |10⟩ states,
/// and e^(i(2φ+θ)) on the |11⟩ state.
///
/// The matrix form is:
/// PSCP = [ [ 1,      0,      0,           0 ],
///          [ 0, e^(iφ),      0,           0 ],
///          [ 0,      0, e^(iφ),           0 ],
///          [ 0,      0,      0, e^(i(2φ+θ)) ] ]
///
/// # Parameters
/// - `control`: The control qubit index
/// - `target`: The target qubit index
/// - `theta`: The additional phase parameter θ
/// - `phi`: The base phase shift parameter φ
pub struct PhaseShiftedControlledPhase {
    control: Qubit,
    target: Qubit,
    theta: f64,
    phi: f64,
}

impl PhaseShiftedControlledPhase {
    pub fn new(control: Qubit, target: Qubit, theta: f64, phi: f64) -> Self {
        Self {
            control,
            target,
            theta,
            phi,
        }
    }
}

impl QuantumGate for PhaseShiftedControlledPhase {
    /// construct the 4x4 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let phi: f64 = self.phi;
        let theta: f64 = self.theta;
        let c: f64 = phi.cos();
        let s: f64 = phi.sin();
        let cos: f64 = (2.0 * phi + theta).cos();
        let sin: f64 = (2.0 * phi + theta).sin();
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(c, s),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(c, s),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(cos, sin),
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "PSCP(control={}, target={}, theta={}, phi={})",
            self.control, self.target, self.theta, self.phi
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
/// Represents a controlled rotation about the X-axis gate.
///
/// This gate performs a rotation about the X-axis on the target qubit when
/// the control qubit is in the |1⟩ state. The rotation angle is θ.
///
/// The matrix form is:
/// CRX = [ [ 1, 0,           0,           0 ],
///         [ 0, 1,           0,           0 ],
///         [ 0, 0,    cos(θ/2), -i·sin(θ/2) ],
///         [ 0, 0, -i·sin(θ/2),    cos(θ/2) ] ]
///
/// # Parameters
/// - `control`: The control qubit index
/// - `target`: The target qubit index
/// - `theta`: The rotation angle θ in radians
pub struct ControlledRotateX {
    control: Qubit,
    target: Qubit,
    theta: f64,
}

impl ControlledRotateX {
    pub fn new(control: Qubit, target: Qubit, theta: f64) -> Self {
        Self {
            control,
            target,
            theta,
        }
    }
}

impl QuantumGate for ControlledRotateX {
    /// construct the 4x4 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(c, 0.0),
                Complex::new(0.0, -1.0 * s),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -1.0 * s),
                Complex::new(c, 0.0),
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "CRX(control={}, target={}, theta={})",
            self.control, self.target, self.theta
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
/// Represents a controlled rotation in the XY-plane gate.
///
/// This gate performs a rotation in the XY-plane on the target qubit when
/// the control qubit is in the |1⟩ state. The rotation is parameterized by
/// angle θ and axis direction φ in the XY-plane.
///
/// The matrix form involves rotations around the axis (cos(φ), sin(φ), 0):
/// CRXY = [ [ 1, 0,               0,                 0 ],
///          [ 0, 1,               0,                 0 ],
///          [ 0, 0,        cos(θ/2), -sin(θ/2)·e^(-iφ) ],
///          [ 0, 0, sin(θ/2)·e^(iφ),          cos(θ/2) ] ]
///
/// # Parameters
/// - `control`: The control qubit index
/// - `target`: The target qubit index
/// - `theta`: The rotation angle θ in radians
/// - `phi`: The axis angle φ in the XY-plane
pub struct ControlledRotateXY {
    control: Qubit,
    target: Qubit,
    theta: f64,
    phi: f64,
}

impl ControlledRotateXY {
    pub fn new(control: Qubit, target: Qubit, theta: f64, phi: f64) -> Self {
        Self {
            control,
            target,
            theta,
            phi,
        }
    }
}

impl QuantumGate for ControlledRotateXY {
    /// construct the 4x4 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        let vx: f64 = (self.phi).cos();
        let vy: f64 = (self.phi).sin();
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(c, 0.0),
                Complex::new(-1.0 * s * vy, -1.0 * s * vx),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(s * vy, -1.0 * s * vx),
                Complex::new(c, 0.0),
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "CRXY(control={}, target={}, theta={}, phi={})",
            self.control, self.target, self.theta, self.phi
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
/// Represents an XY gate (iSWAP-like interaction).
///
/// This gate implements an XY interaction between two qubits, creating
/// entanglement by mixing the |01⟩ and |10⟩ states with a rotation angle θ.
/// It's related to the iSWAP gate family.
///
/// The matrix form is:
/// XY = [ [ 1,          0,          0, 0 ],
///        [ 0,   cos(θ/2), i·sin(θ/2), 0 ],
///        [ 0, i·sin(θ/2),   cos(θ/2), 0 ],
///        [ 0,          0,          0, 1 ] ]
///
/// # Parameters
/// - `control`: The first qubit index
/// - `target`: The second qubit index
/// - `theta`: The interaction strength θ in radians
pub struct XY {
    control: Qubit,
    target: Qubit,
    theta: f64,
}

impl XY {
    pub fn new(control: Qubit, target: Qubit, theta: f64) -> Self {
        Self {
            control,
            target,
            theta,
        }
    }
}

impl QuantumGate for XY {
    /// construct the 4x4 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(c, 0.0),
                Complex::new(0.0, s),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, s),
                Complex::new(c, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "XY(control={}, target={}, theta={})",
            self.control, self.target, self.theta
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
/// Represents a controlled phase shift gate.
///
/// This gate applies a phase shift e^(iθ) to the |11⟩ state when both
/// control and target qubits are in the |1⟩ state, leaving other states unchanged.
///
/// The matrix form is:
/// CPS = [ [ 1, 0, 0,      0 ],
///         [ 0, 1, 0,      0 ],
///         [ 0, 0, 1,      0 ],
///         [ 0, 0, 0, e^(iθ) ] ]
///
/// # Parameters
/// - `control`: The control qubit index
/// - `target`: The target qubit index
/// - `theta`: The phase shift angle θ in radians
pub struct ControlledPhaseShift {
    control: Qubit,
    target: Qubit,
    theta: f64,
}

impl ControlledPhaseShift {
    pub fn new(control: Qubit, target: Qubit, theta: f64) -> Self {
        Self {
            control,
            target,
            theta,
        }
    }
}

impl QuantumGate for ControlledPhaseShift {
    /// construct the 4x4 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta).cos();
        let s: f64 = (self.theta).sin();
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(c, s),
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "CPS(control={}, target={}, theta={})",
            self.control, self.target, self.theta
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
/// Represents a variable Mølmer-Sørensen XX gate.
///
/// This gate implements a parametric version of the Mølmer-Sørensen gate,
/// which creates entanglement through XX interactions. It's commonly used
/// in trapped ion quantum computing.
///
/// The matrix form is:
/// VMSXX = [ [    cos(θ/2),           0,           0, -i·sin(θ/2) ],
///           [           0,    cos(θ/2), -i·sin(θ/2),           0 ],
///           [           0, -i·sin(θ/2),    cos(θ/2),           0 ],
///           [ -i·sin(θ/2),           0,           0,    cos(θ/2) ] ]
///
/// # Parameters
/// - `control`: The first qubit index
/// - `target`: The second qubit index
/// - `theta`: The interaction strength θ in radians
pub struct VariableMSXX {
    control: Qubit,
    target: Qubit,
    theta: f64,
}

impl VariableMSXX {
    pub fn new(control: Qubit, target: Qubit, theta: f64) -> Self {
        Self {
            control,
            target,
            theta,
        }
    }
}

impl QuantumGate for VariableMSXX {
    /// construct the 4x4 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.theta / 2.0).cos();
        let s: f64 = (self.theta / 2.0).sin();
        array![
            [
                Complex::new(c, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -1.0 * s),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(c, 0.0),
                Complex::new(0.0, -1.0 * s),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -1.0 * s),
                Complex::new(c, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, -1.0 * s),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(c, 0.0),
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "VMSXX(control={}, target={}, theta={})",
            self.control, self.target, self.theta
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
/// Represents a Givens rotation gate.
///
/// This gate implements a Givens rotation in the computational basis,
/// which rotates the |01⟩ and |10⟩ subspace by angle θ with an additional
/// phase parameter φ. It's useful for quantum chemistry and optimization algorithms.
///
/// The matrix form is:
/// GR = [ [ 1,              0,      0,      0 ],
///        [ 0,  cos(θ)·e^(iφ), sin(θ),      0 ],
///        [ 0, -sin(θ)·e^(iφ), cos(θ),      0 ],
///        [ 0,              0,      0, e^(iφ) ] ]
///
/// # Parameters
/// - `control`: The first qubit index
/// - `target`: The second qubit index
/// - `theta`: The rotation angle θ in radians
pub struct GivensRotation {
    control: Qubit,
    target: Qubit,
    theta: f64,
    phi: f64,
}

impl GivensRotation {
    pub fn new(control: Qubit, target: Qubit, theta: f64, phi: f64) -> Self {
        Self {
            control,
            target,
            theta,
            phi,
        }
    }
}

impl QuantumGate for GivensRotation {
    /// construct the 4x4 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let ct: f64 = (self.theta).cos();
        let st: f64 = (self.theta).sin();
        let cp: f64 = (self.phi).cos();
        let sp: f64 = (self.phi).sin();
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(ct * cp, ct * sp),
                Complex::new(st, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(-1.0 * st * cp, -1.0 * st * sp),
                Complex::new(ct, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(cp, sp),
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "GR(control={}, target={}, theta={}, phi={})",
            self.control, self.target, self.theta, self.phi
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
/// Represents a Givens rotation gate with little-endian qubit ordering.
///
/// This is similar to the standard Givens rotation but uses little-endian
/// qubit ordering conventions. The rotation acts on the |01⟩ and |10⟩ subspace
/// with different matrix structure to account for the bit ordering.
///
/// The matrix form is:
/// GRLE = [ [ 1,              0,             0,      0 ],
///          [ 0,         cos(θ),        sin(θ),      0 ],
///          [ 0, -sin(θ)·e^(iφ), cos(θ)·e^(iφ),      0 ],
///          [ 0,              0,             0, e^(iφ) ] ]
///
/// # Parameters
/// - `control`: The first qubit index
/// - `target`: The second qubit index
/// - `theta`: The rotation angle θ in radians
/// - `phi`: The phase parameter φ in radians
pub struct GivensRotationLittleEndian {
    control: Qubit,
    target: Qubit,
    theta: f64,
    phi: f64,
}

impl GivensRotationLittleEndian {
    pub fn new(control: Qubit, target: Qubit, theta: f64, phi: f64) -> Self {
        Self {
            control,
            target,
            theta,
            phi,
        }
    }
}

impl QuantumGate for GivensRotationLittleEndian {
    /// construct the 4x4 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let ct: f64 = (self.theta).cos();
        let st: f64 = (self.theta).sin();
        let cp: f64 = (self.phi).cos();
        let sp: f64 = (self.phi).sin();
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(ct, 0.0),
                Complex::new(st, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(-1.0 * st * cp, -1.0 * st * sp),
                Complex::new(ct * cp, ct * sp),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(cp, sp),
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "GRLE(control={}, target={}, theta={}, phi={})",
            self.control, self.target, self.theta, self.phi
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
/// Represents a Bogoliubov transformation gate.
///
/// This gate implements a Bogoliubov transformation used in quantum many-body
/// physics and quantum field theory. It's particularly relevant for fermionic
/// systems and superconductivity simulations.
///
/// The transformation is parameterized by a complex number Δ = δᵣₑ + i·δᵢₘ,
/// where |Δ| determines the mixing strength and arg(Δ) determines the phase.
///
/// The matrix form involves rotations by |Δ| and phases by arg(Δ):
/// B = [ [               cos(|Δ|), 0, 0, -sin(|Δ|)·e^(i·arg(Δ)) ],
///       [                      0, 1, 0,                      0 ],
///       [                      0, 0, 1,                      0 ],
///       [ sin(|Δ|)·e^(-i·arg(Δ)), 0, 0,               cos(|Δ|) ] ]
///
/// # Parameters
/// - `control`: The first qubit index
/// - `target`: The second qubit index
/// - `delta_re`: Real part of the complex parameter Δ
/// - `delta_im`: Imaginary part of the complex parameter Δ
pub struct Bogoliubov {
    control: Qubit,
    target: Qubit,
    delta_re: f64,
    delta_im: f64,
}

impl Bogoliubov {
    pub fn new(control: Qubit, target: Qubit, delta_re: f64, delta_im: f64) -> Self {
        Self {
            control,
            target,
            delta_re,
            delta_im,
        }
    }
}

impl QuantumGate for Bogoliubov {
    /// construct the 4x4 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let delta: Complex = Complex::new(self.delta_re, self.delta_im);
        let dn: f64 = delta.norm();
        let da: f64 = delta.arg();
        array![
            [
                Complex::new(dn.cos(), 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(-1.0 * dn.sin() * da.sin(), dn.sin() * da.cos()),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0)
            ],
            [
                Complex::new(dn.sin() * da.sin(), dn.sin() * da.cos()),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(dn.cos(), 0.0),
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "B(control={}, target={}, delta(Re)={}, delta(Im)={})",
            self.control, self.target, self.delta_re, self.delta_im
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
/// Represents a parametric magnetic interaction gate between two qubits.
///
/// This gate implements a controlled rotation that depends on the state of the control qubit.
/// When the control qubit is in state |1⟩, it applies a rotation in the X-Y plane to the target qubit.
/// The strength parameter determines the angle of rotation.
///
/// The gate acts as:
/// - |00⟩ → |00⟩ (no change)
/// - |01⟩ → cos(θ)|01⟩ - i·sin(θ)|10⟩ (rotation when control is |0⟩, target is |1⟩)
/// - |10⟩ → -i·sin(θ)|01⟩ + cos(θ)|10⟩ (rotation when control is |1⟩, target is |0⟩)
/// - |11⟩ → |11⟩ (no change)
///
/// The 4x4 unitary matrix is:
///
/// PMI = [ [ 1,    0,    0, 0 ],
///         [ 0,    c, -i·s, 0 ],
///         [ 0, -i·s,    c, 0 ],
///         [ 0,    0,    0, 1 ] ]
///
/// # Parameters
/// * `control` - The control qubit
/// * `target` - The target qubit that receives the conditional rotation
/// * `strength` - The rotation angle parameter (in radians)
pub struct PMInteraction {
    control: Qubit,
    target: Qubit,
    strength: f64,
}

impl PMInteraction {
    pub fn new(control: Qubit, target: Qubit, strength: f64) -> Self {
        Self {
            control,
            target,
            strength,
        }
    }
}

impl QuantumGate for PMInteraction {
    /// construct the 4x4 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let c: f64 = (self.strength).cos();
        let s: f64 = (self.strength).sin();
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(c, 0.0),
                Complex::new(0.0, -1.0 * s),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -1.0 * s),
                Complex::new(c, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "PMI(control={}, target={}, strength={})",
            self.control, self.target, self.strength
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
/// Represents a complex parametric magnetic interaction gate between two qubits.
///
/// This is an extension of PMInteraction that allows for complex-valued strength parameters,
/// enabling more general rotations in the complex plane. The gate applies conditional rotations
/// based on both the magnitude and phase of the complex strength parameter.
///
/// The complex strength parameter is decomposed into:
/// - Magnitude (sn): determines the rotation angle
/// - Argument/Phase (sa): determines the phase shift
///
/// The gate transformation involves:
/// - |00⟩ and |11⟩ states remain unchanged
/// - |01⟩ and |10⟩ states undergo complex rotation mixing
///
/// The 4x4 unitary matrix is:
///
/// CPMI = [ [ 1,                     0,                       0, 0 ],
///          [ 0,              cos(|z|), -sin(|z|)·e^(-i·arg(z)), 0 ],
///          [ 0, sin(|z|)·e^(i·arg(z)),                cos(|z|), 0 ],
///          [ 0,                     0,                       0, 1 ] ]
///
/// # Parameters
/// * `control` - The control qubit
/// * `target` - The target qubit that receives the conditional rotation
/// * `strength_re` - Real part of the complex strength parameter
/// * `strength_im` - Imaginary part of the complex strength parameter
pub struct ComplexPMInteraction {
    control: Qubit,
    target: Qubit,
    strength_re: f64,
    strength_im: f64,
}

impl ComplexPMInteraction {
    pub fn new(control: Qubit, target: Qubit, strength_re: f64, strength_im: f64) -> Self {
        Self {
            control,
            target,
            strength_re,
            strength_im,
        }
    }
}

impl QuantumGate for ComplexPMInteraction {
    /// construct the 4x4 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let strength: Complex = Complex::new(self.strength_re, self.strength_im);
        let sn: f64 = strength.norm();
        let sa: f64 = strength.arg();
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(sn.cos(), 0.0),
                Complex::new(-1.0 * sn.sin(), -1.0 * sn.sin() * sa.cos()),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(sn.sin() * sa.sin(), -1.0 * sn.sin() * sa.cos()),
                Complex::new(sn.cos(), 0.0),
                Complex::new(0.0, 0.0)
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ]
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "CPMI(control={}, target={}, strength(Re)={}, strength(Im)={})",
            self.control, self.target, self.strength_re, self.strength_im
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
/// Represents a general spin-spin interaction gate between two qubits.
///
/// This gate models the interaction between two quantum spins with coupling
/// strengths in all three spatial directions (x, y, z). It's commonly used
/// to simulate Heisenberg-type interactions in quantum many-body systems.
///
/// The 4x4 unitary matrix has a complex structure involving trigonometric
/// functions of the coupling parameters:
///
/// SI = [ [    cos(x-y)·e^(-iz),                 0,                 0, -i·sin(x-y)·e^(-iz) ],
///        [                   0,   cos(x+y)·e^(iz), i·sin(x+y)·e^(iz),                   0 ],
///        [                   0, i·sin(x+y)·e^(iz),   cos(x+y)·e^(iz),                   0 ],
///        [ -i·sin(x-y)·e^(-iz),                 0,                 0,    cos(x-y)·e^(-iz) ] ]
///
/// This gate is particularly useful for simulating quantum magnetism and
/// spin chain dynamics.
///
/// # Arguments
/// * `control` - The first qubit in the spin interaction
/// * `target` - The second qubit in the spin interaction
/// * `x` - Coupling strength in x-direction (σₓ ⊗ σₓ term)
/// * `y` - Coupling strength in y-direction (σᵧ ⊗ σᵧ term)
/// * `z` - Coupling strength in z-direction (σᵤ ⊗ σᵤ term)
///
/// The transformation is based on the Hamiltonian:
/// H = x·(σₓ ⊗ σₓ) + y·(σᵧ ⊗ σᵧ) + z·(σᵤ ⊗ σᵤ)
pub struct SpinInteraction {
    control: Qubit,
    target: Qubit,
    x: f64,
    y: f64,
    z: f64,
}

impl SpinInteraction {
    pub fn new(control: Qubit, target: Qubit, x: f64, y: f64, z: f64) -> Self {
        Self {
            control,
            target,
            x,
            y,
            z,
        }
    }
}

impl QuantumGate for SpinInteraction {
    /// construct the 4x4 matrix representing the gate.
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let cm: f64 = (self.x - self.y).cos();
        let cp: f64 = (self.x + self.y).cos();
        let sm: f64 = (self.x - self.y).sin();
        let sp: f64 = (self.x + self.y).sin();
        let cz: f64 = self.z.cos();
        let sz: f64 = self.z.sin();
        array![
            [
                Complex::new(cm * cz, (-1.0) * cm * sz),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new((-1.0) * sm * sz, (-1.0) * sm * cz)
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(cp * cz, cp * sz),
                Complex::new(sp * sz, (-1.0) * sp * cz),
                Complex::new(0.0, 0.0)
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(sp * sz, (-1.0) * sp * cz),
                Complex::new(cp * cz, cp * sz),
                Complex::new(0.0, 0.0)
            ],
            [
                Complex::new((-1.0) * sm * sz, (-1.0) * sm * cz),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(cm * cz, (-1.0) * cm * sz)
            ],
        ]
    }

    /// returns the alias name representing the gate.
    fn name(&self) -> String {
        format!(
            "SI(control={}, target={}, x={}, y={}, z={})",
            self.control, self.target, self.x, self.y, self.z
        )
    }

    /// construct targets for quantum state.
    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}
