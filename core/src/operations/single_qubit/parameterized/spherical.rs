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

    fn name(&self) -> String {
        format!(
            "RAS(target={}, theta={:.4}, sphere_theta={:.4}, sphere_phi{:.4})",
            self.target, self.theta, self.spherical_theta, self.spherical_phi
        )
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::RotateAroundSphericalAxis(Self::new(
            self.target,
            self.theta,
            self.spherical_theta,
            self.spherical_phi,
        )))
    }
}

impl SingleQubitGate for RotateAroundSphericalAxis {
    fn alpha_re(&self) -> f64 {
        (self.theta / 2.0).cos()
    }

    fn alpha_im(&self) -> f64 {
        let s = (self.theta / 2.0).sin();
        let vz = (self.spherical_theta).cos();
        (-1.0) * s * vz
    }

    fn beta_re(&self) -> f64 {
        let s = (self.theta / 2.0).sin();
        let vy = (self.spherical_phi).sin();
        let st = (self.spherical_theta).sin();
        s * vy * st
    }

    fn beta_im(&self) -> f64 {
        let s = (self.theta / 2.0).sin();
        let vx = (self.spherical_phi).cos();
        let st = (self.spherical_theta).sin();
        (-1.0) * s * vx * st
    }

    fn global_phase(&self) -> f64 {
        0.0
    }
}
