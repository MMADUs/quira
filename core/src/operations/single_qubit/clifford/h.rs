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
use crate::{Complex, GateType, Matrix, Qubit, constant::PI};

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

    pub fn power(&self, exponent: f64) -> HPowGate {
        HPowGate::new(self.target, exponent)
    }
}

impl QuantumGate for Hadamard {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let f = 1.0 / ((2.0_f64).sqrt());
        array![
            [Complex::new(f, 0.0), Complex::new(f, 0.0)],
            [Complex::new(f, 0.0), Complex::new(-1.0 * f, 0.0)]
        ]
    }

    fn name(&self) -> String {
        format!("H(target={})", self.target)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::Hadamard(Self::new(self.target)))
    }
}

impl SingleQubitGate for Hadamard {
    fn alpha_re(&self) -> f64 {
        0.0
    }

    fn alpha_im(&self) -> f64 {
        -1.0 / ((2.0_f64).sqrt())
    }

    fn beta_re(&self) -> f64 {
        0.0
    }

    fn beta_im(&self) -> f64 {
        -1.0 / ((2.0_f64).sqrt())
    }

    fn global_phase(&self) -> f64 {
        PI / 2.0
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

    fn enumerated(&self) -> GateType {
        GateType::SingleQubit(SingleQubitType::HPowGate(Self::new(
            self.target,
            self.exponent,
        )))
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
