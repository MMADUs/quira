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
    target: Qubit,
    alpha_re: f64,
    alpha_im: f64,
    beta_re: f64,
    beta_im: f64,
    global_phase: f64,
}

impl Arbitrary {
    pub fn new(
        target: Qubit,
        alpha_re: f64,
        alpha_im: f64,
        beta_re: f64,
        beta_im: f64,
        global_phase: f64,
    ) -> Self {
        if alpha_re == 0.0 && alpha_im == 0.0 && beta_re == 0.0 && beta_im == 0.0 {
            panic!("Invalid gate parameters: all alpha and beta components are zero");
        }

        let norm: f64 =
            alpha_re.powf(2.0) + alpha_im.powf(2.0) + beta_re.powf(2.0) + beta_im.powf(2.0);

        if (norm - 1.0).abs() > 1e-6 {
            panic!(
                "Unitarity violation: |α|² + |β|² must equal 1.0, got {}",
                norm
            );
        }

        Self {
            target,
            alpha_re,
            alpha_im,
            beta_re,
            beta_im,
            global_phase,
        }
    }
}

impl QuantumGate for Arbitrary {
    fn unitary_matrix(&self) -> Matrix<Complex> {
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

    fn name(&self) -> String {
        let alpha = Complex::new(self.alpha_re, self.alpha_im);
        let beta = Complex::new(self.beta_re, self.beta_im);
        format!(
            "AR(target={}, alpha={}, beta={}, global_phase={})",
            self.target, alpha, beta, self.global_phase
        )
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target]
    }

    fn enumerated(&self) -> crate::GateType {
        GateType::SingleQubit(SingleQubitType::Arbitrary(Self::new(
            self.target,
            self.alpha_re,
            self.alpha_im,
            self.beta_re,
            self.beta_im,
            self.global_phase,
        )))
    }
}

impl SingleQubitGate for Arbitrary {
    fn alpha_re(&self) -> f64 {
        self.alpha_re
    }

    fn alpha_im(&self) -> f64 {
        self.alpha_im
    }

    fn beta_re(&self) -> f64 {
        self.beta_re
    }

    fn beta_im(&self) -> f64 {
        self.beta_im
    }

    fn global_phase(&self) -> f64 {
        self.global_phase
    }
}
