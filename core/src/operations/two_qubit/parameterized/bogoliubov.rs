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
use crate::operations::two_qubit::TwoQubitType;
use crate::{Complex, GateType, Matrix, Qubit};

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

    fn name(&self) -> String {
        format!(
            "B(control={}, target={}, delta(Re)={}, delta(Im)={})",
            self.control, self.target, self.delta_re, self.delta_im
        )
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::Bogoliubov(Self::new(
            self.control,
            self.target,
            self.delta_re,
            self.delta_im,
        )))
    }
}
