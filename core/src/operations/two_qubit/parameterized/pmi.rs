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

    fn name(&self) -> String {
        format!(
            "PMI(control={}, target={}, strength={})",
            self.control, self.target, self.strength
        )
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::PMInteraction(Self::new(
            self.control,
            self.target,
            self.strength,
        )))
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

    fn name(&self) -> String {
        format!(
            "CPMI(control={}, target={}, strength(Re)={}, strength(Im)={})",
            self.control, self.target, self.strength_re, self.strength_im
        )
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }

    fn enumerated(&self) -> GateType {
        GateType::TwoQubit(TwoQubitType::ComplexPMInteraction(Self::new(
            self.control,
            self.target,
            self.strength_re,
            self.strength_im,
        )))
    }
}
