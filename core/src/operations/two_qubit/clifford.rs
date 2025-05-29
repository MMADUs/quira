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
    operations::QuantumGate,
    types::{Complex, Matrix, Qubit},
};

use ndarray::array;

#[derive(Debug, Clone)]
pub struct ControlledNot {
    control: Qubit,
    target: Qubit,
}

impl ControlledNot {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for ControlledNot {
    fn unitary_matrix(&self) -> Matrix<Complex> {
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
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
            ]
        ]
    }

    fn name(&self) -> String {
        format!("CNOT(control={}, target={})", self.control, self.target)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
pub struct SWAP {
    control: Qubit,
    target: Qubit,
}

impl SWAP {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for SWAP {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
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
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
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
        format!("SWAP(control={}, target={})", self.control, self.target)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
pub struct ISWAP {
    control: Qubit,
    target: Qubit,
}

impl ISWAP {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for ISWAP {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 1.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 1.0),
                Complex::new(0.0, 0.0),
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
        format!("ISWAP(control={}, target={})", self.control, self.target)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
pub struct FSWAP {
    control: Qubit,
    target: Qubit,
}

impl FSWAP {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for FSWAP {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        array![
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
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
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(-1.0, 0.0)
            ]
        ]
    }

    fn name(&self) -> String {
        format!("FSWAP(control={}, target={})", self.control, self.target)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
pub struct ControlledPauliY {
    control: Qubit,
    target: Qubit,
}

impl ControlledPauliY {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for ControlledPauliY {
    fn unitary_matrix(&self) -> Matrix<Complex> {
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
                Complex::new(0.0, -1.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -1.0),
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 1.0),
                Complex::new(0.0, 0.0),
            ]
        ]
    }

    fn name(&self) -> String {
        format!("CY(control={}, target={})", self.control, self.target)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
pub struct ControlledPauliZ {
    control: Qubit,
    target: Qubit,
}

impl ControlledPauliZ {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for ControlledPauliZ {
    fn unitary_matrix(&self) -> Matrix<Complex> {
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
                Complex::new(-1.0, 0.0)
            ]
        ]
    }

    fn name(&self) -> String {
        format!("CZ(control={}, target={})", self.control, self.target)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}

#[derive(Debug, Clone)]
pub struct EchoCrossResonance {
    control: Qubit,
    target: Qubit,
}

impl EchoCrossResonance {
    pub fn new(control: Qubit, target: Qubit) -> Self {
        Self { control, target }
    }
}

impl QuantumGate for EchoCrossResonance {
    fn unitary_matrix(&self) -> Matrix<Complex> {
        let matrix: Matrix<Complex> = array![
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 1.0)
            ],
            [
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 1.0),
                Complex::new(1.0, 0.0)
            ],
            [
                Complex::new(1.0, 0.0),
                Complex::new(0.0, -1.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0)
            ],
            [
                Complex::new(0.0, -1.0),
                Complex::new(1.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0)
            ]
        ];
        // divide
        matrix / 2.0_f64.sqrt()
    }

    fn name(&self) -> String {
        format!("ECR(control={}, target={})", self.control, self.target)
    }

    fn construct_targets(&self) -> Vec<Qubit> {
        vec![self.target, self.control]
    }
}
