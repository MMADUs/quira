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

#[cfg(test)]
mod tests;

pub mod clifford;
pub use clifford::*;

pub mod non_clifford;
pub use non_clifford::*;

pub mod parameterized;
pub use parameterized::*;

pub mod special;
pub use special::*;

pub mod trivial;
pub use trivial::*;

pub mod universal;
pub use universal::*;

use super::QuantumGate;

/// the composition of kak decomposition
pub struct KakDecomposition {
    pub global_phase: f64,
    pub k_vector: [f64; 3],
    pub u_before: Vec<Box<dyn QuantumGate>>,
    pub u_after: Vec<Box<dyn QuantumGate>>,
}

pub enum TwoQubitOperation {
    // clifford
    ControlledNot,
    SWAP,
    ISWAP,
    FSWAP,
    ControlledPauliY,
    ControlledPauliZ,
    EchoCrossResonance,
    // todo!
    DoubleControlledNot,
    ControlledSqrtX,
    ControlledHadamard,

    // non-clifford
    SqrtISWAP,
    InvSqrtISWAP,

    // parameterized
    PhaseShiftedControlledZ,
    PhaseShiftedControlledPhase,
    ControlledRotateX,
    ControlledRotateXY,
    XY,
    ControlledPhaseShift,
    VariableMSXX,
    GivensRotation,
    GivensRotationLittleEndian,
    Bogoliubov,
    PMInteraction,
    ComplexPMInteraction,
    SpinInteraction,
    // todo!
    RotationXX,
    RotationYY,
    RotationZZ,
    RotationZX,
    SwapPowGate, // check
    XXPower,
    YYPower,
    ZZPower,
    ControlledZPower,
    ControlledNotPower,
    ISwapPowGate, // check

    // universal? todo!
    CU3Gate,
    CU1Gate, // ControlledPhaseShift?
    CUGate,

    // special
    MolmerSorensenXX,
    Qsim,
    Fsim,
}
