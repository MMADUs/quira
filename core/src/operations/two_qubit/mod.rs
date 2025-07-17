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

pub mod parameterized;
pub use parameterized::*;

pub mod extra;
pub use extra::*;

pub enum TwoQubitType {
    // clifford
    ControlledNot(ControlledNot),
    SWAP(SWAP),
    ISWAP(ISWAP),
    FSWAP(FSWAP),
    ControlledPauliY(ControlledPauliY),
    ControlledPauliZ(ControlledPauliZ),
    EchoCrossResonance(EchoCrossResonance),
    // todo!
    DoubleControlledNot,
    ControlledSqrtX,
    ControlledHadamard,

    // non-clifford
    SqrtISWAP(SqrtISWAP),
    InvSqrtISWAP(InvSqrtISWAP),

    // parameterized
    PhaseShiftedControlledZ(PhaseShiftedControlledZ),
    PhaseShiftedControlledPhase(PhaseShiftedControlledPhase),
    ControlledRotateX(ControlledRotateX),
    ControlledRotateXY(ControlledRotateXY),
    XY(XY),
    ControlledPhaseShift(ControlledPhaseShift),
    VariableMSXX(VariableMSXX),
    GivensRotation(GivensRotation),
    GivensRotationLittleEndian(GivensRotationLittleEndian),
    Bogoliubov(Bogoliubov),
    PMInteraction(PMInteraction),
    ComplexPMInteraction(ComplexPMInteraction),
    SpinInteraction(SpinInteraction),
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
    MolmerSorensenXX(MolmerSorensenXX),
    Qsim(Qsim),
    Fsim(Fsim),
}
