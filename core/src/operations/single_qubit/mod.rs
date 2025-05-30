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

/// single qubit gate operation.
pub trait SingleQubitGate: Send + Sync {
    /// alpha real Re(α) of the on-diagonal elements of the unitary matrix.
    fn alpha_re(&self) -> f64;

    /// alpha imaginary Im(α) of the on-diagonal elements of the unitary matrix.
    fn alpha_im(&self) -> f64;

    /// beta real Re(β) of the off-diagonal elements of the unitary matrix.
    fn beta_re(&self) -> f64;

    /// beta imaginary Im(β) of the off-diagonal elements of the unitary matrix.
    fn beta_im(&self) -> f64;

    /// the global phase φ of the unitary.
    fn global_phase(&self) -> f64;
}

/// enum definitions
pub enum SingleQubitOperation {
    // parameterized
    RotateX,
    RotateY,
    RotateZ,
    RotateXY,
    PhaseShift0, // rename to U0Gate
    PhaseShift1, // rename to U1Gate
    U2Gate, // todo!
    U3Gate, // todo!
    XPowGate, // todo!
    YPowGate, // todo!
    ZPowGate, // todo!
    PhasedXPowGate, // todo!
    HPowGate, // todo!
    RotateAroundSphericalAxis,

    // clifford
    PauliX,
    PauliY,
    PauliZ,
    SqrtPauliX,
    InvSqrtPauliX,
    SqrtPauliY,
    InvSqrtPauliY,
    Hadamard,
    SGate,
    InvSGate,
    SXGate,
    InvSXGate,

    // non-clifford
    TGate,
    InvTGate,

    // universal
    UGate,
    Arbitrary,

    // special
    GPi,
    GPi2,

    // trivial
    Identity,
}
