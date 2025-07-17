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

#[cfg(test)]
mod tests;

pub mod clifford;
pub mod parameterized;
pub mod extra;

pub use clifford::*;
pub use parameterized::*;
pub use extra::*;

/// enum definitions
pub enum SingleQubitType {
    // parameterized
    RotateX(RotateX),
    RotateY(RotateY),
    RotateZ(RotateZ),
    RotateXY(RotateXY),
    U0Gate(U0Gate),
    U1Gate(U1Gate),
    U2Gate(U2Gate),
    U3Gate(U3Gate),
    XPowGate(XPowGate),
    YPowGate(YPowGate),
    ZPowGate(ZPowGate),
    PhasedXPowGate(PhasedXPowGate),
    HPowGate(HPowGate),
    RotateAroundSphericalAxis(RotateAroundSphericalAxis),

    // clifford
    PauliX(PauliX),
    PauliY(PauliY),
    PauliZ(PauliZ),
    SqrtPauliX(SqrtPauliX),
    InvSqrtPauliX(InvSqrtPauliX),
    SqrtPauliY(SqrtPauliY),
    InvSqrtPauliY(InvSqrtPauliY),
    Hadamard(Hadamard),
    SGate(SGate),
    InvSGate(InvSGate),

    // non-clifford
    TGate(TGate),
    InvTGate(InvTGate),

    // universal
    UGate(UGate),
    Arbitrary(Arbitrary),

    // special
    GPi(GPi),
    GPi2(GPi2),

    // trivial
    Identity(Identity),
}
