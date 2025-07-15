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

pub mod single_qubit;
pub mod two_qubit;

use single_qubit::SingleQubitType;
use two_qubit::TwoQubitType;

use crate::{
    endian::{QubitIndexing, expand_unitary},
    kernel::QuantumState,
    types::{Complex, Matrix, Qubit},
};

pub enum GateType {
    SingleQubit(SingleQubitType),
    TwoQubit(TwoQubitType),
}

pub trait QuantumGate: Send + Sync {
    /// Get the unitary matrix representation of the gate
    fn unitary_matrix(&self) -> Matrix<Complex>;

    /// Quantum gate name representation
    fn name(&self) -> String;

    /// Get the target qubits this gate operates on
    fn construct_targets(&self) -> Vec<Qubit>;

    /// Categorized the gate into enum types
    fn enumerated(&self) -> GateType;
}

// Helper trait (optional)
pub trait GateApplyExt {
    fn apply<T: QuantumState>(&self, state: &mut T, indexing: &QubitIndexing);
}

// Implement for all references to types that implement QuantumGate
impl<T: QuantumGate + ?Sized> GateApplyExt for T {
    fn apply<U: QuantumState>(&self, state: &mut U, indexing: &QubitIndexing) {
        let n = state.num_qubits();
        let targets = self.construct_targets();
        let unitary = self.unitary_matrix();
        let enumerated = self.enumerated();
        let u = expand_unitary(n, &targets, &unitary, indexing);
        state.apply(u, enumerated);
    }
}
