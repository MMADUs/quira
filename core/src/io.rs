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

use crate::bit::QuantumBit;
use crate::operations::GateType;
use crate::types::{Complex, Matrix, Vector};

/// every kernel (quantum state) implements backend operation
/// for executing circuit in the backend
pub trait BackendOperation {
    /// reset the kernel state
    fn reset_state(&mut self);

    /// expand the kernel state
    fn expand_state(&mut self, state: Vector<Complex>);

    /// measure qubit on the kernel
    fn measure_qubit(&mut self, qubit: usize) -> bool;

    /// measure all qubit on the kernel
    fn measure_all(&mut self) -> Vec<bool>;
}

/// every kernel (quantum state) implements quantum state
/// for applying quantum gate to state
pub trait QuantumState: BackendOperation + Send + Sync {
    /// return number of qubit in state
    fn num_qubits(&self) -> usize;

    /// apply any N quantum gate to kernel
    fn apply(&mut self, u: Matrix<Complex>, enumerated: GateType);

    /// dynamic clone
    fn box_clone(&self) -> Box<dyn QuantumState>;
}

/// implement clone for Box dynamic quantum state
impl Clone for Box<dyn QuantumState> {
    fn clone(&self) -> Self {
        self.box_clone()
    }
}

/// Generic quantum debugger trait
pub trait QuantumDebugger {
    /// Type of internal full state (Vec<Complex> for statevec, Matrix for density)
    type StateType;

    /// Type used to represent a "bit state" — either amplitude or probability
    type BitStateType;

    /// Type used to represent a "qubit state" — either amplitude pair or probability pair
    type QubitStateType;

    /// Returns the entire state representation
    fn entire_state(&self, filter_zero: bool) -> Self::StateType;

    /// Returns info about a specific bitstring
    fn bit_state(&self, bits: &[bool]) -> Self::BitStateType;

    /// Returns info about a specific qubit (e.g., marginal probability)
    fn qubit_state(&self, qubit: &QuantumBit) -> Self::QubitStateType;
}
