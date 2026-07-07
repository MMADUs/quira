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

use rand::Rng;
use std::collections::HashMap;

use crate::bit::QuantumBit;
use crate::io::{BackendOperation, QuantumDebugger, QuantumState};
use crate::operations::GateType;
use crate::operations::single_qubit::SingleQubitType;
use crate::operations::two_qubit::TwoQubitType;
use crate::ops::QuantumGate;
use crate::types::{Complex, Matrix, Vector};

use super::pstring::{PauliOperator, PauliString};

#[derive(Clone)]
/// Standard stabilizer quantum state representation using tableau method
/// Following Qiskit's implementation approach
pub struct StabilizerState {
    /// Number of qubits in the system
    n_qubits: usize,
    /// Stabilizer tableau: 2n x 2n matrix
    /// First n columns are X part, next n columns are Z part
    /// First n rows are destabilizers, next n rows are stabilizers
    tableau: Vec<Vec<bool>>,
    /// Phase vector: false = +1, true = -1
    phases: Vec<bool>,
}

impl StabilizerState {
    /// Initialize empty stabilizer state
    pub fn new() -> Self {
        Self {
            n_qubits: 0,
            tableau: Vec::new(),
            phases: Vec::new(),
        }
    }

    /// Create stabilizer state for n qubits in |0...0⟩ state
    pub fn zeros(n_qubits: usize) -> Self {
        let mut tableau = vec![vec![false; 2 * n_qubits]; 2 * n_qubits];
        let phases = vec![false; 2 * n_qubits];

        // Initialize tableau for |0...0⟩ state
        for i in 0..n_qubits {
            tableau[i][i] = true; // X_i destabilizers
            tableau[n_qubits + i][n_qubits + i] = true; // Z_i stabilizers
        }

        Self {
            n_qubits,
            tableau,
            phases,
        }
    }

    /// Create stabilizer state for n qubits in |1...1⟩ state
    pub fn ones(n_qubits: usize) -> Self {
        let mut state = Self::zeros(n_qubits);
        // Apply X gate to each qubit to flip |0⟩ to |1⟩
        for i in 0..n_qubits {
            state.x(i);
        }
        state
    }

    /// Check if state is empty
    pub fn is_empty(&self) -> bool {
        self.n_qubits == 0
    }

    /// Get number of qubits
    pub fn num_qubits(&self) -> usize {
        self.n_qubits
    }

    /// Check if the state is in a mixed state (not a pure state)
    pub fn is_mixed(&self) -> bool {
        // Check if all stabilizers are independent
        let rank = self.stabilizer_rank();
        rank < self.n_qubits
    }

    /// Get the rank of the stabilizer group
    pub fn stabilizer_rank(&self) -> usize {
        let mut rank = 0;
        for i in self.n_qubits..2 * self.n_qubits {
            let mut has_support = false;
            for j in 0..2 * self.n_qubits {
                if self.tableau[i][j] {
                    has_support = true;
                    break;
                }
            }
            if has_support {
                rank += 1;
            }
        }
        rank
    }

    /// Check if tableau is in canonical form
    pub fn is_valid(&self) -> bool {
        // Check commutation relations
        for i in 0..self.n_qubits {
            for j in 0..self.n_qubits {
                if i != j && !self.commutes(i, j) {
                    return false;
                }
            }
        }
        true
    }

    /// Check if two stabilizers commute
    fn commutes(&self, i: usize, j: usize) -> bool {
        let mut count = 0;
        for k in 0..self.n_qubits {
            if self.tableau[i][k] && self.tableau[j][k + self.n_qubits] {
                count += 1;
            }
            if self.tableau[i][k + self.n_qubits] && self.tableau[j][k] {
                count += 1;
            }
        }
        count % 2 == 0
    }

    // More efficient row multiplication with standard phase calculation
    fn row_mult(&mut self, i: usize, j: usize) {
        let mut phase = 0u8;
        let n = self.n_qubits;

        // Calculate symplectic inner product for phase
        // phase = sum_k (x_i[k] * z_j[k] - z_i[k] * x_j[k]) mod 2
        for k in 0..n {
            let x_i = self.tableau[i][k];
            let z_i = self.tableau[i][k + n];
            let x_j = self.tableau[j][k];
            let z_j = self.tableau[j][k + n];

            // Symplectic inner product contribution
            if x_i && z_j {
                phase ^= 1;
            }
            if z_i && x_j {
                phase ^= 1;
            }
        }

        // Update tableau row i = row i XOR row j
        for k in 0..2 * n {
            self.tableau[i][k] ^= self.tableau[j][k];
        }

        // Update phase: phase_i = phase_i XOR phase_j XOR calculated_phase
        self.phases[i] ^= self.phases[j] ^ (phase == 1);
    }

    /// Apply Pauli-X gate to qubit
    pub fn x(&mut self, qubit: usize) {
        if qubit >= self.n_qubits {
            return;
        }

        // X gate: add phase to rows that have Z on this qubit
        for i in 0..2 * self.n_qubits {
            if self.tableau[i][qubit + self.n_qubits] {
                self.phases[i] = !self.phases[i];
            }
        }
    }

    /// Apply Pauli-Y gate to qubit
    pub fn y(&mut self, qubit: usize) {
        if qubit >= self.n_qubits {
            return;
        }

        // Y gate: add phase to rows that have X or Z on this qubit
        for i in 0..2 * self.n_qubits {
            if self.tableau[i][qubit] || self.tableau[i][qubit + self.n_qubits] {
                self.phases[i] = !self.phases[i];
            }
        }
    }

    /// Apply Pauli-Z gate to qubit
    pub fn z(&mut self, qubit: usize) {
        if qubit >= self.n_qubits {
            return;
        }

        // Z gate: add phase to rows that have X on this qubit
        for i in 0..2 * self.n_qubits {
            if self.tableau[i][qubit] {
                self.phases[i] = !self.phases[i];
            }
        }
    }

    /// Apply Hadamard gate to qubit
    pub fn h(&mut self, qubit: usize) {
        if qubit >= self.n_qubits {
            return;
        }

        // H gate: swap X and Z columns, add phase if both were true
        for i in 0..2 * self.n_qubits {
            let x_val = self.tableau[i][qubit];
            let z_val = self.tableau[i][qubit + self.n_qubits];

            // Swap X and Z
            self.tableau[i][qubit] = z_val;
            self.tableau[i][qubit + self.n_qubits] = x_val;

            // Add phase if both were true (Y -> -Y under H)
            if x_val && z_val {
                self.phases[i] = !self.phases[i];
            }
        }
    }

    /// Apply S gate (phase gate) to qubit
    pub fn s(&mut self, qubit: usize) {
        if qubit >= self.n_qubits {
            return;
        }

        // S gate: X -> Y, Y -> -X, Z -> Z
        for i in 0..2 * self.n_qubits {
            let x_val = self.tableau[i][qubit];
            let z_val = self.tableau[i][qubit + self.n_qubits];

            if x_val && z_val {
                // Y -> -X: remove Z component, add phase
                self.tableau[i][qubit + self.n_qubits] = false;
                self.phases[i] = !self.phases[i];
            } else if x_val && !z_val {
                // X -> Y: add Z component
                self.tableau[i][qubit + self.n_qubits] = true;
            }
        }
    }

    /// Apply CNOT gate (control, target)
    pub fn cnot(&mut self, control: usize, target: usize) {
        if control >= self.n_qubits || target >= self.n_qubits || control == target {
            return; // Add same-qubit check
        }

        // CNOT gate tableau update rules
        for i in 0..2 * self.n_qubits {
            let x_ctrl = self.tableau[i][control];
            let x_tgt = self.tableau[i][target];
            let z_ctrl = self.tableau[i][control + self.n_qubits];
            let z_tgt = self.tableau[i][target + self.n_qubits];

            // Phase update from commutation
            if x_ctrl && z_tgt && !z_ctrl && !x_tgt {
                self.phases[i] = !self.phases[i];
            }

            // Update X components
            self.tableau[i][target] = x_tgt ^ x_ctrl;

            // Update Z components
            self.tableau[i][control + self.n_qubits] = z_ctrl ^ z_tgt;
        }
    }

    /// Apply CZ gate (control, target)
    pub fn cz(&mut self, control: usize, target: usize) {
        // CZ = H_target * CNOT(control, target) * H_target
        self.h(target);
        self.cnot(control, target);
        self.h(target);
    }

    /// Apply a Pauli string to the state
    pub fn apply_pauli_string(&mut self, pauli_string: &PauliString) -> Result<(), String> {
        if pauli_string.len() != self.n_qubits {
            return Err("Pauli string length doesn't match number of qubits".to_string());
        }

        // Apply each Pauli operator
        for (&qubit, &op) in &pauli_string.operators {
            match op {
                PauliOperator::I => {} // Do nothing
                PauliOperator::X => self.x(qubit),
                PauliOperator::Y => self.y(qubit),
                PauliOperator::Z => self.z(qubit),
            }
        }

        // Apply global phase if needed
        if pauli_string.phase() {
            self.apply_global_phase();
        }

        Ok(())
    }

    /// Apply global phase (-1)
    fn apply_global_phase(&mut self) {
        // Global phase affects all stabilizers
        for i in self.n_qubits..2 * self.n_qubits {
            self.phases[i] = !self.phases[i];
        }
    }

    /// Get stabilizer generators as Pauli strings
    pub fn get_stabilizers(&self) -> Vec<PauliString> {
        let mut stabilizers = Vec::new();

        for i in self.n_qubits..2 * self.n_qubits {
            let mut pauli_string = PauliString::new(self.n_qubits);
            pauli_string.set_phase(self.phases[i]);

            for j in 0..self.n_qubits {
                let x_part = self.tableau[i][j];
                let z_part = self.tableau[i][j + self.n_qubits];

                let op = match (x_part, z_part) {
                    (false, false) => PauliOperator::I,
                    (true, false) => PauliOperator::X,
                    (false, true) => PauliOperator::Z,
                    (true, true) => PauliOperator::Y,
                };

                if !op.is_identity() {
                    pauli_string.set(j, op).unwrap();
                }
            }

            stabilizers.push(pauli_string);
        }

        stabilizers
    }

    /// Get destabilizer generators as Pauli strings
    pub fn get_destabilizers(&self) -> Vec<PauliString> {
        let mut destabilizers = Vec::new();

        for i in 0..self.n_qubits {
            let mut pauli_string = PauliString::new(self.n_qubits);
            pauli_string.set_phase(self.phases[i]);

            for j in 0..self.n_qubits {
                let x_part = self.tableau[i][j];
                let z_part = self.tableau[i][j + self.n_qubits];

                let op = match (x_part, z_part) {
                    (false, false) => PauliOperator::I,
                    (true, false) => PauliOperator::X,
                    (false, true) => PauliOperator::Z,
                    (true, true) => PauliOperator::Y,
                };

                if !op.is_identity() {
                    pauli_string.set(j, op).unwrap();
                }
            }

            destabilizers.push(pauli_string);
        }

        destabilizers
    }

    /// Check if a Pauli string is stabilized by this state
    pub fn is_stabilized_by(&self, pauli_string: &PauliString) -> bool {
        // Check if the Pauli string commutes with all stabilizers
        let stabilizers = self.get_stabilizers();
        stabilizers
            .iter()
            .all(|stab| stab.commutes_with(pauli_string))
    }

    /// Add a new stabilizer constraint (for building custom stabilizer states)
    pub fn add_stabilizer(&mut self, pauli_string: &PauliString) -> Result<(), String> {
        if pauli_string.len() != self.n_qubits {
            return Err("Pauli string length doesn't match number of qubits".to_string());
        }

        // Check if it commutes with existing stabilizers
        if !self.is_stabilized_by(pauli_string) {
            return Err("New stabilizer doesn't commute with existing stabilizers".to_string());
        }

        // Find an empty row in the stabilizer part
        for i in self.n_qubits..2 * self.n_qubits {
            let mut is_empty = true;
            for j in 0..2 * self.n_qubits {
                if self.tableau[i][j] {
                    is_empty = false;
                    break;
                }
            }

            if is_empty {
                // Add the new stabilizer
                self.phases[i] = pauli_string.phase();

                for j in 0..self.n_qubits {
                    let op = pauli_string.get(j);
                    match op {
                        PauliOperator::I => {}
                        PauliOperator::X => self.tableau[i][j] = true,
                        PauliOperator::Z => self.tableau[i][j + self.n_qubits] = true,
                        PauliOperator::Y => {
                            self.tableau[i][j] = true;
                            self.tableau[i][j + self.n_qubits] = true;
                        }
                    }
                }
                return Ok(());
            }
        }

        Err("No space for new stabilizer".to_string())
    }

    /// Create a stabilizer state from a list of stabilizer generators
    pub fn from_stabilizers(stabilizers: Vec<PauliString>) -> Result<Self, String> {
        if stabilizers.is_empty() {
            return Ok(Self::new());
        }

        let n_qubits = stabilizers[0].len();
        let mut state = Self::zeros(n_qubits);

        // Clear the default stabilizers
        for i in 0..2 * n_qubits {
            state.phases[i] = false;
            for j in 0..2 * n_qubits {
                state.tableau[i][j] = false;
            }
        }

        // Add each stabilizer
        for (idx, stabilizer) in stabilizers.iter().enumerate() {
            if idx >= n_qubits {
                return Err("Too many stabilizers provided".to_string());
            }

            let row_idx = n_qubits + idx;
            state.phases[row_idx] = stabilizer.phase();

            for j in 0..n_qubits {
                let op = stabilizer.get(j);
                match op {
                    PauliOperator::I => {}
                    PauliOperator::X => state.tableau[row_idx][j] = true,
                    PauliOperator::Z => state.tableau[row_idx][j + n_qubits] = true,
                    PauliOperator::Y => {
                        state.tableau[row_idx][j] = true;
                        state.tableau[row_idx][j + n_qubits] = true;
                    }
                }
            }
        }

        // Set up corresponding destabilizers
        for i in 0..stabilizers.len() {
            state.tableau[i][i] = true; // X_i destabilizer
        }

        Ok(state)
    }
}

impl QuantumState for StabilizerState {
    fn num_qubits(&self) -> usize {
        self.n_qubits
    }

    fn apply(&mut self, _u: Matrix<Complex>, enumerated: GateType) {
        match enumerated {
            GateType::SingleQubit(single_qubit) => match single_qubit {
                SingleQubitType::PauliX(x) => {
                    let targets = x.construct_targets();
                    self.x(targets[0]);
                }
                SingleQubitType::PauliY(y) => {
                    let targets = y.construct_targets();
                    self.y(targets[0]);
                }
                SingleQubitType::PauliZ(z) => {
                    let targets = z.construct_targets();
                    self.z(targets[0]);
                }
                SingleQubitType::Hadamard(h) => {
                    let targets = h.construct_targets();
                    self.h(targets[0]);
                }
                SingleQubitType::SGate(s) => {
                    let targets = s.construct_targets();
                    self.s(targets[0]);
                }
                _ => unreachable!("Gate not compatible with Clifford stabilizers."),
            },
            GateType::TwoQubit(two_qubit) => match two_qubit {
                TwoQubitType::ControlledNot(cnot) => {
                    let targets = cnot.construct_targets();
                    self.cnot(targets[0], targets[1]);
                }
                TwoQubitType::ControlledPauliZ(cz) => {
                    let targets = cz.construct_targets();
                    self.cz(targets[0], targets[1]);
                }
                _ => unreachable!("Gate not compatible with Clifford stabilizers."),
            },
        }
    }

    fn box_clone(&self) -> Box<dyn QuantumState> {
        Box::new(self.clone())
    }
}

impl BackendOperation for StabilizerState {
    fn expand_state(&mut self, state: Vector<Complex>) {
        // Expand state by one qubit with given 2D vector
        if state.len() != 2 {
            return; // Invalid input
        }

        let old_n_qubits = self.n_qubits;
        self.n_qubits += 1;

        // Resize tableau
        let new_size = 2 * self.n_qubits;
        for row in &mut self.tableau {
            row.resize(new_size, false);
        }
        self.tableau.resize(new_size, vec![false; new_size]);
        self.phases.resize(new_size, false);

        // Determine if we're adding |0⟩ or |1⟩
        let is_one_state = state[1].norm() > state[0].norm();

        // Add new destabilizer (X_new)
        self.tableau[old_n_qubits][old_n_qubits] = true;

        // Add new stabilizer (Z_new)
        self.tableau[self.n_qubits + old_n_qubits][self.n_qubits + old_n_qubits] = true;
        if is_one_state {
            self.phases[self.n_qubits + old_n_qubits] = true; // -Z for |1⟩ state
        }
    }

    /// Perform measurement on qubit following Aaronson-Gottesman algorithm
    fn measure_qubit(&mut self, qubit: usize) -> bool {
        if qubit >= self.n_qubits {
            return false;
        }

        // Check if measurement is deterministic
        // Find if any stabilizer has X_qubit = 1
        let mut p = None;
        for i in self.n_qubits..2 * self.n_qubits {
            if self.tableau[i][qubit] {
                p = Some(i);
                break;
            }
        }

        match p {
            None => {
                // Deterministic measurement
                // Find stabilizer with Z_qubit = 1 to determine outcome
                for i in self.n_qubits..2 * self.n_qubits {
                    if self.tableau[i][qubit + self.n_qubits] {
                        // Measurement outcome is determined by phase
                        // If phase is false (+), measuring |0⟩ gives 0, |1⟩ gives 1
                        // If phase is true (-), measuring |0⟩ gives 1, |1⟩ gives 0
                        return self.phases[i];
                    }
                }
                // If no Z found, default to 0
                false
            }
            Some(p) => {
                // Random measurement - need to update tableau
                let mut rng = rand::rng();
                let result = rng.random_bool(0.5);

                // Step 1: Eliminate X_qubit from all other rows
                for i in 0..2 * self.n_qubits {
                    if i != p && self.tableau[i][qubit] {
                        self.row_mult(i, p);
                    }
                }

                // Step 2: Set row p to correspond to measurement result
                // Clear the row first
                for k in 0..2 * self.n_qubits {
                    self.tableau[p][k] = false;
                }
                // Set Z_qubit = 1
                self.tableau[p][qubit + self.n_qubits] = true;
                // Set phase based on measurement result
                self.phases[p] = result;

                // Step 3: Update corresponding destabilizer (row p - n_qubits)
                let destab_idx = p - self.n_qubits;
                for k in 0..2 * self.n_qubits {
                    self.tableau[destab_idx][k] = false;
                }
                // Set X_qubit = 1 for destabilizer
                self.tableau[destab_idx][qubit] = true;
                self.phases[destab_idx] = false;

                result
            }
        }
    }

    fn measure_all(&mut self) -> Vec<bool> {
        let mut results = Vec::with_capacity(self.n_qubits);
        for qubit in 0..self.n_qubits {
            results.push(self.measure_qubit(qubit));
        }
        results
    }
}

impl QuantumDebugger for StabilizerState {
    type StateType = HashMap<String, Complex>;
    type BitStateType = Complex;
    type QubitStateType = (Complex, Complex);

    fn entire_state(&self, filter_zero: bool) -> Self::StateType {
        todo!()
    }

    fn bit_state(&self, bits: &[bool]) -> Self::BitStateType {
        todo!()
    }

    fn qubit_state(&self, qubit: &QuantumBit) -> Self::QubitStateType {
        todo!()
    }
}
