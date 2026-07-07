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

use std::collections::HashMap;
use std::fmt;

/// Represents a single Pauli operator
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PauliOperator {
    I, // Identity
    X, // Pauli-X
    Y, // Pauli-Y
    Z, // Pauli-Z
}

impl fmt::Display for PauliOperator {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PauliOperator::I => write!(f, "I"),
            PauliOperator::X => write!(f, "X"),
            PauliOperator::Y => write!(f, "Y"),
            PauliOperator::Z => write!(f, "Z"),
        }
    }
}

impl PauliOperator {
    /// Create from character
    pub fn from_char(c: char) -> Option<Self> {
        match c.to_ascii_uppercase() {
            'I' => Some(PauliOperator::I),
            'X' => Some(PauliOperator::X),
            'Y' => Some(PauliOperator::Y),
            'Z' => Some(PauliOperator::Z),
            _ => None,
        }
    }

    /// Check if operator is identity
    pub fn is_identity(&self) -> bool {
        matches!(self, PauliOperator::I)
    }

    /// Commutation relation with another Pauli operator
    pub fn commutes_with(&self, other: &PauliOperator) -> bool {
        match (self, other) {
            (PauliOperator::I, _) | (_, PauliOperator::I) => true,
            (a, b) if a == b => true,
            _ => false,
        }
    }

    /// Anticommutation relation
    pub fn anticommutes_with(&self, other: &PauliOperator) -> bool {
        !self.commutes_with(other) && !self.is_identity() && !other.is_identity()
    }
}

/// Represents a tensor product of Pauli operators with phase
#[derive(Clone, Debug)]
pub struct PauliString {
    /// Map from qubit index to Pauli operator
    pub(crate) operators: HashMap<usize, PauliOperator>,
    /// Phase: false = +1, true = -1
    phase: bool,
    /// Number of qubits this string can act on
    n_qubits: usize,
}

impl PauliString {
    /// Create empty Pauli string (identity)
    pub fn new(n_qubits: usize) -> Self {
        Self {
            operators: HashMap::new(),
            phase: false,
            n_qubits,
        }
    }

    /// Create from string representation like "XIZY"
    pub fn from_str(s: &str) -> Result<Self, String> {
        let mut operators = HashMap::new();

        for (i, c) in s.chars().enumerate() {
            if let Some(op) = PauliOperator::from_char(c) {
                if !op.is_identity() {
                    operators.insert(i, op);
                }
            } else {
                return Err(format!("Invalid Pauli operator: {}", c));
            }
        }

        Ok(Self {
            operators,
            phase: false,
            n_qubits: s.len(),
        })
    }

    /// Create from vector of operators
    pub fn from_ops(ops: Vec<PauliOperator>) -> Self {
        let mut operators = HashMap::new();

        for (i, op) in ops.into_iter().enumerate() {
            if !op.is_identity() {
                operators.insert(i, op);
            }
        }

        let len = operators.len();

        Self {
            operators,
            phase: false,
            n_qubits: len,
        }
    }

    /// Set operator at specific qubit
    pub fn set(&mut self, qubit: usize, op: PauliOperator) -> Result<(), String> {
        if qubit >= self.n_qubits {
            return Err(format!("Qubit index {} out of bounds", qubit));
        }

        if op.is_identity() {
            self.operators.remove(&qubit);
        } else {
            self.operators.insert(qubit, op);
        }
        Ok(())
    }

    /// Get operator at specific qubit
    pub fn get(&self, qubit: usize) -> PauliOperator {
        self.operators
            .get(&qubit)
            .copied()
            .unwrap_or(PauliOperator::I)
    }

    /// Set phase
    pub fn set_phase(&mut self, negative: bool) {
        self.phase = negative;
    }

    /// Get phase
    pub fn phase(&self) -> bool {
        self.phase
    }

    /// Flip phase
    pub fn flip_phase(&mut self) {
        self.phase = !self.phase;
    }

    /// Get number of qubits
    pub fn len(&self) -> usize {
        self.n_qubits
    }

    /// Check if string is empty (all identity)
    pub fn is_empty(&self) -> bool {
        self.operators.is_empty()
    }

    /// Get weight (number of non-identity operators)
    pub fn weight(&self) -> usize {
        self.operators.len()
    }

    /// Check if this string commutes with another
    pub fn commutes_with(&self, other: &PauliString) -> bool {
        let mut anticommute_count = 0;

        for (&qubit, &op1) in &self.operators {
            if let Some(&op2) = other.operators.get(&qubit) {
                if op1.anticommutes_with(&op2) {
                    anticommute_count += 1;
                }
            }
        }

        anticommute_count % 2 == 0
    }

    /// Multiply with another Pauli string
    pub fn multiply(&self, other: &PauliString) -> Result<PauliString, String> {
        if self.n_qubits != other.n_qubits {
            return Err("Cannot multiply Pauli strings of different lengths".to_string());
        }

        let mut result = PauliString::new(self.n_qubits);
        result.phase = self.phase ^ other.phase;

        // Collect all qubits that have operators in either string
        let mut all_qubits: Vec<usize> = self
            .operators
            .keys()
            .chain(other.operators.keys())
            .copied()
            .collect();
        all_qubits.sort();
        all_qubits.dedup();

        for qubit in all_qubits {
            let op1 = self.get(qubit);
            let op2 = other.get(qubit);

            let (new_op, phase_change) = multiply_pauli_ops(op1, op2);

            if phase_change {
                result.flip_phase();
            }

            if !new_op.is_identity() {
                result.operators.insert(qubit, new_op);
            }
        }

        Ok(result)
    }

    /// Create single-qubit Pauli strings
    pub fn x(qubit: usize, n_qubits: usize) -> Self {
        let mut string = Self::new(n_qubits);
        string.set(qubit, PauliOperator::X).unwrap();
        string
    }

    pub fn y(qubit: usize, n_qubits: usize) -> Self {
        let mut string = Self::new(n_qubits);
        string.set(qubit, PauliOperator::Y).unwrap();
        string
    }

    pub fn z(qubit: usize, n_qubits: usize) -> Self {
        let mut string = Self::new(n_qubits);
        string.set(qubit, PauliOperator::Z).unwrap();
        string
    }

    /// Convert to string representation
    pub fn to_string(&self) -> String {
        let mut result = String::new();

        if self.phase {
            result.push('-');
        }

        for i in 0..self.n_qubits {
            result.push_str(&format!("{}", self.get(i)));
        }

        result
    }
}

impl fmt::Display for PauliString {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

/// Multiply two Pauli operators, returning (result, phase_flip)
fn multiply_pauli_ops(op1: PauliOperator, op2: PauliOperator) -> (PauliOperator, bool) {
    use PauliOperator::*;

    match (op1, op2) {
        (I, op) | (op, I) => (op, false),
        (X, X) | (Y, Y) | (Z, Z) => (I, false),
        (X, Y) => (Z, false),
        (Y, X) => (Z, true),
        (Y, Z) => (X, false),
        (Z, Y) => (X, true),
        (Z, X) => (Y, false),
        (X, Z) => (Y, true),
    }
}
