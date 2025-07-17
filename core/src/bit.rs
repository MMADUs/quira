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

use core::fmt;

/// Index type for qubit
#[derive(Clone, Debug)]
pub struct QuantumBit {
    index: usize,
    label: Option<String>,
}

impl QuantumBit {
    pub fn new(index: usize) -> Self {
        Self { index, label: None }
    }

    pub fn label(&self, label: &str) -> Self {
        if let Some(current) = &self.label {
            println!("Unable to rename Qubit {} to {}", current, label);
            return self.clone();
        }

        let mut new_self = self.clone();
        new_self.label = Some(label.to_string());
        new_self
    }

    pub fn index(&self) -> usize {
        self.index
    }
}

impl fmt::Display for QuantumBit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Qubit={}, Label={}",
            self.index,
            self.label.as_deref().unwrap_or("_")
        )
    }
}

/// Type of classical bit
#[derive(Clone, Debug)]
pub struct ClassicalBit {
    index: usize,
    outcome: Option<bool>,
    label: Option<String>,
}

impl ClassicalBit {
    pub fn new(index: usize) -> Self {
        Self {
            index,
            outcome: None,
            label: None,
        }
    }

    pub fn label(&self, label: &str) -> Self {
        let mut new_self = self.clone();
        new_self.label = Some(label.to_string());
        new_self
    }

    pub fn set_outcome(&mut self, bit: bool) {
        self.outcome = Some(bit)
    }

    pub fn raw_outcome(&self) -> Option<bool> {
        self.outcome
    }

    pub fn outcome(&self) -> bool {
        self.raw_outcome().expect(&format!(
            "No outcome found in classical bit: {}",
            self.index()
        ))
    }

    pub fn bit(&self) -> usize {
        self.outcome() as usize
    }

    pub fn index(&self) -> usize {
        self.index
    }
}

impl fmt::Display for ClassicalBit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let outcome_str = match self.outcome {
            Some(true) => "1",
            Some(false) => "0",
            None => "_",
        };
        let label_str = self.label.as_deref().unwrap_or("_");
        write!(f, "Bit={}, Label={}", outcome_str, label_str)
    }
}
