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

use std::ops::Index;

use crate::bit::{ClassicalBit, QuantumBit};

#[derive(Clone, Debug)]
pub struct QuantumRegister {
    register: Vec<QuantumBit>,
    capacity: usize,
    label: Option<String>,
}

impl QuantumRegister {
    pub fn new(capacity: usize) -> Self {
        let mut register = Vec::new();
        for i in 0..capacity {
            register.push(QuantumBit::new(i));
        }
        Self {
            register,
            capacity,
            label: None,
        }
    }

    pub fn label(&self, label: &str) -> Self {
        let mut new_self = self.clone();
        new_self.label = Some(label.to_string());
        new_self
    }

    pub fn push(&mut self, qubit: QuantumBit) {
        self.register.push(qubit);
        self.capacity += 1;
    }

    pub fn push_iter<I>(&mut self, qubits: I)
    where
        I: IntoIterator<Item = QuantumBit>,
    {
        let qubit_vec: Vec<QuantumBit> = qubits.into_iter().collect();
        self.capacity += qubit_vec.len();
        self.register.extend(qubit_vec);
    }

    pub fn num_qubits(&self) -> usize {
        self.capacity
    }

    pub fn extend(&mut self, qreg: QuantumRegister) {
        self.register.extend(qreg.register);
        self.capacity += qreg.capacity;
    }
}

impl Index<usize> for QuantumRegister {
    type Output = QuantumBit;

    fn index(&self, index: usize) -> &Self::Output {
        &self.register[index]
    }
}

impl IntoIterator for QuantumRegister {
    type Item = QuantumBit;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.register.into_iter()
    }
}

#[derive(Clone, Debug)]
pub struct ClassicalRegister {
    register: Vec<ClassicalBit>,
    capacity: usize,
    label: Option<String>,
}

impl ClassicalRegister {
    pub fn new(capacity: usize) -> Self {
        let mut register = Vec::new();
        for i in 0..capacity {
            register.push(ClassicalBit::new(i));
        }
        Self {
            register,
            capacity,
            label: None,
        }
    }

    pub fn label(&self, label: &str) -> Self {
        let mut new_self = self.clone();
        new_self.label = Some(label.to_string());
        new_self
    }

    pub fn get_raw(&self, index: usize) -> Option<bool> {
        self.register[index].raw_outcome()
    }

    pub fn get(&self, index: usize) -> bool {
        self.register[index].outcome()
    }

    pub fn reg(&self, bit: &ClassicalBit) -> usize {
        self.get(bit.index()) as usize
    }

    pub fn set(&mut self, index: usize, value: bool) {
        self.register[index].set_outcome(value)
    }

    pub fn extend(&mut self, creg: ClassicalRegister) {
        self.register.extend(creg.register);
        self.capacity += creg.capacity;
    }

    pub fn capacity(&self) -> usize {
        self.register.len()
    }
}

impl Index<usize> for ClassicalRegister {
    type Output = ClassicalBit;

    fn index(&self, index: usize) -> &Self::Output {
        &self.register[index]
    }
}

impl IntoIterator for ClassicalRegister {
    type Item = ClassicalBit;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.register.into_iter()
    }
}
