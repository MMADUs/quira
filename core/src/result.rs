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

use rayon::prelude::*;
use std::collections::HashMap;

use crate::register::ClassicalRegister;

/// Result after running a quantum circuit
#[derive(Debug, Clone)]
pub struct RunResult {
    pub classical_bit: Vec<ClassicalRegister>,
    pub shots: usize,
}

impl RunResult {
    /// Get counts of measurement outcomes (similar to Qiskit's get_counts())
    /// Returns empty HashMap if no measurements were performed
    pub fn get_counts(&self) -> HashMap<String, usize> {
        // Check if any measurements were actually performed
        if !self.has_measurements() {
            return HashMap::new();
        }

        // Use parallel processing for large shot counts
        if self.classical_bit.len() > 1000 {
            // Parallel reduce for large datasets
            self.classical_bit
                .par_iter()
                .map(|measurement| {
                    let bit_string = self.classical_register_to_string(measurement);
                    let mut local_counts = HashMap::new();
                    *local_counts.entry(bit_string).or_insert(0) += 1;
                    local_counts
                })
                .reduce(
                    || HashMap::new(),
                    |mut acc, local_counts| {
                        for (key, count) in local_counts {
                            *acc.entry(key).or_insert(0) += count;
                        }
                        acc
                    },
                )
        } else {
            // Sequential processing for smaller datasets to avoid overhead
            let mut counts = HashMap::new();
            for measurement in &self.classical_bit {
                let bit_string = self.classical_register_to_string(measurement);
                *counts.entry(bit_string).or_insert(0) += 1;
            }
            counts
        }
    }

    /// Check if any measurements were performed in the circuit
    /// This checks if any classical bits have been set (have outcomes)
    fn has_measurements(&self) -> bool {
        if self.classical_bit.is_empty() {
            return false;
        }

        // Use parallel processing for large shot counts
        if self.classical_bit.len() > 100 {
            self.classical_bit
                .par_iter()
                .any(|register| self.register_has_measurements(register))
        } else {
            // Sequential for smaller datasets
            self.classical_bit
                .iter()
                .any(|register| self.register_has_measurements(register))
        }
    }

    /// Check if a classical register has any measured bits
    fn register_has_measurements(&self, register: &ClassicalRegister) -> bool {
        for i in 0..register.capacity() {
            // Use get_raw to check if the bit has been measured (not None)
            if register.get_raw(i).is_some() {
                return true;
            }
        }
        false
    }

    /// Helper function to convert ClassicalRegister to binary string
    /// Only includes bits that were actually measured
    fn classical_register_to_string(&self, register: &ClassicalRegister) -> String {
        let mut bit_string = String::new();
        for i in 0..register.capacity() {
            if let Some(bit_value) = register.get_raw(i) {
                bit_string.push(if bit_value { '1' } else { '0' });
            }
            // Skip bits that are None (not measured)
        }
        bit_string
    }

    /// Get probabilities of measurement outcomes
    /// Returns empty HashMap if no measurements were performed
    pub fn get_probabilities(&self) -> HashMap<String, f64> {
        let counts = self.get_counts();

        if counts.is_empty() {
            return HashMap::new();
        }

        let total_shots = self.shots as f64;

        // Use parallel processing for large result sets
        if counts.len() > 50 {
            counts
                .into_par_iter()
                .map(|(outcome, count)| (outcome, count as f64 / total_shots))
                .collect()
        } else {
            counts
                .into_iter()
                .map(|(outcome, count)| (outcome, count as f64 / total_shots))
                .collect()
        }
    }

    /// Get the most frequent measurement outcome
    /// Returns None if no measurements were performed
    pub fn get_most_frequent(&self) -> Option<(String, usize)> {
        self.get_counts()
            .into_iter()
            .max_by_key(|(_, count)| *count)
    }

    /// Print results in a formatted way (similar to Qiskit's display)
    pub fn print_results(&self) {
        let counts = self.get_counts();

        if counts.is_empty() {
            println!(
                "No measurements performed in this circuit ({} shots)",
                self.shots
            );
            return;
        }

        let probabilities = self.get_probabilities();

        println!("Measurement Results ({} shots):", self.shots);
        println!("{:-<50}", "");
        println!("{:<10} {:<10} {:<15}", "Outcome", "Count", "Probability");
        println!("{:-<50}", "");

        // Sort by outcome for consistent display
        let mut sorted_outcomes: Vec<_> = counts.iter().collect();
        sorted_outcomes.sort_by(|a, b| a.0.cmp(b.0));

        for (outcome, count) in sorted_outcomes {
            let probability = probabilities[outcome];
            println!("{:<10} {:<10} {:<15.4}", outcome, count, probability);
        }
        println!("{:-<50}", "");
    }

    /// Get raw measurement data for a specific shot
    pub fn get_shot(&self, shot_index: usize) -> Option<&ClassicalRegister> {
        self.classical_bit.get(shot_index)
    }

    /// Get all raw measurement data
    pub fn get_all_shots(&self) -> &Vec<ClassicalRegister> {
        &self.classical_bit
    }

    /// Create a histogram-like display string
    pub fn histogram(&self, max_width: usize) -> String {
        let counts = self.get_counts();

        if counts.is_empty() {
            return "No measurements performed".to_string();
        }

        let max_count = counts.values().max().copied().unwrap_or(0);

        if max_count == 0 {
            return "No measurements found".to_string();
        }

        let mut result = String::new();
        let mut sorted_outcomes: Vec<_> = counts.iter().collect();
        sorted_outcomes.sort_by(|a, b| a.0.cmp(b.0));

        for (outcome, count) in sorted_outcomes {
            let bar_length = (count * max_width) / max_count.max(1);
            let bar = "█".repeat(bar_length);
            result.push_str(&format!("{}: {} ({})\n", outcome, bar, count));
        }

        result
    }
}
