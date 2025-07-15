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

use std::collections::HashMap;

/// Result after running a quantum circuit
#[derive(Debug, Clone)]
pub struct RunResult {
    pub classical_bit: Vec<Vec<bool>>,
    pub shots: usize,
}

impl RunResult {
    /// Get counts of measurement outcomes (similar to Qiskit's get_counts())
    pub fn get_counts(&self) -> HashMap<String, usize> {
        let mut counts = HashMap::new();

        for measurement in &self.classical_bit {
            // Convert Vec<bool> to binary string (e.g., [true, false] -> "10")
            let bit_string = measurement
                .iter()
                .map(|&bit| if bit { '1' } else { '0' })
                .collect::<String>();

            *counts.entry(bit_string).or_insert(0) += 1;
        }

        counts
    }

    /// Get probabilities of measurement outcomes
    pub fn get_probabilities(&self) -> HashMap<String, f64> {
        let counts = self.get_counts();
        let total_shots = self.shots as f64;

        counts
            .into_iter()
            .map(|(outcome, count)| (outcome, count as f64 / total_shots))
            .collect()
    }

    /// Get the most frequent measurement outcome
    pub fn get_most_frequent(&self) -> Option<(String, usize)> {
        self.get_counts()
            .into_iter()
            .max_by_key(|(_, count)| *count)
    }

    /// Print results in a formatted way (similar to Qiskit's display)
    pub fn print_results(&self) {
        let counts = self.get_counts();
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
    pub fn get_shot(&self, shot_index: usize) -> Option<&Vec<bool>> {
        self.classical_bit.get(shot_index)
    }

    /// Get all raw measurement data
    pub fn get_all_shots(&self) -> &Vec<Vec<bool>> {
        &self.classical_bit
    }

    /// Create a histogram-like display string
    pub fn histogram(&self, max_width: usize) -> String {
        let counts = self.get_counts();
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
