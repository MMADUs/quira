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

use ndarray::{Array1, Array2, Array3};
use num_complex::Complex64;

/// Index type for qubit
pub type Qubit = usize;

/// Mathematical type for complex number
pub type Complex = Complex64;

/// Mathematical structure for 1d tensor (Vector)
pub type Vector<T> = Array1<T>;

/// Mathematical structure for 2d tensor (Matrix)
pub type Matrix<T> = Array2<T>;

/// Mathematical structure for 3d tensor (Tensor)
pub type Tensor<T> = Array3<T>;
