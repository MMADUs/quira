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

//! Quira ~ Strong Type System Quantum Computing Package.
//! Start building your quantum algorithms!
//!
//! # Core Modules
//! * `include` - The main quantum module.
//! * `operation` - The quantum operation module.
//! * `common` - The common quantum type module.
//! * `provider` - The extra quantum provider module.

mod bit;
mod circuit;
mod constant;
mod endian;
mod io;
mod job;
mod kernel;
mod operations;
mod ops;
mod prebuilt;
mod prelude;
mod register;
mod result;
mod simulator;
mod types;
mod utils;

/// The main quantum module
pub mod include {
    pub use crate::circuit::*;
    pub use crate::job::*;
    pub use crate::kernel::*;
    pub use crate::register::*;
    pub use crate::simulator::*;
    pub use crate::bit::*;
}

/// The quantum operation module
pub mod operation {
    pub use crate::operations::*;
}

/// The common type module
pub mod common {
    pub use crate::constant::*;
    pub use crate::prelude::*;
    pub use crate::result::*;
    pub use crate::types::*;
}

/// The extra provider quantum module
pub mod provider {
    pub use crate::endian::*;
    pub use crate::io::*;
    pub use crate::ops::*;
    pub use crate::prebuilt::*;
    pub use crate::utils::*;
}
