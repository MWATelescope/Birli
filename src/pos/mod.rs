// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! Positional coordinate types.
//!
//! Most of this was blatently stolen (with permission) from [Chris Jordan](https://github.com/cjordan)

pub mod azel;
pub mod enh;
pub mod hadec;
pub mod lmn;
pub mod radec;
pub mod uvw;
pub mod xyz;
pub mod pal;
pub mod precess;
pub mod earth;
pub mod error;

// Re-exports.
pub use azel::AzEl;
pub use enh::ENH;
pub use hadec::HADec;
pub use lmn::LMN;
pub use radec::RADec;
pub use uvw::UVW;
pub use xyz::{XyzBaseline, XyzGeocentric, XyzGeodetic};
