//! Birli is a library of common preprocessing tasks performed in the data pipeline of the Murchison
//! Widefield Array (MWA) Telescope.
//!
//! # Examples
//!
//! Here's an example of how to flag some visibility files
//!
//! ```rust
//! use birli::{
//!     write_flags,
//!     mwalib::CorrelatorContext, write_uvfits,
//!     get_weight_factor, flag_to_weight_array,
//!     FlagContext, VisSelection
//! };
//! use tempfile::tempdir;
//!
//! // define our input files
//! let metafits_path = "tests/data/1297526432_mwax/1297526432.metafits";
//! let gpufits_paths = vec![
//!     "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits",
//!     "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_001.fits",
//!     "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_000.fits",
//!     "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_001.fits",
//! ];
//!
//! // define a temporary directory for output files
//! let tmp_dir = tempdir().unwrap();
//!
//! // define our output paths
//! let flag_template = tmp_dir.path().join("Flagfile%%%.mwaf");
//! let uvfits_out = tmp_dir.path().join("1297526432.uvfits");
//!
//! // Create an mwalib::CorrelatorContext for accessing visibilities.
//! let corr_ctx = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
//!
//! // Determine which timesteps and coarse channels we want to use
//! let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
//!
//! // Prepare our flagmasks with known bad antennae
//! let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
//!
//! // Create a blank array to store flags and visibilities
//! let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
//! let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
//! flag_ctx.set_flags(
//!     &mut flag_array,
//!     &vis_sel.timestep_range,
//!     &vis_sel.coarse_chan_range,
//!     &vis_sel.get_ant_pairs(&corr_ctx.metafits_context)
//! );
//! let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
//!
//! // read visibilities out of the gpubox files
//! vis_sel
//!     .read_mwalib(&corr_ctx, jones_array.view_mut(), flag_array.view_mut(), false)
//!     .unwrap();
//!
//! // write the flags to disk as .mwaf
//! write_flags(&corr_ctx, &flag_array, flag_template.to_str().unwrap(), &vis_sel.coarse_chan_range).unwrap();
//! // write the visibilities to disk as .uvfits
//!
//! let num_pols = corr_ctx.metafits_context.num_visibility_pols;
//! let weight_factor = get_weight_factor(&corr_ctx);
//! let weight_array = flag_to_weight_array(&flag_array.view(), weight_factor);
//! write_uvfits(
//!     uvfits_out.as_path(),
//!     &corr_ctx,
//!     jones_array.view(),
//!     weight_array.view(),
//!     &vis_sel.timestep_range,
//!     &vis_sel.coarse_chan_range,
//!     &vis_sel.baseline_idxs,
//!     None,
//!     None,
//!     1,
//!     1,
//!     false,
//! ).unwrap();
//! ```
//!
//! # Details
//!
//! Birli reads visibilities with [`MWALib`] and uses CXX to bind to the [`AOFlagger`] C++ library.
//! For more details its interface, check out the [`aoflagger::AOFlagger`]
//! documentation
//!
//! [`MWALib`]: https://github.com/MWATelescope/mwalib
//! [`AOFlagger`]: https://gitlab.com/aroffringa/aoflagger
//! [`aoflagger::AOFlagger`]: http://www.andreoffringa.org/aoflagger/doxygen/classaoflagger_1_1AOFlagger.html

#![warn(missing_docs)]
#![warn(rustdoc::missing_doc_code_examples)]
#![deny(clippy::all)]
// //////// //
// pedantic //
// //////// //

// #![warn(clippy::pedantic)]

// useful:
#![warn(clippy::missing_safety_doc)]
#![warn(clippy::missing_errors_doc)]
#![warn(clippy::float_cmp)]
#![warn(clippy::if_not_else)]
#![warn(clippy::explicit_iter_loop)]
#![warn(clippy::wildcard_imports)]
#![warn(clippy::cloned_instead_of_copied)]
#![warn(clippy::cognitive_complexity)]
#![warn(clippy::type_repetition_in_bounds)]
#![warn(clippy::redundant_closure_for_method_calls)]
#![warn(clippy::manual_assert)]
#![warn(clippy::doc_markdown)]
#![warn(clippy::match_same_arms)]
#![warn(clippy::default_trait_access)]
#![warn(clippy::semicolon_if_nothing_returned)]
#![warn(clippy::explicit_into_iter_loop)]
#![warn(clippy::inefficient_to_string)]
#![warn(clippy::needless_pass_by_value)]
#![warn(clippy::used_underscore_binding)]
// TODO: Look at these later:
// #![allow(clippy::too_many_lines)]
// #![allow(clippy::missing_panics_doc)]
// #![allow(clippy::must_use_candidate)]

// won't fix:
// #![allow(clippy::module_name_repetitions)]
// #![allow(clippy::cast_precision_loss)]
// #![allow(clippy::cast_sign_loss)]
// #![allow(clippy::cast_possible_truncation)]
// #![allow(clippy::cast_lossless)]
// #![allow(clippy::unreadable_literal)]
// #![allow(clippy::similar_names)]
// #![allow(clippy::struct_excessive_bools)]
// #![allow(clippy::range_plus_one)]
// #![allow(clippy::cast_possible_wrap)]

// /////// //
// nursery //
// /////// //

// #![warn(clippy::nursery)]

// useful:
#![warn(clippy::use_self)]
#![warn(clippy::missing_const_for_fn)]
#![warn(clippy::option_if_let_else)]
#![warn(clippy::equatable_if_let)]
// TODO: Look at these later:
#![allow(clippy::suboptimal_flops)]
#![allow(clippy::redundant_pub_crate)]
// won't fix:
#![allow(clippy::debug_assert_with_mut_call)]

// ///// //
// cargo //
// ///// //

// #![warn(clippy::cargo)]

use cfg_if::cfg_if;
use log::warn;

pub mod io;
pub use io::{mwaf::FlagFileSet, write_ms, write_uvfits};
pub mod corrections;
pub use corrections::{correct_cable_lengths, correct_geometry, ScrunchType};
pub mod calibration;
pub mod flags;
#[cfg(test)]
pub use approx;
#[cfg(test)]
pub(crate) mod types;
pub use flags::{flag_to_weight_array, get_weight_factor, write_flags, FlagContext};
pub mod passband_gains;
pub use marlu;
pub use marlu::{
    mwalib,
    mwalib::{fitsio, fitsio_sys, CorrelatorContext},
    ndarray,
    ndarray::{parallel::prelude::*, Array3, Axis},
    selection::VisSelection,
    Complex, Jones,
};
#[cfg(test)]
pub(crate) use types::TestJones;

mod error;
pub use error::BirliError;

pub mod preprocessing;
pub use preprocessing::PreprocessContext;

pub mod cli;
pub use cli::BirliContext;

cfg_if! {
    if #[cfg(feature = "aoflagger")] {
        pub use flags::{flag_jones_array, flag_jones_array_existing};
        pub use aoflagger_sys::{cxx_aoflagger_new, CxxAOFlagger, CxxFlagMask, UniquePtr, CxxImageSet};
    }
}

#[cfg(test)]
pub mod test_common;

#[macro_export]
/// Time a statement and increment the timer given by name in the hashmap of durations
macro_rules! with_increment_duration {
    ($durs:expr, $name:literal, $($s:stmt);+ $(;)?) => {
        {
            let _now = std::time::Instant::now();
            let _res = {
                $(
                    $s
                )*
            };
            *$durs.entry($name.into()).or_insert(Duration::default()) += _now.elapsed();
            _res
        }
    };
}
