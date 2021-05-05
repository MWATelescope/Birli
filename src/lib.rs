#![warn(missing_docs)]
#![warn(missing_doc_code_examples)]

//! Birli is a library of common preprocessing tasks performed in the data pipeline of the Murchison
//! Widefield Array (MWA) Telescope.
//!
//! # Examples
//!
//! Here's an example of how to flag some visibility files
//!
//! ```rust
//! use birli::{context_to_baseline_imgsets, flag_imgsets, write_flags, cxx_aoflagger_new};
//! use mwalib::CorrelatorContext;
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
//! // define our output flag file template
//! let flag_template = tmp_dir.path().join("Flagfile%%%.mwaf");
//!
//! // Create an mwalib::CorrelatorContext for accessing visibilities.
//! let context = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
//!
//! // create a CxxAOFlagger object to perform AOFlagger operations
//! let aoflagger = unsafe { cxx_aoflagger_new() };
//!
//! // generate imagesets for each baseline in the format required by aoflagger
//! let baseline_imgsets = context_to_baseline_imgsets(&aoflagger, &context);
//!
//! // use the default strategy file location for MWA
//! let strategy_filename = &aoflagger.FindStrategyFileMWA();
//!
//! // run the strategy on the imagesets, and get the resulting flagmasks for each baseline
//! let baseline_flagmasks = flag_imgsets(&aoflagger, &strategy_filename, baseline_imgsets);
//!
//! // Get a list of all gpubox IDs
//! let gpubox_ids: Vec<usize> = context
//!             .coarse_chans
//!             .iter()
//!             .map(|chan| chan.gpubox_number)
//!             .collect();
//!
//! // write the flags to disk as .mwaf
//! write_flags(&context, baseline_flagmasks, flag_template.to_str().unwrap(), &gpubox_ids);
//! ```
//!
//! # Details
//!
//! Birli reads visibilities with [`MWALib`] and uses CXX to bind to the [`AOFlagger`] C++ library.
//! For more details about AOFlagger's interface, check out the [`aoflagger::AOFlagger`]
//! documentation
//!
//! [`MWALib`]: https://github.com/MWATelescope/mwalib
//! [`AOFlagger`]: https://gitlab.com/aroffringa/aoflagger
//! [`aoflagger::AOFlagger`]: http://www.andreoffringa.org/aoflagger/doxygen/classaoflagger_1_1AOFlagger.html

mod cxx_aoflagger;
use cxx::UniquePtr;
pub use cxx_aoflagger::ffi::{
    cxx_aoflagger_new, CxxAOFlagger, CxxFlagMask, CxxImageSet, CxxStrategy,
};

use mwalib::CorrelatorContext;
use std::collections::BTreeMap;
use std::os::raw::c_short;

pub mod flag_io;
use flag_io::FlagFileSet;
mod error;

/// Get the version of the AOFlagger library from the library itself.
///
/// # Examples
///
/// ```rust
/// use birli::get_aoflagger_version_string;
/// use regex::Regex;
///
/// let aoflagger_version = get_aoflagger_version_string();
/// // This ensures we're using aoflagger 3.*
/// let version_regex = Regex::new(r"3\.\d+\.\d+").unwrap();
/// assert!(version_regex.is_match(&aoflagger_version));
/// ```
pub fn get_aoflagger_version_string() -> String {
    let mut major: c_short = -1;
    let mut minor: c_short = -1;
    let mut sub_minor: c_short = -1;

    unsafe {
        let aoflagger = cxx_aoflagger_new();
        aoflagger.GetVersion(&mut major, &mut minor, &mut sub_minor);
    }

    return format!("{}.{}.{}", major, minor, sub_minor);
}

/// Read aan observation's visibilities into a [`std::collections::BTreeMap`] mapping each baseline
/// in the observation to a [`CxxImageSet`], given a [`CxxAOFlagger`] instance and that
/// observation's [`mwalib::CorrelatorContext`].
///
/// [`mwalib::CorrelatorContext`]: https://docs.rs/mwalib/0.7.0/mwalib/struct.CorrelatorContext.html
///
/// # Examples
///
/// ```
/// use birli::{context_to_baseline_imgsets, flag_imgsets, write_flags, cxx_aoflagger_new};
/// use mwalib::CorrelatorContext;
/// use tempfile::tempdir;
///
/// // define our input files
/// let metafits_path = "tests/data/1297526432_mwax/1297526432.metafits";
/// let gpufits_paths = vec![
///     "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits",
///     "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_001.fits",
///     "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_000.fits",
///     "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_001.fits",
/// ];
///
/// // Create an mwalib::CorrelatorContext for accessing visibilities.
/// let context = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
///
/// let baseline_imgsets = unsafe {
///     let aoflagger = cxx_aoflagger_new();
///     context_to_baseline_imgsets(&aoflagger, &context)
/// };
/// ```
pub fn context_to_baseline_imgsets(
    aoflagger: &CxxAOFlagger,
    context: &CorrelatorContext,
) -> BTreeMap<usize, UniquePtr<CxxImageSet>> {
    let coarse_chan_arr = context.coarse_chans.clone();
    let timestep_arr = context.timesteps.clone();

    let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
    let floats_per_finechan = context.metafits_context.num_visibility_pols * 2;
    let floats_per_baseline = fine_chans_per_coarse * floats_per_finechan;
    let height = context.num_coarse_chans * fine_chans_per_coarse;
    let width = context.num_timesteps;
    let img_stride = (((width - 1) / 8) + 1) * 8;

    let mut baseline_imgsets: BTreeMap<usize, UniquePtr<CxxImageSet>> = context
        .metafits_context
        .baselines
        .iter()
        .enumerate()
        .map(|(baseline_idx, _)| unsafe {
            (
                baseline_idx,
                aoflagger.MakeImageSet(width, height, 8, 0 as f32, width),
            )
        })
        .collect();

    for (coarse_chan_idx, _) in coarse_chan_arr.iter().enumerate() {
        for (timestep_idx, _) in timestep_arr.iter().enumerate() {
            let img_buf = context
                .read_by_baseline(timestep_idx, coarse_chan_idx)
                .unwrap();
            for (baseline_idx, baseline_chunk) in img_buf.chunks(floats_per_baseline).enumerate() {
                let imgset = baseline_imgsets.get_mut(&baseline_idx).unwrap();

                for float_idx in 0..8 {
                    let imgset_buf = imgset.pin_mut().ImageBufferMut(float_idx);
                    for (fine_chan_idx, fine_chan_chunk) in
                        baseline_chunk.chunks(floats_per_finechan).enumerate()
                    {
                        let x = timestep_idx;
                        let y = fine_chans_per_coarse * coarse_chan_idx + fine_chan_idx;
                        imgset_buf[y * img_stride + x] = fine_chan_chunk[float_idx];
                    }
                }
            }
        }
    }

    return baseline_imgsets;
}

/// Flag an observation's visibilities, given a [`CxxAOFlagger`] instance, a [`CxxStrategy`]
/// filename, and a [`std::collections::BTreeMap`] mapping each baseline in an observation to a
/// [`CxxImageSet`], returning a [`std::collections::BTreeMap`] mapping from each baseline in the
/// observation to a [`CxxFlagMask`].
///
/// # Examples
///
/// ```
/// use birli::{context_to_baseline_imgsets, flag_imgsets, write_flags, cxx_aoflagger_new};
/// use mwalib::CorrelatorContext;
/// use tempfile::tempdir;
///
/// // define our input files
/// let metafits_path = "tests/data/1297526432_mwax/1297526432.metafits";
/// let gpufits_paths = vec![
///     "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits",
///     "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_001.fits",
///     "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_000.fits",
///     "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_001.fits",
/// ];
///
/// // Create an mwalib::CorrelatorContext for accessing visibilities.
/// let context = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
///
/// let baseline_imgsets = unsafe {
///     let aoflagger = cxx_aoflagger_new();
///     context_to_baseline_imgsets(&aoflagger, &context)
/// };
/// ```
pub fn flag_imgsets(
    aoflagger: &CxxAOFlagger,
    strategy_filename: &String,
    baseline_imgsets: BTreeMap<usize, UniquePtr<CxxImageSet>>,
) -> BTreeMap<usize, UniquePtr<CxxFlagMask>> {
    // TODO: figure out how to parallelize with Rayon, into_iter(). You'll probably need to convert between UniquePtr and Box

    return baseline_imgsets
        .iter()
        .map(|(&baseline, imgset)| {
            (
                baseline,
                aoflagger.LoadStrategyFile(strategy_filename).Run(&imgset),
            )
        })
        .collect();
}

/// Write flags to disk, given an observation's [`mwalib::CorrelatorContext`], a
/// [`std::collections::BTreeMap`] mapping from each baseline in the observation to a
/// [`CxxFlagMask`], a filename template and a vector of gpubox IDs, .
///
/// The filename template should contain two or 3 percentage (`%`) characters which will be replaced
/// by the gpubox id or channel number (depending on correlator type). See [`flag_io::FlagFileSet::new`]
///
/// # Examples
///
/// Here's an example of how to flag some visibility files
///
/// ```rust
/// use birli::{context_to_baseline_imgsets, flag_imgsets, write_flags, cxx_aoflagger_new};
/// use mwalib::CorrelatorContext;
/// use tempfile::tempdir;
///
/// // define our input files
/// let metafits_path = "tests/data/1297526432_mwax/1297526432.metafits";
/// let gpufits_paths = vec![
///     "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits",
///     "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_001.fits",
///     "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_000.fits",
///     "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_001.fits",
/// ];
///
/// // define a temporary directory for output files
/// let tmp_dir = tempdir().unwrap();
///
/// // define our output flag file template
/// let flag_template = tmp_dir.path().join("Flagfile%%%.mwaf");
///
/// // Create an mwalib::CorrelatorContext for accessing visibilities.
/// let context = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
///
/// // create a CxxAOFlagger object to perform AOFlagger operations
/// let aoflagger = unsafe { cxx_aoflagger_new() };
///
/// // generate imagesets for each baseline in the format required by aoflagger
/// let baseline_imgsets = context_to_baseline_imgsets(&aoflagger, &context);
///
/// // use the default strategy file location for MWA
/// let strategy_filename = &aoflagger.FindStrategyFileMWA();
///
/// // run the strategy on the imagesets, and get the resulting flagmasks for each baseline
/// let baseline_flagmasks = flag_imgsets(&aoflagger, &strategy_filename, baseline_imgsets);
///
/// // Get a list of all gpubox IDs
/// let gpubox_ids: Vec<usize> = context
///             .coarse_chans
///             .iter()
///             .map(|chan| chan.gpubox_number)
///             .collect();
///
/// // write the flags to disk as .mwaf
/// write_flags(&context, baseline_flagmasks, flag_template.to_str().unwrap(), &gpubox_ids);
/// ```
pub fn write_flags(
    context: &CorrelatorContext,
    baseline_flagmasks: BTreeMap<usize, UniquePtr<CxxFlagMask>>,
    filename_template: &str,
    gpubox_ids: &Vec<usize>,
) {
    let mut flag_file_set = FlagFileSet::new(context, filename_template, &gpubox_ids).unwrap();
    flag_file_set
        .write_baseline_flagmasks(&context, baseline_flagmasks)
        .unwrap();
}

#[cfg(test)]
mod tests {
    use super::{context_to_baseline_imgsets, flag_imgsets, write_flags, FlagFileSet};
    use crate::cxx_aoflagger::ffi::{cxx_aoflagger_new, CxxFlagMask, CxxImageSet};
    use cxx::UniquePtr;
    use glob::glob;
    use mwalib::CorrelatorContext;
    use std::collections::BTreeMap;
    use tempfile::tempdir;

    fn get_mwax_context() -> CorrelatorContext {
        let metafits_path = "tests/data/1297526432_mwax/1297526432.metafits";
        let gpufits_paths = vec![
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits",
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_001.fits",
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_000.fits",
            "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_001.fits",
        ];
        CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap()
    }

    fn get_mwa_ord_context() -> CorrelatorContext {
        let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
        let gpufits_paths = vec![
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits",
        ];
        CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap()
    }

    #[test]
    fn test_context_to_baseline_imgsets_mwax() {
        let context = get_mwax_context();
        let width = context.num_timesteps;
        let img_stride = (((width - 1) / 8) + 1) * 8;

        let baseline_imgsets = unsafe {
            let aoflagger = cxx_aoflagger_new();
            context_to_baseline_imgsets(&aoflagger, &context)
        };

        let imgset0 = baseline_imgsets.get(&0).unwrap();

        assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 0], 0x410000 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 1], 0x410100 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 2], 0x410200 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 3], 0x410300 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 4], 0.0);
        assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 5], 0.0);
        assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 6], 0.0);
        assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 7], 0.0);

        assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 0], 0x410008 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 1], 0x410108 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 2], 0x410208 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 3], 0x410308 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 4], 0.0);
        assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 5], 0.0);
        assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 6], 0.0);
        assert_eq!(imgset0.ImageBuffer(0)[1 * img_stride + 7], 0.0);

        assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 0], 0x410400 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 1], 0x410500 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 2], 0x410600 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 3], 0x410700 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 4], 0.0);
        assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 5], 0.0);
        assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 6], 0.0);
        assert_eq!(imgset0.ImageBuffer(0)[2 * img_stride + 7], 0.0);

        assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 0], 0x410408 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 1], 0x410508 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 2], 0x410608 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 3], 0x410708 as f32);
        assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 4], 0.0);
        assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 5], 0.0);
        assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 6], 0.0);
        assert_eq!(imgset0.ImageBuffer(0)[3 * img_stride + 7], 0.0);

        assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 0], 0x410001 as f32);
        assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 1], 0x410101 as f32);
        assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 2], 0x410201 as f32);
        assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 3], 0x410301 as f32);
        assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 4], 0.0);
        assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 5], 0.0);
        assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 6], 0.0);
        assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 7], 0.0);

        /* ... */
        assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 0], 0x410007 as f32);
        assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 1], 0x410107 as f32);
        assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 2], 0x410207 as f32);
        assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 3], 0x410307 as f32);
        assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 4], 0.0);
        assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 5], 0.0);
        assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 6], 0.0);
        assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 7], 0.0);

        let imgset2 = baseline_imgsets.get(&2).unwrap();

        assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 0], 0x410020 as f32);
        assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 1], 0x410120 as f32);
        assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 2], 0x410220 as f32);
        assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 3], 0x410320 as f32);
        assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 4], 0.0);
        assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 5], 0.0);
        assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 6], 0.0);
        assert_eq!(imgset2.ImageBuffer(0)[0 * img_stride + 7], 0.0);

        assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 0], 0x410028 as f32);
        assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 1], 0x410128 as f32);
        assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 2], 0x410228 as f32);
        assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 3], 0x410328 as f32);
        assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 4], 0.0);
        assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 5], 0.0);
        assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 6], 0.0);
        assert_eq!(imgset2.ImageBuffer(0)[1 * img_stride + 7], 0.0);
    }

    #[test]
    fn test_context_to_baseline_imgsets_mwa_ord() {
        let context = get_mwa_ord_context();
        let width = context.num_timesteps;
        let img_stride = (((width - 1) / 8) + 1) * 8;

        let baseline_imgsets = unsafe {
            let aoflagger = cxx_aoflagger_new();
            context_to_baseline_imgsets(&aoflagger, &context)
        };

        let imgset0 = baseline_imgsets.get(&0).unwrap();

        assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 0], 0x10c5be as f32);
        assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 1], 0x14c5be as f32);
        assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 2], 0x18c5be as f32);
        assert_eq!(imgset0.ImageBuffer(0)[0 * img_stride + 3], 0x1cc5be as f32);

        assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 0], 0x10c5bf as f32);
        assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 1], 0x14c5bf as f32);
        assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 2], 0x18c5bf as f32);
        assert_eq!(imgset0.ImageBuffer(1)[0 * img_stride + 3], 0x1cc5bf as f32);

        assert_eq!(imgset0.ImageBuffer(2)[0 * img_stride + 0], 0x10c5ae as f32);
        assert_eq!(imgset0.ImageBuffer(2)[0 * img_stride + 1], 0x14c5ae as f32);
        assert_eq!(imgset0.ImageBuffer(2)[0 * img_stride + 2], 0x18c5ae as f32);
        assert_eq!(imgset0.ImageBuffer(2)[0 * img_stride + 3], 0x1cc5ae as f32);

        assert_eq!(imgset0.ImageBuffer(3)[0 * img_stride + 0], -0x10c5af as f32);
        assert_eq!(imgset0.ImageBuffer(3)[0 * img_stride + 1], -0x14c5af as f32);
        assert_eq!(imgset0.ImageBuffer(3)[0 * img_stride + 2], -0x18c5af as f32);
        assert_eq!(imgset0.ImageBuffer(3)[0 * img_stride + 3], -0x1cc5af as f32);

        assert_eq!(imgset0.ImageBuffer(4)[0 * img_stride + 0], 0x10c5ae as f32);
        assert_eq!(imgset0.ImageBuffer(4)[0 * img_stride + 1], 0x14c5ae as f32);
        assert_eq!(imgset0.ImageBuffer(4)[0 * img_stride + 2], 0x18c5ae as f32);
        assert_eq!(imgset0.ImageBuffer(4)[0 * img_stride + 3], 0x1cc5ae as f32);

        assert_eq!(imgset0.ImageBuffer(5)[0 * img_stride + 0], 0x10c5af as f32);
        assert_eq!(imgset0.ImageBuffer(5)[0 * img_stride + 1], 0x14c5af as f32);
        assert_eq!(imgset0.ImageBuffer(5)[0 * img_stride + 2], 0x18c5af as f32);
        assert_eq!(imgset0.ImageBuffer(5)[0 * img_stride + 3], 0x1cc5af as f32);

        assert_eq!(imgset0.ImageBuffer(6)[0 * img_stride + 0], 0x10bec6 as f32);
        assert_eq!(imgset0.ImageBuffer(6)[0 * img_stride + 1], 0x14bec6 as f32);
        assert_eq!(imgset0.ImageBuffer(6)[0 * img_stride + 2], 0x18bec6 as f32);
        assert_eq!(imgset0.ImageBuffer(6)[0 * img_stride + 3], 0x1cbec6 as f32);

        assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 0], 0x10bec7 as f32);
        assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 1], 0x14bec7 as f32);
        assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 2], 0x18bec7 as f32);
        assert_eq!(imgset0.ImageBuffer(7)[0 * img_stride + 3], 0x1cbec7 as f32);

        let imgset5 = baseline_imgsets.get(&5).unwrap();

        assert_eq!(imgset5.ImageBuffer(0)[0 * img_stride + 0], 0x10f1ce as f32);
        assert_eq!(imgset5.ImageBuffer(0)[0 * img_stride + 1], 0x14f1ce as f32);
        assert_eq!(imgset5.ImageBuffer(0)[0 * img_stride + 2], 0x18f1ce as f32);
        assert_eq!(imgset5.ImageBuffer(0)[0 * img_stride + 3], 0x1cf1ce as f32);

        assert_eq!(imgset5.ImageBuffer(1)[0 * img_stride + 0], -0x10f1cf as f32);
        assert_eq!(imgset5.ImageBuffer(1)[0 * img_stride + 1], -0x14f1cf as f32);
        assert_eq!(imgset5.ImageBuffer(1)[0 * img_stride + 2], -0x18f1cf as f32);
        assert_eq!(imgset5.ImageBuffer(1)[0 * img_stride + 3], -0x1cf1cf as f32);

        assert_eq!(imgset5.ImageBuffer(2)[0 * img_stride + 0], 0x10ea26 as f32);
        assert_eq!(imgset5.ImageBuffer(2)[0 * img_stride + 1], 0x14ea26 as f32);
        assert_eq!(imgset5.ImageBuffer(2)[0 * img_stride + 2], 0x18ea26 as f32);
        assert_eq!(imgset5.ImageBuffer(2)[0 * img_stride + 3], 0x1cea26 as f32);

        assert_eq!(imgset5.ImageBuffer(3)[0 * img_stride + 0], -0x10ea27 as f32);
        assert_eq!(imgset5.ImageBuffer(3)[0 * img_stride + 1], -0x14ea27 as f32);
        assert_eq!(imgset5.ImageBuffer(3)[0 * img_stride + 2], -0x18ea27 as f32);
        assert_eq!(imgset5.ImageBuffer(3)[0 * img_stride + 3], -0x1cea27 as f32);

        assert_eq!(imgset5.ImageBuffer(4)[0 * img_stride + 0], 0x10f1be as f32);
        assert_eq!(imgset5.ImageBuffer(4)[0 * img_stride + 1], 0x14f1be as f32);
        assert_eq!(imgset5.ImageBuffer(4)[0 * img_stride + 2], 0x18f1be as f32);
        assert_eq!(imgset5.ImageBuffer(4)[0 * img_stride + 3], 0x1cf1be as f32);

        assert_eq!(imgset5.ImageBuffer(5)[0 * img_stride + 0], -0x10f1bf as f32);
        assert_eq!(imgset5.ImageBuffer(5)[0 * img_stride + 1], -0x14f1bf as f32);
        assert_eq!(imgset5.ImageBuffer(5)[0 * img_stride + 2], -0x18f1bf as f32);
        assert_eq!(imgset5.ImageBuffer(5)[0 * img_stride + 3], -0x1cf1bf as f32);

        assert_eq!(imgset5.ImageBuffer(6)[0 * img_stride + 0], 0x10ea16 as f32);
        assert_eq!(imgset5.ImageBuffer(6)[0 * img_stride + 1], 0x14ea16 as f32);
        assert_eq!(imgset5.ImageBuffer(6)[0 * img_stride + 2], 0x18ea16 as f32);
        assert_eq!(imgset5.ImageBuffer(6)[0 * img_stride + 3], 0x1cea16 as f32);

        assert_eq!(imgset5.ImageBuffer(7)[0 * img_stride + 0], -0x10ea17 as f32);
        assert_eq!(imgset5.ImageBuffer(7)[0 * img_stride + 1], -0x14ea17 as f32);
        assert_eq!(imgset5.ImageBuffer(7)[0 * img_stride + 2], -0x18ea17 as f32);
        assert_eq!(imgset5.ImageBuffer(7)[0 * img_stride + 3], -0x1cea17 as f32);
    }

    #[test]
    fn test_flag_imgsets_minimal() {
        let mut baseline_imgsets: BTreeMap<usize, UniquePtr<CxxImageSet>> = BTreeMap::new();
        let width = 64;
        let height = 64;
        let img_stride = (((width - 1) / 8) + 1) * 8;

        let noise_x = 32;
        let noise_y = 32;
        let noise_z = 1;
        let noise_val = 0xffffff as f32;
        let baseline_flagmasks = unsafe {
            let aoflagger = cxx_aoflagger_new();
            let imgset0 = aoflagger.MakeImageSet(width, height, 8, 0 as f32, width);
            baseline_imgsets.insert(0, imgset0);
            let mut imgset1 = aoflagger.MakeImageSet(width, height, 8, 0 as f32, width);
            imgset1.pin_mut().ImageBufferMut(noise_z)[noise_y * img_stride + noise_x] = noise_val;
            baseline_imgsets.insert(1, imgset1);

            let strategy_file_minimal = aoflagger.FindStrategyFileGeneric(&String::from("minimal"));
            flag_imgsets(&aoflagger, &strategy_file_minimal, baseline_imgsets)
        };

        let flagmask0 = baseline_flagmasks.get(&0).unwrap();
        let flagmask1 = baseline_flagmasks.get(&1).unwrap();
        let flag_stride = flagmask0.HorizontalStride();
        assert!(!flagmask0.Buffer()[0]);
        assert!(!flagmask0.Buffer()[noise_y * flag_stride + noise_x]);
        assert!(!flagmask1.Buffer()[0]);
        assert!(flagmask1.Buffer()[noise_y * flag_stride + noise_x]);
    }

    #[test]
    fn test_write_flags_mwax_minimal() {
        let context = get_mwax_context();
        let mut baseline_flagmasks: BTreeMap<usize, UniquePtr<CxxFlagMask>> = BTreeMap::new();

        let height =
            context.num_coarse_chans * context.metafits_context.num_corr_fine_chans_per_coarse;
        let width = context.num_timesteps;

        let flag_timestep = 1;
        let flag_channel = 1;
        let flag_baseline = 1;

        unsafe {
            let aoflagger = cxx_aoflagger_new();
            for (baseline_idx, _) in context.metafits_context.baselines.iter().enumerate() {
                baseline_flagmasks
                    .insert(baseline_idx, aoflagger.MakeFlagMask(width, height, false));
            }
        }
        let flagmask = baseline_flagmasks.get_mut(&flag_baseline).unwrap();
        let flag_stride = flagmask.HorizontalStride();
        flagmask.pin_mut().BufferMut()[flag_channel * flag_stride + flag_timestep] = true;

        let tmp_dir = tempdir().unwrap();

        let gpubox_ids: Vec<usize> = context
            .coarse_chans
            .iter()
            .map(|chan| chan.gpubox_number)
            .collect();

        let filename_template = tmp_dir.path().join("Flagfile%%%.mwaf");
        let selected_gpuboxes = gpubox_ids[..1].to_vec();

        write_flags(
            &context,
            baseline_flagmasks,
            filename_template.to_str().unwrap(),
            &selected_gpuboxes,
        );

        let flag_files = glob(&tmp_dir.path().join("Flagfile*.mwaf").to_str().unwrap()).unwrap();

        assert_eq!(flag_files.count(), selected_gpuboxes.len());

        let mut flag_file_set = FlagFileSet::open(
            &context,
            filename_template.to_str().unwrap(),
            &selected_gpuboxes,
        )
        .unwrap();
        let chan_header_flags_raw = flag_file_set.read_chan_header_flags_raw().unwrap();
        let (chan1_header, chan1_flags_raw) =
            chan_header_flags_raw.get(&selected_gpuboxes[0]).unwrap();
        assert_eq!(chan1_header.gpubox_id, gpubox_ids[0]);
        let num_fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;

        let num_baselines = chan1_header.num_ants * (chan1_header.num_ants + 1) / 2;
        assert_eq!(chan1_header.num_timesteps, context.num_timesteps);
        assert_eq!(num_baselines, context.metafits_context.num_baselines);
        assert_eq!(chan1_header.num_channels, num_fine_chans_per_coarse);
        assert_eq!(
            chan1_flags_raw.len(),
            chan1_header.num_timesteps * num_baselines * chan1_header.num_channels
        );
        dbg!(&chan1_flags_raw);

        let tests = [
            (0, 0, 0, i8::from(false)),
            (0, 0, 1, i8::from(false)),
            (0, 1, 0, i8::from(false)),
            (0, 1, 1, i8::from(false)),
            (0, 2, 0, i8::from(false)),
            (0, 2, 1, i8::from(false)),
            (1, 0, 0, i8::from(false)),
            (1, 0, 1, i8::from(false)),
            (1, 1, 0, i8::from(false)),
            (1, 1, 1, i8::from(true)),
            (1, 2, 0, i8::from(false)),
            (1, 2, 1, i8::from(false)),
        ];
        for (timestep_idx, baseline_idx, fine_chan_idx, expected_flag) in tests.iter() {
            let row_idx = timestep_idx * num_baselines + baseline_idx;
            let offset = row_idx * num_fine_chans_per_coarse + fine_chan_idx;
            assert_eq!(
                &chan1_flags_raw[offset], expected_flag,
                "with timestep {}, baseline {}, fine_chan {}, expected {} at row_idx {}, offset {}",
                timestep_idx, baseline_idx, fine_chan_idx, expected_flag, row_idx, offset
            );
        }
    }
}
