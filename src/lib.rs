#![warn(missing_docs)]
#![warn(missing_doc_code_examples)]
#![warn(clippy::missing_safety_doc)]
#![warn(clippy::missing_errors_doc)]

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
use std::os::raw::c_short;

pub mod flag_io;
use flag_io::FlagFileSet;

pub mod error;

use crossbeam_channel::{bounded, unbounded};
use crossbeam_utils::thread;

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

/// Read aan observation's visibilities into a vector containing a [`CxxImageSet`]s for each
/// baseline in the observation, given a [`CxxAOFlagger`] instance and that observation's
/// [`mwalib::CorrelatorContext`].
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
) -> Vec<UniquePtr<CxxImageSet>> {
    let coarse_chan_idxs: Vec<usize> = context
        .coarse_chans
        .iter()
        .enumerate()
        .map(|(idx, _)| idx)
        .collect();
    let timestep_idxs: Vec<usize> = context
        .timesteps
        .iter()
        .enumerate()
        .map(|(idx, _)| idx)
        .collect();

    let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
    let floats_per_finechan = context.metafits_context.num_visibility_pols * 2;
    let floats_per_baseline = fine_chans_per_coarse * floats_per_finechan;
    let height = context.num_coarse_chans * fine_chans_per_coarse;
    let width = context.num_timesteps;
    let img_stride = (((width - 1) / 8) + 1) * 8;

    let mut baseline_imgsets: Vec<UniquePtr<CxxImageSet>> = context
        .metafits_context
        .baselines
        .iter()
        .map(|_| unsafe { aoflagger.MakeImageSet(width, height, 8, 0 as f32, width) })
        .collect();

    let num_workers = coarse_chan_idxs.len();

    // A queue of coarse channel indices for the producers to work through
    let (tx_coarse_chan_idx, rx_coarse_chan_idx) = unbounded();
    // A queue of image buffers to process
    let (tx_img, rx_img) = bounded(num_workers);

    thread::scope(|scope| {
        // Queue up coarse channels to do
        for coarse_chan_idx in coarse_chan_idxs {
            tx_coarse_chan_idx.send(coarse_chan_idx).unwrap();
        }

        // Create a producer thread for worker
        for _ in 0..num_workers {
            let (_, rx_coarse_chan_idx_worker) =
                (tx_coarse_chan_idx.clone(), rx_coarse_chan_idx.clone());
            let (tx_img_worker, _) = (tx_img.clone(), rx_img.clone());
            let timestep_idxs_worker = timestep_idxs.to_owned();
            scope.spawn(move |_| {
                for coarse_chan_idx in rx_coarse_chan_idx_worker.iter() {
                    for &timestep_idx in timestep_idxs_worker.iter() {
                        let img_buf = context
                            .read_by_baseline(timestep_idx, coarse_chan_idx)
                            .unwrap();
                        let msg = (coarse_chan_idx, timestep_idx, img_buf);
                        tx_img_worker.send(msg).unwrap();
                    }
                }
            });
        }

        drop(tx_coarse_chan_idx);
        drop(tx_img);

        // consume the rx_img queue
        for (coarse_chan_idx, timestep_idx, img_buf) in rx_img.iter() {
            for (baseline_idx, baseline_chunk) in img_buf.chunks_exact(floats_per_baseline).enumerate() {
                let imgset = &mut baseline_imgsets[baseline_idx];

                for float_idx in 0..8 {
                    let imgset_buf = imgset.pin_mut().ImageBufferMut(float_idx);
                    for (fine_chan_idx, fine_chan_chunk) in
                        baseline_chunk.chunks_exact(floats_per_finechan).enumerate()
                    {
                        let img_x = timestep_idx;
                        let img_y = fine_chans_per_coarse * coarse_chan_idx + fine_chan_idx;
                        imgset_buf[img_y * img_stride + img_x] = fine_chan_chunk[float_idx];
                    }
                }
            }
        }
    })
    .unwrap();

    baseline_imgsets
}

/// Flag an observation's visibilities, given a [`CxxAOFlagger`] instance, a [`CxxStrategy`]
/// filename, and a vector of [`CxxImageSet`]s for each baseline in the observation returning a
/// vector of [`CxxFlagMask`]s.
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
    strategy_filename: &str,
    baseline_imgsets: Vec<UniquePtr<CxxImageSet>>,
) -> Vec<UniquePtr<CxxFlagMask>> {
    // TODO: figure out how to parallelize with Rayon, into_iter(). You'll probably need to convert between UniquePtr and Box

    baseline_imgsets
        .iter()
        .map(|imgset| {
            aoflagger
                .LoadStrategyFile(&strategy_filename.to_string())
                .Run(&imgset)
        })
        .collect()
}

/// Write flags to disk, given an observation's [`mwalib::CorrelatorContext`], a vector of
/// [`CxxFlagMask`]s for each baseline in the observation, a filename template and a vector of
/// gpubox IDs.
///
/// The filename template should contain two or 3 percentage (`%`) characters which will be replaced
/// by the gpubox id or channel number (depending on correlator type) provided in `gpubox_ids`. See
/// [`flag_io::FlagFileSet::new`] for more details.
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
    baseline_flagmasks: Vec<UniquePtr<CxxFlagMask>>,
    filename_template: &str,
    gpubox_ids: &[usize],
) {
    let mut flag_file_set = FlagFileSet::new(context, filename_template, &gpubox_ids).unwrap();
    flag_file_set
        .write_baseline_flagmasks(&context, baseline_flagmasks)
        .unwrap();
}

#[cfg(test)]
mod tests {
    // TODO: Why does clippy think CxxImageSet.ImageBuffer() is &[f64]?
    #![allow(clippy::float_cmp)]

    use super::{context_to_baseline_imgsets, flag_imgsets, write_flags, FlagFileSet};
    use crate::cxx_aoflagger::ffi::{cxx_aoflagger_new, CxxFlagMask, CxxImageSet};
    use cxx::UniquePtr;
    use glob::glob;
    use mwalib::CorrelatorContext;
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

    macro_rules! test_imgset_val {
        ($imgset:expr, $imgset_idx:expr, $img_stride:expr, $x:expr, $y:expr, $val:expr) => {
            assert_eq!(
                $imgset.ImageBuffer($imgset_idx)[$x * $img_stride + $y],
                $val
            )
        };
    }

    #[test]
    fn test_context_to_baseline_imgsets_mwax() {
        let context = get_mwax_context();
        let width = context.num_timesteps;
        let img_stride = (((width - 1) / 8) + 1) * 8;

        let aoflagger = unsafe { cxx_aoflagger_new() };
        let baseline_imgsets = context_to_baseline_imgsets(&aoflagger, &context);

        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 0, 0x410000 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 1, 0x410100 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 2, 0x410200 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 3, 0x410300 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 4, 0.0);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 5, 0.0);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 6, 0.0);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 7, 0.0);

        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 1, 0, 0x410008 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 1, 1, 0x410108 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 1, 2, 0x410208 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 1, 3, 0x410308 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 1, 4, 0.0);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 1, 5, 0.0);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 1, 6, 0.0);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 1, 7, 0.0);

        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 2, 0, 0x410400 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 2, 1, 0x410500 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 2, 2, 0x410600 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 2, 3, 0x410700 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 2, 4, 0.0);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 2, 5, 0.0);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 2, 6, 0.0);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 2, 7, 0.0);

        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 3, 0, 0x410408 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 3, 1, 0x410508 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 3, 2, 0x410608 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 3, 3, 0x410708 as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 3, 4, 0.0);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 3, 5, 0.0);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 3, 6, 0.0);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 3, 7, 0.0);

        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 0, 0x410001 as f32);
        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 1, 0x410101 as f32);
        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 2, 0x410201 as f32);
        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 3, 0x410301 as f32);
        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 4, 0.0);
        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 5, 0.0);
        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 6, 0.0);
        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 7, 0.0);

        /* ... */
        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 0, 0x410007 as f32);
        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 1, 0x410107 as f32);
        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 2, 0x410207 as f32);
        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 3, 0x410307 as f32);
        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 4, 0.0);
        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 5, 0.0);
        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 6, 0.0);
        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 7, 0.0);

        /* ... */

        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 0, 0, 0x410020 as f32);
        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 0, 1, 0x410120 as f32);
        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 0, 2, 0x410220 as f32);
        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 0, 3, 0x410320 as f32);
        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 0, 4, 0.0);
        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 0, 5, 0.0);
        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 0, 6, 0.0);
        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 0, 7, 0.0);

        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 1, 0, 0x410028 as f32);
        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 1, 1, 0x410128 as f32);
        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 1, 2, 0x410228 as f32);
        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 1, 3, 0x410328 as f32);
        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 1, 4, 0.0);
        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 1, 5, 0.0);
        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 1, 6, 0.0);
        test_imgset_val!(baseline_imgsets[2], 0, img_stride, 1, 7, 0.0);
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

        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 0, 0x10c5be as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 1, 0x14c5be as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 2, 0x18c5be as f32);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 3, 0x1cc5be as f32);

        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 0, 0x10c5bf as f32);
        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 1, 0x14c5bf as f32);
        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 2, 0x18c5bf as f32);
        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 3, 0x1cc5bf as f32);

        test_imgset_val!(baseline_imgsets[0], 2, img_stride, 0, 0, 0x10c5ae as f32);
        test_imgset_val!(baseline_imgsets[0], 2, img_stride, 0, 1, 0x14c5ae as f32);
        test_imgset_val!(baseline_imgsets[0], 2, img_stride, 0, 2, 0x18c5ae as f32);
        test_imgset_val!(baseline_imgsets[0], 2, img_stride, 0, 3, 0x1cc5ae as f32);

        test_imgset_val!(baseline_imgsets[0], 3, img_stride, 0, 0, -0x10c5af as f32);
        test_imgset_val!(baseline_imgsets[0], 3, img_stride, 0, 1, -0x14c5af as f32);
        test_imgset_val!(baseline_imgsets[0], 3, img_stride, 0, 2, -0x18c5af as f32);
        test_imgset_val!(baseline_imgsets[0], 3, img_stride, 0, 3, -0x1cc5af as f32);

        test_imgset_val!(baseline_imgsets[0], 4, img_stride, 0, 0, 0x10c5ae as f32);
        test_imgset_val!(baseline_imgsets[0], 4, img_stride, 0, 1, 0x14c5ae as f32);
        test_imgset_val!(baseline_imgsets[0], 4, img_stride, 0, 2, 0x18c5ae as f32);
        test_imgset_val!(baseline_imgsets[0], 4, img_stride, 0, 3, 0x1cc5ae as f32);

        test_imgset_val!(baseline_imgsets[0], 5, img_stride, 0, 0, 0x10c5af as f32);
        test_imgset_val!(baseline_imgsets[0], 5, img_stride, 0, 1, 0x14c5af as f32);
        test_imgset_val!(baseline_imgsets[0], 5, img_stride, 0, 2, 0x18c5af as f32);
        test_imgset_val!(baseline_imgsets[0], 5, img_stride, 0, 3, 0x1cc5af as f32);

        test_imgset_val!(baseline_imgsets[0], 6, img_stride, 0, 0, 0x10bec6 as f32);
        test_imgset_val!(baseline_imgsets[0], 6, img_stride, 0, 1, 0x14bec6 as f32);
        test_imgset_val!(baseline_imgsets[0], 6, img_stride, 0, 2, 0x18bec6 as f32);
        test_imgset_val!(baseline_imgsets[0], 6, img_stride, 0, 3, 0x1cbec6 as f32);

        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 0, 0x10bec7 as f32);
        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 1, 0x14bec7 as f32);
        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 2, 0x18bec7 as f32);
        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 3, 0x1cbec7 as f32);

        /* ... */

        test_imgset_val!(baseline_imgsets[5], 0, img_stride, 0, 0, 0x10f1ce as f32);
        test_imgset_val!(baseline_imgsets[5], 0, img_stride, 0, 1, 0x14f1ce as f32);
        test_imgset_val!(baseline_imgsets[5], 0, img_stride, 0, 2, 0x18f1ce as f32);
        test_imgset_val!(baseline_imgsets[5], 0, img_stride, 0, 3, 0x1cf1ce as f32);

        test_imgset_val!(baseline_imgsets[5], 1, img_stride, 0, 0, -0x10f1cf as f32);
        test_imgset_val!(baseline_imgsets[5], 1, img_stride, 0, 1, -0x14f1cf as f32);
        test_imgset_val!(baseline_imgsets[5], 1, img_stride, 0, 2, -0x18f1cf as f32);
        test_imgset_val!(baseline_imgsets[5], 1, img_stride, 0, 3, -0x1cf1cf as f32);

        test_imgset_val!(baseline_imgsets[5], 2, img_stride, 0, 0, 0x10ea26 as f32);
        test_imgset_val!(baseline_imgsets[5], 2, img_stride, 0, 1, 0x14ea26 as f32);
        test_imgset_val!(baseline_imgsets[5], 2, img_stride, 0, 2, 0x18ea26 as f32);
        test_imgset_val!(baseline_imgsets[5], 2, img_stride, 0, 3, 0x1cea26 as f32);

        test_imgset_val!(baseline_imgsets[5], 3, img_stride, 0, 0, -0x10ea27 as f32);
        test_imgset_val!(baseline_imgsets[5], 3, img_stride, 0, 1, -0x14ea27 as f32);
        test_imgset_val!(baseline_imgsets[5], 3, img_stride, 0, 2, -0x18ea27 as f32);
        test_imgset_val!(baseline_imgsets[5], 3, img_stride, 0, 3, -0x1cea27 as f32);

        test_imgset_val!(baseline_imgsets[5], 4, img_stride, 0, 0, 0x10f1be as f32);
        test_imgset_val!(baseline_imgsets[5], 4, img_stride, 0, 1, 0x14f1be as f32);
        test_imgset_val!(baseline_imgsets[5], 4, img_stride, 0, 2, 0x18f1be as f32);
        test_imgset_val!(baseline_imgsets[5], 4, img_stride, 0, 3, 0x1cf1be as f32);

        test_imgset_val!(baseline_imgsets[5], 5, img_stride, 0, 0, -0x10f1bf as f32);
        test_imgset_val!(baseline_imgsets[5], 5, img_stride, 0, 1, -0x14f1bf as f32);
        test_imgset_val!(baseline_imgsets[5], 5, img_stride, 0, 2, -0x18f1bf as f32);
        test_imgset_val!(baseline_imgsets[5], 5, img_stride, 0, 3, -0x1cf1bf as f32);

        test_imgset_val!(baseline_imgsets[5], 6, img_stride, 0, 0, 0x10ea16 as f32);
        test_imgset_val!(baseline_imgsets[5], 6, img_stride, 0, 1, 0x14ea16 as f32);
        test_imgset_val!(baseline_imgsets[5], 6, img_stride, 0, 2, 0x18ea16 as f32);
        test_imgset_val!(baseline_imgsets[5], 6, img_stride, 0, 3, 0x1cea16 as f32);

        test_imgset_val!(baseline_imgsets[5], 7, img_stride, 0, 0, -0x10ea17 as f32);
        test_imgset_val!(baseline_imgsets[5], 7, img_stride, 0, 1, -0x14ea17 as f32);
        test_imgset_val!(baseline_imgsets[5], 7, img_stride, 0, 2, -0x18ea17 as f32);
        test_imgset_val!(baseline_imgsets[5], 7, img_stride, 0, 3, -0x1cea17 as f32);
    }

    #[test]
    fn test_flag_imgsets_minimal() {
        let width = 64;
        let height = 64;
        let img_stride = (((width - 1) / 8) + 1) * 8;

        let noise_x = 32;
        let noise_y = 32;
        let noise_z = 1;
        let noise_val = 0xffffff as f32;

        let aoflagger = unsafe { cxx_aoflagger_new() };

        let mut baseline_imgsets: Vec<UniquePtr<CxxImageSet>> = (0..2)
            .into_iter()
            .map(|_| unsafe { aoflagger.MakeImageSet(width, height, 8, 0 as f32, width) })
            .collect();

        // modify the imageset to add some synthetic noise
        baseline_imgsets[1].pin_mut().ImageBufferMut(noise_z)[noise_y * img_stride + noise_x] =
            noise_val;

        let strategy_file_minimal = aoflagger.FindStrategyFileGeneric(&String::from("minimal"));

        let baseline_flagmasks = flag_imgsets(&aoflagger, &strategy_file_minimal, baseline_imgsets);

        let flag_stride = baseline_flagmasks[0].HorizontalStride();
        assert!(!baseline_flagmasks[0].Buffer()[0]);
        assert!(!baseline_flagmasks[0].Buffer()[noise_y * flag_stride + noise_x]);
        assert!(!baseline_flagmasks[1].Buffer()[0]);
        assert!(baseline_flagmasks[1].Buffer()[noise_y * flag_stride + noise_x]);
    }

    #[test]
    fn test_write_flags_mwax_minimal() {
        let context = get_mwax_context();

        let height =
            context.num_coarse_chans * context.metafits_context.num_corr_fine_chans_per_coarse;
        let width = context.num_timesteps;

        let flag_timestep = 1;
        let flag_channel = 1;
        let flag_baseline = 1;

        let aoflagger = unsafe { cxx_aoflagger_new() };
        let mut baseline_flagmasks: Vec<UniquePtr<CxxFlagMask>> = context
            .metafits_context
            .baselines
            .iter()
            .map(|_| unsafe { aoflagger.MakeFlagMask(width, height, false) })
            .collect();
        let flagmask = baseline_flagmasks.get_mut(flag_baseline).unwrap();
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
