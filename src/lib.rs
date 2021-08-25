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
//! use birli::{
//!     context_to_jones_array, flag_jones_array_existing, write_flags,
//!     cxx_aoflagger_new, get_flaggable_timesteps, init_flag_array,
//!     get_antenna_flags, mwalib,
//! };
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
//! // Determine which timesteps and coarse channels we want to use
//! let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
//! let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
//!
//! let img_timestep_range =
//!     *img_timestep_idxs.first().unwrap()..(*img_timestep_idxs.last().unwrap() + 1);
//! let img_coarse_chan_range =
//!     *img_coarse_chan_idxs.first().unwrap()..(*img_coarse_chan_idxs.last().unwrap() + 1);
//!
//! // Prepare our flagmasks with known bad antennae
//! let flag_array = init_flag_array(
//!     &context,
//!     &img_timestep_range,
//!     &img_coarse_chan_range,
//!     Some(get_antenna_flags(&context)),
//! );
//!
//! // load visibilities into our array of jones matrices
//! let (mut jones_array, flag_array) = context_to_jones_array(
//!     &context,
//!     &img_timestep_range,
//!     &img_coarse_chan_range,
//!     Some(flag_array),
//! );
//!
//! // use the default strategy file location for MWA
//! let strategy_filename = &aoflagger.FindStrategyFileMWA();
//!
//! // run the strategy on the imagesets, and get the resulting flagmasks for each baseline
//! let flag_array = flag_jones_array_existing(
//!     &aoflagger,
//!     &strategy_filename,
//!     &jones_array,
//!     Some(flag_array),
//!     true,
//! );
//!
//! // write the flags to disk as .mwaf
//! write_flags(&context, &flag_array, flag_template.to_str().unwrap(), &img_coarse_chan_range).unwrap();
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
use error::BirliError;
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};
use itertools::izip;
use ndarray::{parallel::prelude::*, Array3, Axis};
use ndarray::{ArrayBase, Dim, ViewRepr};
use std::{ops::Range, os::raw::c_short};

use mwalib::CorrelatorContext;

pub mod error;

// pub mod math;

pub mod io;
pub use io::{mwaf::FlagFileSet, uvfits::UvfitsWriter};

pub mod corrections;
pub use corrections::{correct_cable_lengths, correct_geometry};

pub mod flags;
pub use flags::{
    flag_jones_array_existing, flagmask_or, get_antenna_flags, get_flaggable_timesteps,
    init_flag_array, write_flags,
};

// pub mod util;

use log::{trace, warn};

use crossbeam_channel::unbounded;
use crossbeam_utils::thread;
// use rayon::prelude::*;

pub use mwa_rust_core;
pub use mwa_rust_core::{mwalib, Complex, Jones};
pub use mwalib::{fitsio, fitsio_sys};

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

/// Initialize a vector of length [`num_baselines`] containing [`CxxImageSet`] of dimensions
/// [`width`] and [`height`],
///
/// # Examples
///
/// ```rust
/// use birli::{init_baseline_imgsets, cxx_aoflagger_new, mwalib};
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
/// let context = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
///
/// let num_baselines = context.metafits_context.num_baselines;
/// let width = context.num_common_timesteps;
/// let height = context.num_common_coarse_chans * context.metafits_context.num_corr_fine_chans_per_coarse;
/// let baseline_imgsets = unsafe {
///     let aoflagger = cxx_aoflagger_new();
///     init_baseline_imgsets(&aoflagger, num_baselines, width, height)
/// };
/// ```
pub fn init_baseline_imgsets(
    aoflagger: &CxxAOFlagger,
    num_baselines: usize,
    width: usize,
    height: usize,
) -> Vec<UniquePtr<CxxImageSet>> {
    trace!("start init_baseline_imgsets");

    // Create a progress bar to show the status of allocating
    let allocation_progress = ProgressBar::new(num_baselines as u64);
    allocation_progress.set_style(
        ProgressStyle::default_bar()
            .template(
                "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
            )
            .progress_chars("=> "),
    );
    allocation_progress.set_message("allocating imgs");

    // Allocate vector of [`CxxImageSet`]s for each baseline
    let baseline_imgsets: Vec<UniquePtr<CxxImageSet>> = (0..num_baselines)
        .into_par_iter()
        .map(|_| {
            let imgset = unsafe { aoflagger.MakeImageSet(width, height, 8, 0 as f32, width) };
            allocation_progress.inc(1);
            imgset
        })
        .collect();

    allocation_progress.finish();

    trace!("end init_baseline_imgsets");
    baseline_imgsets
}

/// Given a buffer from mwalib::CorrelatorContext.read_by_baseline, which is in the order
/// [baseline][chan][pol][complex]
/// Write these visibilities to a jones matrix array view for a given timestep, coarse channel
/// [chan][baseline]
macro_rules! _write_hdu_buffer_to_jones_view {
    ($hdu_buffer:expr, $jones_hdu_view:expr, $floats_per_baseline:expr, $floats_per_chan:expr) => {
        // jones_hdu_view is [chan][baseline] for the coarse channel and timestep (HDU)
        // hdu_buffer is [baseline][chan][pol][complex] for the coarse channel and timestep (HDU)
        for (mut jones_baseline_view, hdu_baseline_chunk) in izip!(
            $jones_hdu_view.axis_iter_mut(Axis(1)),
            $hdu_buffer.chunks($floats_per_baseline)
        ) {
            // jones_baseline_view is [chan] for the given baseline in the hdu
            // hdu_baseline_chunk is [chan][pol][complex]
            for (mut jones_chan_view, hdu_chan_chunk) in izip!(
                jones_baseline_view.outer_iter_mut(),
                hdu_baseline_chunk.chunks($floats_per_chan)
            ) {
                // Sanity check that our ndarray view is a single cell.
                assert_eq!(jones_chan_view.ndim(), 0);
                let jones = Jones::from([
                    Complex::new(hdu_chan_chunk[0], hdu_chan_chunk[1]),
                    Complex::new(hdu_chan_chunk[2], hdu_chan_chunk[3]),
                    Complex::new(hdu_chan_chunk[4], hdu_chan_chunk[5]),
                    Complex::new(hdu_chan_chunk[6], hdu_chan_chunk[7]),
                ]);
                jones_chan_view.fill(jones);
            }
        }
    };
}

/// generate a 3 dimensional array of Jones matrices from an observation's 
/// [`mwalib::CorrelatorContext`], for all baselines, over a given range of 
/// mwalib timestep and coarse channel indices.
///
/// The dimensions of the array are:
///  - timestep
///  - channel
///  - baseline
///
/// An equally sized flag array is also returned with flags indicating when reading via mwalib
/// causes a GPUBoxError.
///
/// # Details:
///
/// Observations can sometimes be too large to fit in memory. This method will only load
/// visibilities from the provided timesteps and coarse channels, in order to enable the visibility to
/// be read in user-defined "chunks" of time or frequency.
///
/// The timesteps are indices in the [`mwalib::CorrelatorContext`]'s timestep array, which should be a contiguous
/// superset of times from all provided coarse gpubox files. A similar concept applies to coarse channels.
/// Instead of reading visibilities for all known timesteps / coarse channels, it is recommended to use
/// `common_coarse_chan_indices` and `common_timestep_indices`, as these ignore timesteps and coarse channels
/// which are missing contiguous data. `common_good_timestep_indices` is also a good choice to avoid quack time.
///
/// For more details, see the [documentation](https://docs.rs/mwalib/latest/mwalib/struct.CorrelatorContext.html).
///
/// Note: it doesn't make sense to ask aoflagger to flag non-contiguous timesteps
/// or coarse channels, and so this interface only allows to ranges to be used. 
/// For flagging an obeservation with "picket fence" coarse channels or timesteps, 
/// contiguous ranges should be flagged separately.
///
/// [`mwalib::CorrelatorContext`]: https://docs.rs/mwalib/latest/mwalib/struct.CorrelatorContext.html
///
/// # Examples
///
/// ```rust
/// use birli::{context_to_jones_array, cxx_aoflagger_new, mwalib};
/// use mwalib::CorrelatorContext;
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
/// let aoflagger = unsafe { cxx_aoflagger_new() };
///
/// // Determine which timesteps and coarse channels we want to use
/// let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
/// let img_timestep_idxs = &context.common_timestep_indices;
/// let good_timestep_idxs = &context.common_good_timestep_indices;
///
/// let img_timestep_range =
///     *img_timestep_idxs.first().unwrap()..(*img_timestep_idxs.last().unwrap() + 1);
/// let good_timestep_range =
///     *good_timestep_idxs.first().unwrap()..(*good_timestep_idxs.last().unwrap() + 1);
/// let img_coarse_chan_range =
///     *img_coarse_chan_idxs.first().unwrap()..(*img_coarse_chan_idxs.last().unwrap() + 1);
///
///
/// // read visibilities out of the gpubox files
/// let (jones_array, _) = context_to_jones_array(
///     &context,
///     &img_timestep_range,
///     &img_coarse_chan_range,
///     None
/// );
///
/// let dims_common = jones_array.dim();
///
/// // read visibilities out of the gpubox files
/// let (jones_array, _) = context_to_jones_array(
///     &context,
///     &good_timestep_range,
///     &img_coarse_chan_range,
///     None
/// );
///
/// let dims_good = jones_array.dim();
///
/// assert_ne!(dims_common, dims_good);
/// ```
///
/// # Assumptions
/// - img_coarse_chan_idxs and img_timestep_idxs are contiguous
pub fn context_to_jones_array(
    context: &CorrelatorContext,
    mwalib_timestep_range: &Range<usize>,
    mwalib_coarse_chan_range: &Range<usize>,
    // TODO: allow subset of baselines
    // mwalib_baseline_idxs: &[usize],
    flag_array: Option<Array3<bool>>,
) -> (Array3<Jones<f32>>, Array3<bool>) {
    trace!("start context_to_jones_array");

    // allocate our result

    let num_timesteps = mwalib_timestep_range.len();
    let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
    let num_coarse_chans = mwalib_coarse_chan_range.len();
    let num_chans = num_coarse_chans * fine_chans_per_coarse;
    let num_baselines = context.metafits_context.num_baselines;

    let shape = (num_timesteps, num_chans, num_baselines);

    let mut jones_array: Array3<Jones<f32>> = Array3::from_elem(shape, Jones::default());

    let mut flag_array: Array3<bool> = if let Some(flag_array_) = flag_array {
        assert_eq!(
            flag_array_.dim(),
            shape,
            "jones array and flag array should be the same dimensions"
        );
        flag_array_
    } else {
        Array3::from_elem(shape, false)
    };

    // since we are using read_by_by basline into buffer, the visibilities are read in order:
    // baseline,frequency,pol,r,i

    let floats_per_chan = context.metafits_context.num_visibility_pols * 2;
    let floats_per_baseline = floats_per_chan * fine_chans_per_coarse;
    let floats_per_hdu = floats_per_baseline * num_baselines;

    // A queue of errors
    let (tx_error, rx_error) = unbounded();

    // a progress bar containing the progress bars associated with loading the observation's HDUs
    let multi_progress = MultiProgress::with_draw_target(ProgressDrawTarget::stderr());
    // a vector of progress bars for the visibility reading progress of each channel.
    let read_progress: Vec<ProgressBar> = mwalib_coarse_chan_range
        .to_owned()
        .map(|mwalib_coarse_chan_idx| {
            let channel_progress = multi_progress.add(
                ProgressBar::new(num_timesteps as _)
                    .with_style(
                        ProgressStyle::default_bar()
                            .template("{msg:16}: [{wide_bar:.blue}] {pos:4}/{len:4}")
                            .progress_chars("=> "),
                    )
                    .with_position(0)
                    .with_message(format!("coarse_chan {:03}", mwalib_coarse_chan_idx)),
            );
            channel_progress.set_position(0);
            channel_progress
        })
        .collect();

    // The total reading progress.
    let total_progress = multi_progress.add(
        ProgressBar::new((num_timesteps * num_coarse_chans) as _)
            .with_style(
                ProgressStyle::default_bar()
                    .template(
                        "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
                    )
                    .progress_chars("=> "),
            )
            .with_position(0)
            .with_message("loading hdus"),
    );

    thread::scope(|scope| {
        // Spawn a thread to draw the progress bars.
        scope.spawn(|_| {
            multi_progress.join().unwrap();
        });

        total_progress.set_position(0);

        scope.spawn(|_| {
            for (mwalib_timestep_idx, mwalib_coarse_chan_idx, err) in rx_error {
                warn!(
                    "could not read hdu ts={}, cc={} {:?}",
                    mwalib_timestep_idx, mwalib_coarse_chan_idx, err
                );

                assert!(mwalib_timestep_range.contains(&mwalib_timestep_idx));
                let img_timestep_idx = mwalib_timestep_idx - mwalib_timestep_range.start;
                assert!(mwalib_coarse_chan_range.contains(&mwalib_coarse_chan_idx));
                let img_coarse_chan_idx = mwalib_coarse_chan_idx - mwalib_coarse_chan_range.start;

                // TODO: there's probably a better way of doing this.
                for mut coarse_chan_flags_view in flag_array
                    .axis_chunks_iter_mut(Axis(1), fine_chans_per_coarse)
                    .skip(img_coarse_chan_idx)
                    .take(1)
                {
                    for mut hdu_flags_view in coarse_chan_flags_view
                        .outer_iter_mut()
                        .skip(img_timestep_idx)
                        .take(1)
                    {
                        hdu_flags_view.fill(true);
                    }
                }
            }
        });

        // Load HDUs from each coarse channel in parallel.
        jones_array
            .axis_chunks_iter_mut(Axis(1), fine_chans_per_coarse)
            .into_par_iter()
            .zip(mwalib_coarse_chan_range.clone())
            .zip(read_progress)
            .for_each(
                |((mut jones_coarse_chan_view, mwalib_coarse_chan_idx), progress)| {
                    progress.set_position(0);

                    let mut hdu_buffer: Vec<f32> = vec![0.0; floats_per_hdu];

                    // jones_coarse_chan_view is [timestep][chan][baseline] for all chans in the coarse channel
                    for (mwalib_timestep_idx, mut jones_hdu_view) in izip!(
                        mwalib_timestep_range.clone(),
                        jones_coarse_chan_view.outer_iter_mut()
                    ) {
                        match context.read_by_baseline_into_buffer(
                            mwalib_timestep_idx,
                            mwalib_coarse_chan_idx,
                            hdu_buffer.as_mut_slice(),
                        ) {
                            Ok(()) => {
                                _write_hdu_buffer_to_jones_view!(
                                    hdu_buffer,
                                    jones_hdu_view,
                                    floats_per_baseline,
                                    floats_per_chan
                                );
                            }
                            Err(err) => {
                                tx_error
                                    .send((mwalib_timestep_idx, mwalib_coarse_chan_idx, err))
                                    .unwrap();
                            }
                        }

                        progress.inc(1);
                        total_progress.inc(1);
                    }
                    progress.finish();
                },
            );

        drop(tx_error);

        // We're done!
        total_progress.finish();
    })
    .unwrap();

    trace!("end context_to_jones_array");

    (jones_array, flag_array)
}

/// Create an aoflagger [`CxxImageSet`] for a particular baseline from the given jones array
///
/// # Assumptions
///
/// - `baseline_jones_view` is [timestep][channel] for one baseline
/// - imageset is timesteps wide, and channels high
/// - jones matrics are always XX, YY, XY, YX
///
/// # Errors
///
/// Will throw [`BirliError`] if there was an error reading from context.
pub fn jones_baseline_view_to_imageset(
    aoflagger: &CxxAOFlagger,
    // jones_array: &Array3<Jones<f32>>,
    baseline_jones_view: &ArrayBase<ViewRepr<&Jones<f32>>, Dim<[usize; 2]>>,
) -> Result<UniquePtr<CxxImageSet>, BirliError> {
    let array_dims = baseline_jones_view.dim();
    let img_count = 8;
    let imgset = unsafe {
        aoflagger.MakeImageSet(
            array_dims.0,
            array_dims.1,
            img_count,
            0 as f32,
            array_dims.0,
        )
    };
    let img_stride = imgset.HorizontalStride();
    let mut img_bufs: Vec<&mut [f32]> = (0..img_count)
        .map(|img_idx| unsafe { imgset.ImageBufferMutUnsafe(img_idx) })
        .collect();

    // TODO: benchmark if iterate over pol first
    for (img_timestep_idx, timestep_jones_view) in baseline_jones_view.outer_iter().enumerate() {
        for (img_chan_idx, singular_jones_view) in timestep_jones_view.outer_iter().enumerate() {
            let jones = singular_jones_view.get(()).unwrap();
            for (img_idx, img_buf) in img_bufs.iter_mut().enumerate() {
                let pol_idx = img_idx / 2;
                img_buf[img_chan_idx * img_stride + img_timestep_idx] = if img_idx % 2 == 0 {
                    jones[pol_idx].re
                } else {
                    jones[pol_idx].im
                };
            }
        }
    }

    Ok(imgset)
}

/// Create an aoflagger [`CxxFlagMask`] for a from the given flag array view
///
/// # Assumptions
///
/// - flag array view is [timestep][channel] for one baseline
/// - flagmask is timesteps wide, and channels high
///
/// # Errors
///
/// Will throw [`BirliError`] if there was an error reading from context.
pub fn flag_baseline_view_to_flagmask(
    aoflagger: &CxxAOFlagger,
    baseline_flag_view: &ArrayBase<ViewRepr<&bool>, Dim<[usize; 2]>>,
) -> Result<UniquePtr<CxxFlagMask>, BirliError> {
    let array_dims = baseline_flag_view.dim();
    let mut flag_mask = unsafe { aoflagger.MakeFlagMask(array_dims.0, array_dims.1, false) };
    let stride = flag_mask.HorizontalStride();
    let flag_buf = flag_mask.pin_mut().BufferMut();

    // TODO: assign by slice
    for (img_timestep_idx, timestep_flag_view) in baseline_flag_view.outer_iter().enumerate() {
        for (img_chan_idx, singular_flag_view) in timestep_flag_view.outer_iter().enumerate() {
            flag_buf[img_chan_idx * stride + img_timestep_idx] =
                *singular_flag_view.get(()).unwrap();
        }
    }
    Ok(flag_mask)
}

#[cfg(test)]
mod tests {
    // TODO: Why does clippy think CxxImageSet.ImageBuffer() is &[f64]?
    #![allow(clippy::float_cmp)]

    // use core::slice::SlicePattern;

    use super::{context_to_jones_array, get_flaggable_timesteps};
    use crate::Jones;
    use approx::assert_abs_diff_eq;
    use mwa_rust_core::{mwalib, Complex};
    use mwalib::CorrelatorContext;
    // use num_complex::Complex;

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

    /// Get a dummy MWA Ord context with multiple holes in the data
    ///
    /// The gpubox (batch, hdu) tuples look like this:
    /// - ts is according to [`mwalib::correlatorContext`]
    ///
    /// |                   | ts=0   | 1      | 2      | 3      | 4      |
    /// | ----------------- | ------ | ------ | ------ | ------ | ------ |
    /// | gpubox=00         | (0, 0) | (0, 1) | .      | (1, 0) | .      |
    /// | 01                | .      | (0, 0) | (0, 1) | (1, 0) | (1, 1) |
    fn get_mwa_ord_dodgy_context() -> CorrelatorContext {
        let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
        let gpufits_paths = vec![
            "tests/data/1196175296_mwa_ord/adjusted_-1/1196175296_20171201145440_gpubox01_00.fits",
            "tests/data/1196175296_mwa_ord/limited_1/1196175296_20171201145540_gpubox01_01.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits",
        ];
        CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap()
    }

    #[test]
    /// We expect coarse channel 0 ( fine channels 0,1 ) to be the same as in get_mwa_ord_context,
    /// but coarse channel 0 (fine channels 2, 3 ) should be shifted.
    fn test_context_to_jones_array_mwax_flags_missing_hdus() {
        let context = get_mwa_ord_dodgy_context();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(img_timestep_idxs.len(), 4);
        let img_timestep_range =
            *img_timestep_idxs.first().unwrap()..(*img_timestep_idxs.last().unwrap() + 1);
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let img_coarse_chan_range =
            *img_coarse_chan_idxs.first().unwrap()..(*img_coarse_chan_idxs.last().unwrap() + 1);

        // let img_baseline_idxs: Vec<usize> = (0..context.metafits_context.num_baselines).collect();

        let (jones_array, flags_array) = context_to_jones_array(
            &context,
            &img_timestep_range,
            &img_coarse_chan_range,
            // img_baseline_idxs.as_slice(),
            None,
        );

        // ts 0, chan 0, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((0, 0, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x10c5be as f32, -0x10c5bf as f32),
                Complex::new(0x10c5ae as f32, 0x10c5af as f32),
                Complex::new(0x10c5ae as f32, -0x10c5af as f32),
                Complex::new(0x10bec6 as f32, -0x10bec7 as f32),
            ])
        );
        // ts 1, chan 0, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((1, 0, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x14c5be as f32, -0x14c5bf as f32),
                Complex::new(0x14c5ae as f32, 0x14c5af as f32),
                Complex::new(0x14c5ae as f32, -0x14c5af as f32),
                Complex::new(0x14bec6 as f32, -0x14bec7 as f32),
            ])
        );
        // ts 2, chan 0, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((2, 0, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x18c5be as f32, -0x18c5bf as f32),
                Complex::new(0x18c5ae as f32, 0x18c5af as f32),
                Complex::new(0x18c5ae as f32, -0x18c5af as f32),
                Complex::new(0x18bec6 as f32, -0x18bec7 as f32),
            ])
        );
        // ts 3, chan 0, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((3, 0, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x1cc5be as f32, -0x1cc5bf as f32),
                Complex::new(0x1cc5ae as f32, 0x1cc5af as f32),
                Complex::new(0x1cc5ae as f32, -0x1cc5af as f32),
                Complex::new(0x1cbec6 as f32, -0x1cbec7 as f32),
            ])
        );

        // Fine channel 2 is in the modified coarse channel, so should have drifted

        // ts 0, chan 2, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((0, 2, 0)).unwrap(),
            // Normally this would be:
            // &Jones::from([
            //     Complex::new(0x00c5be as f32, -0x00c5bf as f32),
            //     Complex::new(0x00c5ae as f32, 0x00c5af as f32),
            //     Complex::new(0x00c5ae as f32, -0x00c5af as f32),
            //     Complex::new(0x00bec6 as f32, -0x00bec7 as f32),
            // ])
            &Jones::from([
                Complex::new(0x04c5be as f32, -0x04c5bf as f32),
                Complex::new(0x04c5ae as f32, 0x04c5af as f32),
                Complex::new(0x04c5ae as f32, -0x04c5af as f32),
                Complex::new(0x04bec6 as f32, -0x04bec7 as f32),
            ])
        );
        assert!(!flags_array.get((0, 2, 0)).unwrap());
        // ts 1, chan 2, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((1, 2, 0)).unwrap(),
            &Jones::from([
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ])
        );
        assert!(flags_array.get((1, 2, 0)).unwrap());
        // ts 2, chan 2, baseline 0 - unchanged
        assert_abs_diff_eq!(
            jones_array.get((2, 2, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x08c5be as f32, -0x08c5bf as f32),
                Complex::new(0x08c5ae as f32, 0x08c5af as f32),
                Complex::new(0x08c5ae as f32, -0x08c5af as f32),
                Complex::new(0x08bec6 as f32, -0x08bec7 as f32),
            ])
        );
        assert!(!flags_array.get((2, 2, 0)).unwrap());
        // ts 3, chan 2, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((3, 2, 0)).unwrap(),
            &Jones::from([
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ])
        );
        assert!(flags_array.get((3, 2, 0)).unwrap());
    }

    #[test]
    fn test_context_to_jones_array_mwax() {
        let context = get_mwax_context();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(img_timestep_idxs.len(), 4);
        let img_timestep_range =
            *img_timestep_idxs.first().unwrap()..(*img_timestep_idxs.last().unwrap() + 1);
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let img_coarse_chan_range =
            *img_coarse_chan_idxs.first().unwrap()..(*img_coarse_chan_idxs.last().unwrap() + 1);

        // let img_baseline_idxs: Vec<usize> = (0..context.metafits_context.num_baselines).collect();

        let (jones_array, _) = context_to_jones_array(
            &context,
            &img_timestep_range,
            &img_coarse_chan_range,
            // img_baseline_idxs.as_slice(),
            None,
        );

        // ts 0, chan 0, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((0, 0, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x410000 as f32, 0x410001 as f32),
                Complex::new(0x410002 as f32, 0x410003 as f32),
                Complex::new(0x410004 as f32, 0x410005 as f32),
                Complex::new(0x410006 as f32, 0x410007 as f32),
            ])
        );

        // ts 1, chan 0, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((1, 0, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x410100 as f32, 0x410101 as f32),
                Complex::new(0x410102 as f32, 0x410103 as f32),
                Complex::new(0x410104 as f32, 0x410105 as f32),
                Complex::new(0x410106 as f32, 0x410107 as f32),
            ])
        );

        // ts 2, chan 0, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((2, 0, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x410200 as f32, 0x410201 as f32),
                Complex::new(0x410202 as f32, 0x410203 as f32),
                Complex::new(0x410204 as f32, 0x410205 as f32),
                Complex::new(0x410206 as f32, 0x410207 as f32),
            ])
        );

        // ts 0, chan 0, baseline 2
        assert_abs_diff_eq!(
            jones_array.get((0, 0, 2)).unwrap(),
            &Jones::from([
                Complex::new(0x410020 as f32, 0x410021 as f32),
                Complex::new(0x410022 as f32, 0x410023 as f32),
                Complex::new(0x410024 as f32, 0x410025 as f32),
                Complex::new(0x410026 as f32, 0x410027 as f32),
            ])
        );

        // ts 0, chan 1, baseline 2
        assert_abs_diff_eq!(
            jones_array.get((0, 1, 2)).unwrap(),
            &Jones::from([
                Complex::new(0x410028 as f32, 0x410029 as f32),
                Complex::new(0x41002a as f32, 0x41002b as f32),
                Complex::new(0x41002c as f32, 0x41002d as f32),
                Complex::new(0x41002e as f32, 0x41002f as f32),
            ])
        );
    }

    #[test]
    fn test_context_to_jones_array_mwa_ord() {
        let context = get_mwa_ord_context();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(img_timestep_idxs.len(), 4);
        let img_timestep_range =
            *img_timestep_idxs.first().unwrap()..(*img_timestep_idxs.last().unwrap() + 1);
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let img_coarse_chan_range =
            *img_coarse_chan_idxs.first().unwrap()..(*img_coarse_chan_idxs.last().unwrap() + 1);

        // let img_baseline_idxs: Vec<usize> = (0..context.metafits_context.num_baselines).collect();

        let (jones_array, _) = context_to_jones_array(
            &context,
            &img_timestep_range,
            &img_coarse_chan_range,
            // img_baseline_idxs.as_slice(),
            None,
        );

        // ts 0, chan 0, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((0, 0, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x10c5be as f32, -0x10c5bf as f32),
                Complex::new(0x10c5ae as f32, 0x10c5af as f32),
                Complex::new(0x10c5ae as f32, -0x10c5af as f32),
                Complex::new(0x10bec6 as f32, -0x10bec7 as f32),
            ])
        );
        // ts 1, chan 0, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((1, 0, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x14c5be as f32, -0x14c5bf as f32),
                Complex::new(0x14c5ae as f32, 0x14c5af as f32),
                Complex::new(0x14c5ae as f32, -0x14c5af as f32),
                Complex::new(0x14bec6 as f32, -0x14bec7 as f32),
            ])
        );
        // ts 2, chan 0, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((2, 0, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x18c5be as f32, -0x18c5bf as f32),
                Complex::new(0x18c5ae as f32, 0x18c5af as f32),
                Complex::new(0x18c5ae as f32, -0x18c5af as f32),
                Complex::new(0x18bec6 as f32, -0x18bec7 as f32),
            ])
        );
        // ts 3, chan 0, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((3, 0, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x1cc5be as f32, -0x1cc5bf as f32),
                Complex::new(0x1cc5ae as f32, 0x1cc5af as f32),
                Complex::new(0x1cc5ae as f32, -0x1cc5af as f32),
                Complex::new(0x1cbec6 as f32, -0x1cbec7 as f32),
            ])
        );

        // ts 0, chan 0, baseline 5
        assert_abs_diff_eq!(
            jones_array.get((0, 0, 5)).unwrap(),
            &Jones::from([
                Complex::new(0x10f1ce as f32, -0x10f1cf as f32),
                Complex::new(0x10ea26 as f32, -0x10ea27 as f32),
                Complex::new(0x10f1be as f32, -0x10f1bf as f32),
                Complex::new(0x10ea16 as f32, -0x10ea17 as f32),
            ])
        );
        // ts 1, chan 0, baseline 5
        assert_abs_diff_eq!(
            jones_array.get((1, 0, 5)).unwrap(),
            &Jones::from([
                Complex::new(0x14f1ce as f32, -0x14f1cf as f32),
                Complex::new(0x14ea26 as f32, -0x14ea27 as f32),
                Complex::new(0x14f1be as f32, -0x14f1bf as f32),
                Complex::new(0x14ea16 as f32, -0x14ea17 as f32),
            ])
        );
        // ts 2, chan 0, baseline 5
        assert_abs_diff_eq!(
            jones_array.get((2, 0, 5)).unwrap(),
            &Jones::from([
                Complex::new(0x18f1ce as f32, -0x18f1cf as f32),
                Complex::new(0x18ea26 as f32, -0x18ea27 as f32),
                Complex::new(0x18f1be as f32, -0x18f1bf as f32),
                Complex::new(0x18ea16 as f32, -0x18ea17 as f32),
            ])
        );
        // ts 3, chan 0, baseline 5
        assert_abs_diff_eq!(
            jones_array.get((3, 0, 5)).unwrap(),
            &Jones::from([
                Complex::new(0x1cf1ce as f32, -0x1cf1cf as f32),
                Complex::new(0x1cea26 as f32, -0x1cea27 as f32),
                Complex::new(0x1cf1be as f32, -0x1cf1bf as f32),
                Complex::new(0x1cea16 as f32, -0x1cea17 as f32),
            ])
        );

        // ts 0, chan 2, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((0, 2, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x00c5be as f32, -0x00c5bf as f32),
                Complex::new(0x00c5ae as f32, 0x00c5af as f32),
                Complex::new(0x00c5ae as f32, -0x00c5af as f32),
                Complex::new(0x00bec6 as f32, -0x00bec7 as f32),
            ])
        );
        // ts 1, chan 2, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((1, 2, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x04c5be as f32, -0x04c5bf as f32),
                Complex::new(0x04c5ae as f32, 0x04c5af as f32),
                Complex::new(0x04c5ae as f32, -0x04c5af as f32),
                Complex::new(0x04bec6 as f32, -0x04bec7 as f32),
            ])
        );
        // ts 2, chan 2, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((2, 2, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x08c5be as f32, -0x08c5bf as f32),
                Complex::new(0x08c5ae as f32, 0x08c5af as f32),
                Complex::new(0x08c5ae as f32, -0x08c5af as f32),
                Complex::new(0x08bec6 as f32, -0x08bec7 as f32),
            ])
        );
        // ts 3, chan 2, baseline 0
        assert_abs_diff_eq!(
            jones_array.get((3, 2, 0)).unwrap(),
            &Jones::from([
                Complex::new(0x0cc5be as f32, -0x0cc5bf as f32),
                Complex::new(0x0cc5ae as f32, 0x0cc5af as f32),
                Complex::new(0x0cc5ae as f32, -0x0cc5af as f32),
                Complex::new(0x0cbec6 as f32, -0x0cbec7 as f32),
            ])
        );
    }
}
