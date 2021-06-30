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
//!     context_to_baseline_imgsets, flag_imgsets_existing, write_flags,
//!     cxx_aoflagger_new, get_flaggable_timesteps, init_baseline_flagmasks,
//!     get_antenna_flags
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
//! // Prepare our flagmasks with known bad antennae
//! let mut baseline_flagmasks = init_baseline_flagmasks(
//!     &aoflagger,
//!     &context,
//!     &img_coarse_chan_idxs,
//!     &img_timestep_idxs,
//!     Some(get_antenna_flags(&context)),
//! );
//!
//! // generate imagesets for each baseline in the format required by aoflagger
//! let baseline_imgsets = context_to_baseline_imgsets(
//!     &aoflagger,
//!     &context,
//!     &img_coarse_chan_idxs,
//!     &img_timestep_idxs,
//!     Some(&mut baseline_flagmasks),
//! );
//!
//! // use the default strategy file location for MWA
//! let strategy_filename = &aoflagger.FindStrategyFileMWA();
//!
//! // run the strategy on the imagesets, and get the resulting flagmasks for each baseline
//! flag_imgsets_existing(
//!     &aoflagger,
//!     &strategy_filename,
//!     &baseline_imgsets,
//!     &mut baseline_flagmasks,
//!     true
//! );
//!
//! // Get a list of all gpubox IDs
//! let gpubox_ids: Vec<usize> = context
//!     .common_coarse_chan_indices
//!     .iter()
//!     .map(|&chan| context.coarse_chans[chan].gpubox_number)
//!     .collect();
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
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};
use std::os::raw::c_short;

use mwalib::CorrelatorContext;
// use std::collections::BTreeMap;

use std::sync::Arc;

pub mod flags;
pub use flags::{
    flag_imgsets, flag_imgsets_existing, get_antenna_flags, get_flaggable_timesteps,
    init_baseline_flagmasks, write_flags,
};

pub mod io;

pub use io::mwaf::FlagFileSet;

pub mod corrections;
pub use corrections::correct_cable_lengths;

pub mod error;

pub mod util;

use log::{trace, warn};

use crossbeam_channel::{bounded, unbounded, Receiver, Sender};
use crossbeam_utils::thread;
use rayon::prelude::*;

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
/// use birli::{init_baseline_imgsets, cxx_aoflagger_new};
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

/// Read an observation's visibilities for the provided coarse channel and timestep indices into a
/// vector containing a [`CxxImageSet`]s for each baseline in the observation, given a
/// [`CxxAOFlagger`] instance and that observation's [`mwalib::CorrelatorContext`]
///
/// [`mwalib::CorrelatorContext`]: https://docs.rs/mwalib/latest/mwalib/struct.CorrelatorContext.html
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
/// or coarse channels. For flagging an obeservation with "picket fence"
/// coarse channels or timesteps, contiguous ranges should be flagged separately.
///
/// # TODO: ensure coarse channels and timesteps are contiguous, or use ranges.
///
/// # Examples
///
/// ```rust
/// use birli::{context_to_baseline_imgsets, cxx_aoflagger_new};
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
/// let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
///
/// let baseline_imgsets = context_to_baseline_imgsets(
///     &aoflagger,
///     &context,
///     &img_coarse_chan_idxs,
///     &context.common_timestep_indices.clone(),
///     None,
/// );
///
/// let width_common = baseline_imgsets[0].Width();
///
/// let baseline_imgsets_good = context_to_baseline_imgsets(
///     &aoflagger,
///     &context,
///     &img_coarse_chan_idxs,
///     &context.common_good_timestep_indices.clone(),
///     None,
/// );
///
/// let width_common_good = baseline_imgsets_good[0].Width();
///
/// assert_ne!(width_common, width_common_good);
/// ```
pub fn context_to_baseline_imgsets(
    aoflagger: &CxxAOFlagger,
    context: &CorrelatorContext,
    img_coarse_chan_idxs: &[usize],
    img_timestep_idxs: &[usize],
    baseline_flagmasks: Option<&mut Vec<UniquePtr<CxxFlagMask>>>,
) -> Vec<UniquePtr<CxxImageSet>> {
    trace!("start context_to_baseline_imgsets");

    // TODO: although collect::Vec<_> is more readable, it uses the stack less
    // efficiently than Vec
    // TODO: error handling instead of unwrap.

    let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
    let floats_per_finechan = context.metafits_context.num_visibility_pols * 2;
    let floats_per_baseline = fine_chans_per_coarse * floats_per_finechan;

    let num_img_coarse_chans = img_coarse_chan_idxs.len();
    let num_img_timesteps = img_timestep_idxs.len();

    let mut baseline_imgsets = init_baseline_imgsets(
        &aoflagger,
        context.metafits_context.num_baselines,
        num_img_timesteps,
        num_img_coarse_chans * fine_chans_per_coarse,
    );
    let baseline_imgsets_arc = Arc::new(&mut baseline_imgsets);

    // TODO: this sucks.
    let default_baseline_flagmasks: &mut Vec<UniquePtr<CxxFlagMask>> = &mut vec![];
    let baseline_flagmasks = baseline_flagmasks.unwrap_or(default_baseline_flagmasks);
    // let baseline_flagmasks_arc = Arc::new();

    let num_producers = num_img_coarse_chans;
    let num_pols_complex = context.metafits_context.num_visibility_pols * 2;

    // A queue of common coarse channel indices for the producers to work through
    let (tx_img_coarse_chan_idx, rx_img_coarse_chan_idx) = unbounded();
    // a queue of raw gpufits visibility image buffers for each complex polarization.
    let pol_img_queues: Vec<(Sender<_>, Receiver<_>)> = (0..num_pols_complex)
        .map(|_| bounded(num_producers))
        .collect();

    let (tx_error, rx_error) = unbounded();

    // a progress bar containing the progress bars associated with loading the
    // observation's HDUs
    let multi_progress = MultiProgress::new();
    // a vector of progress bars for the visibility reading progress of each
    // channel.
    let read_progress: Vec<ProgressBar> = (0..num_img_coarse_chans)
        .map(|img_coarse_chan_idx| {
            multi_progress.add(
                ProgressBar::new(num_img_timesteps as _)
                    .with_style(
                        ProgressStyle::default_bar()
                            .template("{msg:16}: [{wide_bar}] {pos:4}/{len:4}")
                            .progress_chars("=> "),
                    )
                    .with_position(0)
                    .with_message(format!(
                        "coarse chan {:3}",
                        img_coarse_chan_idxs[img_coarse_chan_idx]
                    )),
            )
        })
        .collect();
    // a vector of progress bars for the visibility writing progress of each
    // complex polarization,
    let write_progress: Vec<ProgressBar> = (0..num_pols_complex)
        .map(|pol_idx| {
            multi_progress.add(
                ProgressBar::new((num_img_timesteps * num_img_coarse_chans) as _)
                    .with_style(
                        ProgressStyle::default_bar()
                            .template("{msg:16}: [{wide_bar}] {pos:4}/{len:4}")
                            .progress_chars("=> "),
                    )
                    .with_position(0)
                    .with_message(format!("pol {:3}", pol_idx)),
            )
        })
        .collect();
    // A progress bar for the total progress of visibility loading, with time
    // elapsed, percentage and ETA.
    let total_progress = multi_progress.add(
        ProgressBar::new((num_img_timesteps * num_img_coarse_chans * num_pols_complex) as _)
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
    // These [`std::sync::Arc`]s provide concurrent access to the inner values.
    let read_progress_arc = Arc::new(read_progress);
    let write_progress_arc = Arc::new(write_progress);
    let total_progress_arc = Arc::new(total_progress);

    thread::scope(|scope| {
        // Spawn a thread to deal with the progress bar
        scope.spawn(|_| {
            multi_progress.join().unwrap();
        });

        // Queue up coarse channels to do
        (0..num_img_coarse_chans).for_each(|img_coarse_chan_idx| {
            tx_img_coarse_chan_idx.send(img_coarse_chan_idx).unwrap();
        });
        // This indicates the the producer threads that there are no more coarse
        // channels to process.
        drop(tx_img_coarse_chan_idx);

        // Create multiple producer threads
        for _ in 0..num_producers {
            let rx_img_coarse_chan_idx_worker = rx_img_coarse_chan_idx.clone();
            let pol_img_queues_worker = pol_img_queues.clone();
            let read_progress_worker = read_progress_arc.clone();
            let tx_error_worker = tx_error.clone();
            scope.spawn(move |_| {
                // Each producer thread consumes the HDUs from one coarse
                // channel at a time.
                for img_coarse_chan_idx in rx_img_coarse_chan_idx_worker.iter() {
                    let coarse_chan_idx = img_coarse_chan_idxs[img_coarse_chan_idx];
                    (0..num_img_timesteps).for_each(|img_timestep_idx| {
                        let timestep_idx = img_timestep_idxs[img_timestep_idx];
                        match context.read_by_baseline(timestep_idx, coarse_chan_idx) {
                            Ok(img_buf) => {
                                let img_buf_arc = Arc::new(img_buf);
                                // Producer sends the visibility buffer to a separate
                                // queue for each complex polarization.
                                pol_img_queues_worker.iter().for_each(|(tx_img, _)| {
                                    tx_img
                                        .send((
                                            img_coarse_chan_idx,
                                            img_timestep_idx,
                                            img_buf_arc.clone(),
                                        ))
                                        .unwrap();
                                });
                            }
                            Err(err) => {
                                tx_error_worker
                                    .send((img_coarse_chan_idx, img_timestep_idx, err))
                                    .unwrap();
                            }
                        }
                        read_progress_worker[img_coarse_chan_idx].inc(1);
                    });
                    read_progress_worker[img_coarse_chan_idx].finish_and_clear();
                }
                pol_img_queues_worker.into_iter().for_each(|(tx_img, _)| {
                    drop(tx_img);
                });
                drop(tx_error_worker);
            });
        }
        drop(tx_error);

        // If an error message comes in on this channel, print it and flag if flagmasks provided.
        scope.spawn(move |_| {
            // let mut baseline_flagmasks = baseline_flagmasks_arc.clone();
            for (img_coarse_chan_idx, img_timestep_idx, err) in rx_error.iter() {
                warn!("could not read hdu {:?}", err);

                baseline_flagmasks.iter_mut().for_each(|flagmask| {
                    let flag_stride = flagmask.HorizontalStride();
                    let flag_buf: &mut [bool] = flagmask.pin_mut().BufferMut();
                    let freq_low = fine_chans_per_coarse * img_coarse_chan_idx;
                    for fine_chan in 0..fine_chans_per_coarse {
                        let freq_idx = freq_low + fine_chan;
                        flag_buf[freq_idx * flag_stride + img_timestep_idx] = true;
                    }
                });
            }
        });

        // create a consumer thread for each complex polarization
        for (pol_idx, (tx_img, rx_img)) in pol_img_queues.into_iter().enumerate() {
            // This ensures that
            drop(tx_img);
            let rx_img_worker = rx_img.clone();
            let baseline_imgsets_worker = baseline_imgsets_arc.clone();
            let write_progress_worker = write_progress_arc.clone();
            let total_progess_worker = total_progress_arc.clone();
            scope.spawn(move |_| {
                // The thread consumes an image from it's img queue
                for (img_coarse_chan_idx, img_timestep_idx, img_buf) in rx_img_worker.iter() {
                    img_buf
                        .chunks_exact(floats_per_baseline)
                        .zip(Arc::as_ref(&baseline_imgsets_worker).iter())
                        .for_each(|(baseline_chunk, imgset)| {
                            // The thread writes visibilities from each baseline
                            // into the image buffer from that baseline's
                            // imageset corresponding to this thread's
                            // complex polarization. This ensures that no two
                            // threads access the same image buffer at the same
                            // time.
                            let img_stride = imgset.HorizontalStride();
                            let imgset_buf = unsafe { imgset.ImageBufferMutUnsafe(pol_idx) };
                            baseline_chunk
                                .chunks_exact(floats_per_finechan)
                                .enumerate()
                                .for_each(|(fine_chan_idx, fine_chan_chunk)| {
                                    let img_x = img_timestep_idx;
                                    let img_y =
                                        fine_chans_per_coarse * img_coarse_chan_idx + fine_chan_idx;
                                    unsafe {
                                        *imgset_buf.get_unchecked_mut(img_y * img_stride + img_x) =
                                            *fine_chan_chunk.get_unchecked(pol_idx)
                                    };
                                });
                        });
                    write_progress_worker[pol_idx].inc(1);
                    total_progess_worker.inc(1);
                }
                write_progress_worker[pol_idx].finish_and_clear();
            });
        }
        total_progress_arc.finish();
    })
    .unwrap();

    trace!("end context_to_baseline_imgsets");

    baseline_imgsets
}

#[cfg(test)]
mod tests {
    // TODO: Why does clippy think CxxImageSet.ImageBuffer() is &[f64]?
    #![allow(clippy::float_cmp)]

    use super::{context_to_baseline_imgsets, get_flaggable_timesteps, init_baseline_flagmasks};
    use crate::cxx_aoflagger::ffi::cxx_aoflagger_new;
    use float_cmp::{approx_eq, F32Margin};
    use mwalib::CorrelatorContext;

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
    ///
    /// | gpubox \ timestep | 0      | 1      | 2      | 3      | 4      |
    /// | ----------------- | ------ | ------ | ------ | ------ | ------ |
    /// | 00                | (0, 0) | (0, 1) | .      | (1, 0) | .      |
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

    macro_rules! test_imgset_val {
        ($imgset:expr, $imgset_idx:expr, $fr:expr, $ts:expr, $val:expr) => {
            let left = $imgset.ImageBuffer($imgset_idx)[$fr * $imgset.HorizontalStride() + $ts];
            let right = $val as f32;
            assert!(
                approx_eq!(f32, left, right, F32Margin::default()),
                "{} (0x{:08x}) != {} (0x{:08x})",
                left,
                left as i32,
                right,
                right as i32
            )
        };
    }

    macro_rules! test_flagmask_val {
        ($flagmask:expr, $fr:expr, $ts:expr, $val:expr) => {
            let result = $flagmask.Buffer()[$fr * $flagmask.HorizontalStride() + $ts];
            assert!(
                result == $val,
                "flag(fr={},ts={}) {} != {}",
                $fr,
                $ts,
                result,
                $val,
            )
        };
    }

    #[test]
    fn test_context_to_baseline_imgsets_mwax_flags_missing_hdus() {
        let context = get_mwa_ord_dodgy_context();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(img_timestep_idxs.len(), 4);
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;

        let aoflagger = unsafe { cxx_aoflagger_new() };

        let mut baseline_flagmasks = init_baseline_flagmasks(
            &aoflagger,
            &context,
            &img_coarse_chan_idxs,
            &img_timestep_idxs,
            None,
        );

        let baseline_imgsets = context_to_baseline_imgsets(
            &aoflagger,
            &context,
            &img_coarse_chan_idxs,
            &img_timestep_idxs,
            Some(&mut baseline_flagmasks),
        );

        test_imgset_val!(baseline_imgsets[0], 0, 0, 0, 0x10c5be);
        test_imgset_val!(baseline_imgsets[0], 0, 0, 1, 0x14c5be);
        test_imgset_val!(baseline_imgsets[0], 0, 0, 2, 0x18c5be);
        test_imgset_val!(baseline_imgsets[0], 0, 0, 3, 0x1cc5be);

        test_imgset_val!(baseline_imgsets[0], 1, 0, 0, 0x10c5bf);
        test_imgset_val!(baseline_imgsets[0], 1, 0, 1, 0x14c5bf);
        test_imgset_val!(baseline_imgsets[0], 1, 0, 2, 0x18c5bf);
        test_imgset_val!(baseline_imgsets[0], 1, 0, 3, 0x1cc5bf);

        test_imgset_val!(baseline_imgsets[0], 2, 0, 0, 0x10c5ae);
        test_imgset_val!(baseline_imgsets[0], 2, 0, 1, 0x14c5ae);
        test_imgset_val!(baseline_imgsets[0], 2, 0, 2, 0x18c5ae);
        test_imgset_val!(baseline_imgsets[0], 2, 0, 3, 0x1cc5ae);

        test_imgset_val!(baseline_imgsets[0], 3, 0, 0, -0x10c5af);
        test_imgset_val!(baseline_imgsets[0], 3, 0, 1, -0x14c5af);
        test_imgset_val!(baseline_imgsets[0], 3, 0, 2, -0x18c5af);
        test_imgset_val!(baseline_imgsets[0], 3, 0, 3, -0x1cc5af);

        test_imgset_val!(baseline_imgsets[0], 4, 0, 0, 0x10c5ae);
        test_imgset_val!(baseline_imgsets[0], 4, 0, 1, 0x14c5ae);
        test_imgset_val!(baseline_imgsets[0], 4, 0, 2, 0x18c5ae);
        test_imgset_val!(baseline_imgsets[0], 4, 0, 3, 0x1cc5ae);

        test_imgset_val!(baseline_imgsets[0], 5, 0, 0, 0x10c5af);
        test_imgset_val!(baseline_imgsets[0], 5, 0, 1, 0x14c5af);
        test_imgset_val!(baseline_imgsets[0], 5, 0, 2, 0x18c5af);
        test_imgset_val!(baseline_imgsets[0], 5, 0, 3, 0x1cc5af);

        test_imgset_val!(baseline_imgsets[0], 6, 0, 0, 0x10bec6);
        test_imgset_val!(baseline_imgsets[0], 6, 0, 1, 0x14bec6);
        test_imgset_val!(baseline_imgsets[0], 6, 0, 2, 0x18bec6);
        test_imgset_val!(baseline_imgsets[0], 6, 0, 3, 0x1cbec6);

        test_imgset_val!(baseline_imgsets[0], 7, 0, 0, 0x10bec7);
        test_imgset_val!(baseline_imgsets[0], 7, 0, 1, 0x14bec7);
        test_imgset_val!(baseline_imgsets[0], 7, 0, 2, 0x18bec7);
        test_imgset_val!(baseline_imgsets[0], 7, 0, 3, 0x1cbec7);

        /* ... */

        test_imgset_val!(baseline_imgsets[5], 0, 0, 0, 0x10f1ce);
        test_imgset_val!(baseline_imgsets[5], 0, 0, 1, 0x14f1ce);
        test_imgset_val!(baseline_imgsets[5], 0, 0, 2, 0x18f1ce);
        test_imgset_val!(baseline_imgsets[5], 0, 0, 3, 0x1cf1ce);

        test_imgset_val!(baseline_imgsets[5], 1, 0, 0, -0x10f1cf);
        test_imgset_val!(baseline_imgsets[5], 1, 0, 1, -0x14f1cf);
        test_imgset_val!(baseline_imgsets[5], 1, 0, 2, -0x18f1cf);
        test_imgset_val!(baseline_imgsets[5], 1, 0, 3, -0x1cf1cf);

        test_imgset_val!(baseline_imgsets[5], 2, 0, 0, 0x10ea26);
        test_imgset_val!(baseline_imgsets[5], 2, 0, 1, 0x14ea26);
        test_imgset_val!(baseline_imgsets[5], 2, 0, 2, 0x18ea26);
        test_imgset_val!(baseline_imgsets[5], 2, 0, 3, 0x1cea26);

        test_imgset_val!(baseline_imgsets[5], 3, 0, 0, -0x10ea27);
        test_imgset_val!(baseline_imgsets[5], 3, 0, 1, -0x14ea27);
        test_imgset_val!(baseline_imgsets[5], 3, 0, 2, -0x18ea27);
        test_imgset_val!(baseline_imgsets[5], 3, 0, 3, -0x1cea27);

        test_imgset_val!(baseline_imgsets[5], 4, 0, 0, 0x10f1be);
        test_imgset_val!(baseline_imgsets[5], 4, 0, 1, 0x14f1be);
        test_imgset_val!(baseline_imgsets[5], 4, 0, 2, 0x18f1be);
        test_imgset_val!(baseline_imgsets[5], 4, 0, 3, 0x1cf1be);

        test_imgset_val!(baseline_imgsets[5], 5, 0, 0, -0x10f1bf);
        test_imgset_val!(baseline_imgsets[5], 5, 0, 1, -0x14f1bf);
        test_imgset_val!(baseline_imgsets[5], 5, 0, 2, -0x18f1bf);
        test_imgset_val!(baseline_imgsets[5], 5, 0, 3, -0x1cf1bf);

        test_imgset_val!(baseline_imgsets[5], 6, 0, 0, 0x10ea16);
        test_imgset_val!(baseline_imgsets[5], 6, 0, 1, 0x14ea16);
        test_imgset_val!(baseline_imgsets[5], 6, 0, 2, 0x18ea16);
        test_imgset_val!(baseline_imgsets[5], 6, 0, 3, 0x1cea16);

        test_imgset_val!(baseline_imgsets[5], 7, 0, 0, -0x10ea17);
        test_imgset_val!(baseline_imgsets[5], 7, 0, 1, -0x14ea17);
        test_imgset_val!(baseline_imgsets[5], 7, 0, 2, -0x18ea17);
        test_imgset_val!(baseline_imgsets[5], 7, 0, 3, -0x1cea17);

        /* ... */

        test_imgset_val!(baseline_imgsets[0], 0, 2, 0, 0x04c5be); //0x00c5be);
        test_imgset_val!(baseline_imgsets[0], 0, 2, 1, 0x0); //0x04c5be);
        test_imgset_val!(baseline_imgsets[0], 0, 2, 2, 0x08c5be);
        test_imgset_val!(baseline_imgsets[0], 0, 2, 3, 0x0); //0x0cc5be);

        test_flagmask_val!(baseline_flagmasks[0], 2, 0, false);
        test_flagmask_val!(baseline_flagmasks[0], 2, 1, true);
        test_flagmask_val!(baseline_flagmasks[0], 2, 2, false);
        test_flagmask_val!(baseline_flagmasks[0], 2, 3, true);

        test_flagmask_val!(baseline_flagmasks[0], 1, 0, false);
        test_flagmask_val!(baseline_flagmasks[0], 1, 1, false);
        test_flagmask_val!(baseline_flagmasks[0], 1, 2, false);
        test_flagmask_val!(baseline_flagmasks[0], 1, 3, false);
    }

    #[test]
    fn test_context_to_baseline_imgsets_mwax() {
        let context = get_mwax_context();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(img_timestep_idxs.len(), 4);
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;

        let aoflagger = unsafe { cxx_aoflagger_new() };

        let baseline_imgsets = context_to_baseline_imgsets(
            &aoflagger,
            &context,
            &img_coarse_chan_idxs,
            &img_timestep_idxs,
            None,
        );

        test_imgset_val!(baseline_imgsets[0], 0, 0, 0, 0x410000);
        test_imgset_val!(baseline_imgsets[0], 0, 0, 1, 0x410100);
        test_imgset_val!(baseline_imgsets[0], 0, 0, 2, 0x410200);
        test_imgset_val!(baseline_imgsets[0], 0, 0, 3, 0x410300);

        test_imgset_val!(baseline_imgsets[0], 0, 1, 0, 0x410008);
        test_imgset_val!(baseline_imgsets[0], 0, 1, 1, 0x410108);
        test_imgset_val!(baseline_imgsets[0], 0, 1, 2, 0x410208);
        test_imgset_val!(baseline_imgsets[0], 0, 1, 3, 0x410308);

        test_imgset_val!(baseline_imgsets[0], 0, 2, 0, 0x410400);
        test_imgset_val!(baseline_imgsets[0], 0, 2, 1, 0x410500);
        test_imgset_val!(baseline_imgsets[0], 0, 2, 2, 0x410600);
        test_imgset_val!(baseline_imgsets[0], 0, 2, 3, 0x410700);

        test_imgset_val!(baseline_imgsets[0], 0, 3, 0, 0x410408);
        test_imgset_val!(baseline_imgsets[0], 0, 3, 1, 0x410508);
        test_imgset_val!(baseline_imgsets[0], 0, 3, 2, 0x410608);
        test_imgset_val!(baseline_imgsets[0], 0, 3, 3, 0x410708);

        test_imgset_val!(baseline_imgsets[0], 1, 0, 0, 0x410001);
        test_imgset_val!(baseline_imgsets[0], 1, 0, 1, 0x410101);
        test_imgset_val!(baseline_imgsets[0], 1, 0, 2, 0x410201);
        test_imgset_val!(baseline_imgsets[0], 1, 0, 3, 0x410301);

        /* ... */

        test_imgset_val!(baseline_imgsets[0], 7, 0, 0, 0x410007);
        test_imgset_val!(baseline_imgsets[0], 7, 0, 1, 0x410107);
        test_imgset_val!(baseline_imgsets[0], 7, 0, 2, 0x410207);
        test_imgset_val!(baseline_imgsets[0], 7, 0, 3, 0x410307);

        /* ... */

        test_imgset_val!(baseline_imgsets[2], 0, 0, 0, 0x410020);
        test_imgset_val!(baseline_imgsets[2], 0, 0, 1, 0x410120);
        test_imgset_val!(baseline_imgsets[2], 0, 0, 2, 0x410220);
        test_imgset_val!(baseline_imgsets[2], 0, 0, 3, 0x410320);

        test_imgset_val!(baseline_imgsets[2], 0, 1, 0, 0x410028);
        test_imgset_val!(baseline_imgsets[2], 0, 1, 1, 0x410128);
        test_imgset_val!(baseline_imgsets[2], 0, 1, 2, 0x410228);
        test_imgset_val!(baseline_imgsets[2], 0, 1, 3, 0x410328);
    }

    #[test]
    fn test_context_to_baseline_imgsets_mwa_ord() {
        let context = get_mwa_ord_context();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(img_timestep_idxs.len(), 4);
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;

        let aoflagger = unsafe { cxx_aoflagger_new() };

        let baseline_imgsets = context_to_baseline_imgsets(
            &aoflagger,
            &context,
            &img_coarse_chan_idxs,
            &img_timestep_idxs,
            None,
        );

        test_imgset_val!(baseline_imgsets[0], 0, 0, 0, 0x10c5be);
        test_imgset_val!(baseline_imgsets[0], 0, 0, 1, 0x14c5be);
        test_imgset_val!(baseline_imgsets[0], 0, 0, 2, 0x18c5be);
        test_imgset_val!(baseline_imgsets[0], 0, 0, 3, 0x1cc5be);

        test_imgset_val!(baseline_imgsets[0], 1, 0, 0, 0x10c5bf);
        test_imgset_val!(baseline_imgsets[0], 1, 0, 1, 0x14c5bf);
        test_imgset_val!(baseline_imgsets[0], 1, 0, 2, 0x18c5bf);
        test_imgset_val!(baseline_imgsets[0], 1, 0, 3, 0x1cc5bf);

        test_imgset_val!(baseline_imgsets[0], 2, 0, 0, 0x10c5ae);
        test_imgset_val!(baseline_imgsets[0], 2, 0, 1, 0x14c5ae);
        test_imgset_val!(baseline_imgsets[0], 2, 0, 2, 0x18c5ae);
        test_imgset_val!(baseline_imgsets[0], 2, 0, 3, 0x1cc5ae);

        test_imgset_val!(baseline_imgsets[0], 3, 0, 0, -0x10c5af);
        test_imgset_val!(baseline_imgsets[0], 3, 0, 1, -0x14c5af);
        test_imgset_val!(baseline_imgsets[0], 3, 0, 2, -0x18c5af);
        test_imgset_val!(baseline_imgsets[0], 3, 0, 3, -0x1cc5af);

        test_imgset_val!(baseline_imgsets[0], 4, 0, 0, 0x10c5ae);
        test_imgset_val!(baseline_imgsets[0], 4, 0, 1, 0x14c5ae);
        test_imgset_val!(baseline_imgsets[0], 4, 0, 2, 0x18c5ae);
        test_imgset_val!(baseline_imgsets[0], 4, 0, 3, 0x1cc5ae);

        test_imgset_val!(baseline_imgsets[0], 5, 0, 0, 0x10c5af);
        test_imgset_val!(baseline_imgsets[0], 5, 0, 1, 0x14c5af);
        test_imgset_val!(baseline_imgsets[0], 5, 0, 2, 0x18c5af);
        test_imgset_val!(baseline_imgsets[0], 5, 0, 3, 0x1cc5af);

        test_imgset_val!(baseline_imgsets[0], 6, 0, 0, 0x10bec6);
        test_imgset_val!(baseline_imgsets[0], 6, 0, 1, 0x14bec6);
        test_imgset_val!(baseline_imgsets[0], 6, 0, 2, 0x18bec6);
        test_imgset_val!(baseline_imgsets[0], 6, 0, 3, 0x1cbec6);

        test_imgset_val!(baseline_imgsets[0], 7, 0, 0, 0x10bec7);
        test_imgset_val!(baseline_imgsets[0], 7, 0, 1, 0x14bec7);
        test_imgset_val!(baseline_imgsets[0], 7, 0, 2, 0x18bec7);
        test_imgset_val!(baseline_imgsets[0], 7, 0, 3, 0x1cbec7);

        /* ... */

        test_imgset_val!(baseline_imgsets[5], 0, 0, 0, 0x10f1ce);
        test_imgset_val!(baseline_imgsets[5], 0, 0, 1, 0x14f1ce);
        test_imgset_val!(baseline_imgsets[5], 0, 0, 2, 0x18f1ce);
        test_imgset_val!(baseline_imgsets[5], 0, 0, 3, 0x1cf1ce);

        test_imgset_val!(baseline_imgsets[5], 1, 0, 0, -0x10f1cf);
        test_imgset_val!(baseline_imgsets[5], 1, 0, 1, -0x14f1cf);
        test_imgset_val!(baseline_imgsets[5], 1, 0, 2, -0x18f1cf);
        test_imgset_val!(baseline_imgsets[5], 1, 0, 3, -0x1cf1cf);

        test_imgset_val!(baseline_imgsets[5], 2, 0, 0, 0x10ea26);
        test_imgset_val!(baseline_imgsets[5], 2, 0, 1, 0x14ea26);
        test_imgset_val!(baseline_imgsets[5], 2, 0, 2, 0x18ea26);
        test_imgset_val!(baseline_imgsets[5], 2, 0, 3, 0x1cea26);

        test_imgset_val!(baseline_imgsets[5], 3, 0, 0, -0x10ea27);
        test_imgset_val!(baseline_imgsets[5], 3, 0, 1, -0x14ea27);
        test_imgset_val!(baseline_imgsets[5], 3, 0, 2, -0x18ea27);
        test_imgset_val!(baseline_imgsets[5], 3, 0, 3, -0x1cea27);

        test_imgset_val!(baseline_imgsets[5], 4, 0, 0, 0x10f1be);
        test_imgset_val!(baseline_imgsets[5], 4, 0, 1, 0x14f1be);
        test_imgset_val!(baseline_imgsets[5], 4, 0, 2, 0x18f1be);
        test_imgset_val!(baseline_imgsets[5], 4, 0, 3, 0x1cf1be);

        test_imgset_val!(baseline_imgsets[5], 5, 0, 0, -0x10f1bf);
        test_imgset_val!(baseline_imgsets[5], 5, 0, 1, -0x14f1bf);
        test_imgset_val!(baseline_imgsets[5], 5, 0, 2, -0x18f1bf);
        test_imgset_val!(baseline_imgsets[5], 5, 0, 3, -0x1cf1bf);

        test_imgset_val!(baseline_imgsets[5], 6, 0, 0, 0x10ea16);
        test_imgset_val!(baseline_imgsets[5], 6, 0, 1, 0x14ea16);
        test_imgset_val!(baseline_imgsets[5], 6, 0, 2, 0x18ea16);
        test_imgset_val!(baseline_imgsets[5], 6, 0, 3, 0x1cea16);

        test_imgset_val!(baseline_imgsets[5], 7, 0, 0, -0x10ea17);
        test_imgset_val!(baseline_imgsets[5], 7, 0, 1, -0x14ea17);
        test_imgset_val!(baseline_imgsets[5], 7, 0, 2, -0x18ea17);
        test_imgset_val!(baseline_imgsets[5], 7, 0, 3, -0x1cea17);

        /* ... */

        test_imgset_val!(baseline_imgsets[0], 0, 2, 0, 0x00c5be);
        test_imgset_val!(baseline_imgsets[0], 0, 2, 1, 0x04c5be);
        test_imgset_val!(baseline_imgsets[0], 0, 2, 2, 0x08c5be);
        test_imgset_val!(baseline_imgsets[0], 0, 2, 3, 0x0cc5be);
    }
}
