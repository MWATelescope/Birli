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
use error::{
    BirliError,
    BirliError::{NoCommonTimesteps, NoProvidedTimesteps},
};
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

use mwalib::CorrelatorContext;
// use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::os::raw::c_short;
use std::sync::Arc;

use itertools::izip;

pub mod flag_io;
use flag_io::FlagFileSet;

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

const SPEED_OF_LIGHT_IN_VACUUM: f64 = 299792458.0; // speed of light in m/s

fn _correct_cable_length_buffers_cotter(
    freq_hz: &u32,
    electrical_length_m: &f64,
    buf_re: &mut [f32],
    buf_im: &mut [f32],
) {
    let angle: f64 = -2.0 * PI * electrical_length_m * (*freq_hz as f64) / SPEED_OF_LIGHT_IN_VACUUM;
    let (sin_angle_f64, cos_angle_f64) = angle.sin_cos();
    let (sin_angle, cos_angle) = (sin_angle_f64 as f32, cos_angle_f64 as f32);

    izip!(buf_re.iter_mut(), buf_im.iter_mut()).for_each(|(re, im)| {
        let vis_re = *re;
        *re = cos_angle * vis_re - sin_angle * *im;
        *im = sin_angle * vis_re + cos_angle * *im;
    })
}

/// Perform cable length corrections, given an observation's
/// [`mwalib::CorrelatorContext`] and a vector of [`CxxImageSet`]s for  each
/// baseline.
///
/// Cable lengths are determined by the difference between a baseline's rfInput
/// electrical lengths in the metafits. Complex visibilities are phase-shifted
/// by an angle determined by the electrical length, and the channel's
/// frequency.
///
/// # Examples
///
/// ```rust
/// use birli::{context_to_baseline_imgsets, cxx_aoflagger_new, correct_cable_lengths};
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
/// // create a CxxAOFlagger object to perform AOFlagger operations
/// let aoflagger = unsafe { cxx_aoflagger_new() };
///
/// // generate imagesets for each baseline in the format required by aoflagger
/// let mut baseline_imgsets = context_to_baseline_imgsets(
///     &aoflagger,
///     &context,
///     &context.common_coarse_chan_indices.clone(),
///     &context.common_timestep_indices.clone(),
///     None,
/// );
///
/// correct_cable_lengths(&context, &mut baseline_imgsets);
/// ```
///
/// # Accuracy
///
/// This follows the Cotter implementation of cable length correction, however
/// there is a slower but more accurate version of the calculation which
/// uses f64 values for the sin_cos. According to benchmarks, the Cotter
/// implementation is about 32% faster (5.9 seconds vs 8.6) than the more
/// precise implementation, and they vary by about three parts in four million.
/// Therefore it was decided that the Cotter implementation was more favourable.
///
/// Here is the more accurate implementation of
/// [`_correct_cable_length_buffers_cotter`]
///
/// ```rust
/// use itertools::izip;
/// use std::f64::consts::PI;
///
/// const SPEED_OF_LIGHT_IN_VACUUM: f64 = 299792458.0; // speed of light in m/s
///
/// fn _correct_cable_length_buffers_precise(
///     freq_hz: &u32,
///     electrical_length_m: &f64,
///     buf_re: &mut [f32],
///     buf_im: &mut [f32],
/// ) {
///     let angle: f64 = -2.0 * PI * electrical_length_m * (*freq_hz as f64) / SPEED_OF_LIGHT_IN_VACUUM;
///     let (sin_angle, cos_angle) = angle.sin_cos();
///     izip!(buf_re.iter_mut(), buf_im.iter_mut()).for_each(|(re, im)| {
///         let vis_re = *re as f64;
///         let vis_im = *im as f64;
///         *re = (cos_angle * vis_re - sin_angle * vis_im) as f32;
///         *im = (sin_angle * vis_re + cos_angle * vis_im) as f32;
///     })
/// }
/// ```

pub fn correct_cable_lengths(
    context: &CorrelatorContext,
    baseline_imgsets: &mut Vec<UniquePtr<CxxImageSet>>,
) {
    trace!("start correct_cable_lengths");

    let baselines = &context.metafits_context.baselines;
    let antennas = &context.metafits_context.antennas;
    let coarse_chans = &context.coarse_chans;

    let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
    let fine_chan_width_hz = &context.metafits_context.corr_fine_chan_width_hz;

    // A vector of all fine channel frequencies for all coarse channels in ascending order.
    let all_freqs_hz: Vec<u32> = coarse_chans
        .iter()
        .flat_map(|coarse_chan| {
            let chan_start_hz = coarse_chan.chan_start_hz;
            (0..fine_chans_per_coarse).map(move |fine_chan_idx| {
                chan_start_hz + (fine_chan_idx as u32 * fine_chan_width_hz)
            })
        })
        .collect();

    izip!(baseline_imgsets.iter_mut(), baselines).for_each(|(imgset, baseline)| {
        let imgset_stride: usize = imgset.HorizontalStride();
        let ant1 = &antennas[baseline.ant1_index];
        let ant2 = &antennas[baseline.ant2_index];

        // TODO: skip if baseline.ant1_index == baseline.ant2_index ?

        let pol_lengths = vec![
            ant2.rfinput_x.electrical_length_m - ant1.rfinput_x.electrical_length_m,
            ant2.rfinput_y.electrical_length_m - ant1.rfinput_x.electrical_length_m,
            ant2.rfinput_x.electrical_length_m - ant1.rfinput_y.electrical_length_m,
            ant2.rfinput_y.electrical_length_m - ant1.rfinput_y.electrical_length_m,
        ];
        pol_lengths
            .iter()
            .enumerate()
            .for_each(|(pol_idx, electrical_length_m)| {
                let imgset_buf_re: &mut [f32] = unsafe { imgset.ImageBufferMutUnsafe(2 * pol_idx) };
                let imgset_buf_im: &mut [f32] =
                    unsafe { imgset.ImageBufferMutUnsafe(2 * pol_idx + 1) };

                izip!(
                    all_freqs_hz.iter(),
                    imgset_buf_re.chunks_mut(imgset_stride),
                    imgset_buf_im.chunks_mut(imgset_stride)
                )
                .for_each(|(freq_hz, imgset_chunk_re, imgset_chunk_im)| {
                    // _correct_cable_length_buffers_precise(
                    _correct_cable_length_buffers_cotter(
                        freq_hz,
                        electrical_length_m,
                        imgset_chunk_re,
                        imgset_chunk_im,
                    )
                });
            });
    });
    trace!("end correct_cable_lengths");
}

/// Produce a vector of timesteps which can be used for creating imagesets and
/// flagmasks for aoflagger flagging, given an [`mwalib::CorrelatorContext`].
///
/// These timesteps:
/// - are contiguous, and are each separated by an integration time
/// - start at the latest start time of visibilities in all provided gpubox files
/// - end at the latest provided visibility
///
/// Start times are determined from the TIME and MILLITIM headers of individual
/// gpubox visibility HDUs
///
/// [`mwalib::CorrelatorContext`]: https://docs.rs/mwalib/latest/mwalib/struct.CorrelatorContext.html
///
/// # Examples
///
/// ```rust
/// use birli::{get_flaggable_timesteps};
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
/// let context = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
///
/// let flaggable_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
///
/// assert_eq!(flaggable_timestep_idxs.len(), 4);
/// ```
///
/// # Errors
///
/// This function uses [`mwalib::CorrelatorContext.common_timestep_indices`] and
/// [`mwalib::CorrelatorContext.provided_timestep_indices`]. It will return
/// [`errors::BirliError::NoProvidedTimesteps`] if mwalib has determined that no
/// timesteps have been provided, and will return
/// [`errors::BirliError::NoCommonTimesteps`] if mwalib has determined
///
/// [`mwalib::CorrelatorContext.common_timestep_indices`]: https://docs.rs/mwalib/latest/mwalib/struct.CorrelatorContext.html#structfield.common_timestep_indices
/// [`mwalib::CorrelatorContext.provided_timestep_indices`]: https://docs.rs/mwalib/latest/mwalib/struct.CorrelatorContext.html#structfield.provided_timestep_indices
///
///
pub fn get_flaggable_timesteps(context: &CorrelatorContext) -> Result<Vec<usize>, BirliError> {
    match (
        context.common_timestep_indices.first(),
        context.provided_timestep_indices.last(),
    ) {
        (Some(&first), Some(&last)) => Ok((first..last + 1).collect()),
        (.., None) => Err(NoProvidedTimesteps {
            timestep_info: format!("{:?}", &context.gpubox_time_map),
        }),
        (None, ..) => Err(NoCommonTimesteps {
            timestep_info: format!("{:?}", &context.gpubox_time_map),
        }),
    }
}

/// Produce a vector of flags for antennas whose imputs are flagged in `context`
///
/// # Examples
///
/// ```rust
/// use birli::{get_antenna_flags, cxx_aoflagger_new};
/// use mwalib::CorrelatorContext;
///
/// // define our input files
/// let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
/// let gpufits_paths = vec![
///     "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits",
///     "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits",
///     "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits",
///     "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits",
/// ];
///
/// // Create an mwalib::CorrelatorContext for accessing visibilities.
/// let context = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
/// let antennta_flags = get_antenna_flags(&context);
///
/// assert!(!antennta_flags[0]);
/// assert!(antennta_flags[63]);
/// assert!(antennta_flags[62]);
/// assert!(antennta_flags[11]);
/// assert!(antennta_flags[111]);
/// assert!(antennta_flags[91]);
/// assert!(antennta_flags[95]);
/// assert!(antennta_flags[93]);
/// assert!(antennta_flags[80]);
/// assert!(antennta_flags[87]);
/// ```
pub fn get_antenna_flags(context: &CorrelatorContext) -> Vec<bool> {
    context
        .metafits_context
        .antennas
        .iter()
        .map(|antenna| antenna.rfinput_x.flagged || antenna.rfinput_y.flagged)
        .collect()
}

/// Produce a vector of flags for baslines whose antennae are flagged in `antenna_flags`
fn get_baseline_flags(context: &CorrelatorContext, antenna_flags: Vec<bool>) -> Vec<bool> {
    assert_eq!(antenna_flags.len(), context.metafits_context.num_ants);
    context
        .metafits_context
        .baselines
        .clone()
        .into_iter()
        .map(|baseline| {
            antenna_flags[baseline.ant1_index] || antenna_flags[baseline.ant2_index]
            // *antenna_flags.get(baseline.ant1_index).unwrap_or(&false) || *antenna_flags.get(baseline.ant2_index).unwrap_or(&false)
        })
        .collect()
}

/// Initialize a vector of [`CxxFlagMask`]s for each baseline and provided
/// coarse channel and timestep indices
///
/// TODO: use:
/// - coarse_chan_flags
/// - fine_chan_flags
/// - timestep_flags
/// to set correlator mask
pub fn init_baseline_flagmasks(
    aoflagger: &CxxAOFlagger,
    context: &CorrelatorContext,
    // num_baselines: usize,
    // fine_chans_per_coarse: usize,
    img_coarse_chan_idxs: &[usize],
    img_timestep_idxs: &[usize],
    antenna_flags: Option<Vec<bool>>,
    // coarse_chan_flags: Option<Vec<bool>>,
    // fine_chan_flags: Option<Vec<bool>>,
    // timestep_flags: Option<Vec<bool>>,
    // flag_bad_inputs: bool,
) -> Vec<UniquePtr<CxxFlagMask>> {
    let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
    let width = img_timestep_idxs.len();
    let height = img_coarse_chan_idxs.len() * fine_chans_per_coarse;

    // create correlator mask
    let correlator_mask = unsafe { aoflagger.MakeFlagMask(width, height, false) };

    let antenna_flags = match antenna_flags {
        Some(antenna_flags) => antenna_flags,
        None => (0..context.metafits_context.num_ants)
            .map(|_| false)
            .collect(),
    };

    get_baseline_flags(&context, antenna_flags)
        .iter()
        .map(|&baseline_flag| {
            if baseline_flag {
                unsafe { aoflagger.MakeFlagMask(width, height, true) }
            } else {
                let mut flag_mask = unsafe { aoflagger.MakeFlagMask(width, height, false) };
                flagmask_set(&mut flag_mask, &correlator_mask);
                flag_mask
            }
        })
        .collect()
}

/// Flag an observation's visibilities, given a [`CxxAOFlagger`] instance, a [`CxxStrategy`]
/// filename, and a vector of [`CxxImageSet`]s for each baseline in the observation returning a
/// vector of [`CxxFlagMask`]s.
///
/// See: [`flag_imgsets_existing`] for working with existing flagmasks
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
/// // create a CxxAOFlagger object to perform AOFlagger operations
/// let aoflagger = unsafe { cxx_aoflagger_new() };
///
/// // generate imagesets for each baseline in the format required by aoflagger
/// let baseline_imgsets = context_to_baseline_imgsets(
///     &aoflagger,
///     &context,
///     &context.common_coarse_chan_indices.clone(),
///     &context.common_timestep_indices.clone(),
///     None,
/// );
///
/// // use the default strategy file location for MWA
/// let strategy_filename = &aoflagger.FindStrategyFileMWA();
///
/// // run the strategy on the imagesets, and get the resulting flagmasks for each baseline
/// let baseline_flagmasks = flag_imgsets(&aoflagger, &strategy_filename, baseline_imgsets);
/// ```
pub fn flag_imgsets(
    aoflagger: &CxxAOFlagger,
    strategy_filename: &str,
    baseline_imgsets: Vec<UniquePtr<CxxImageSet>>,
) -> Vec<UniquePtr<CxxFlagMask>> {
    trace!("start flag_imgsets");

    let flag_progress = ProgressBar::new(baseline_imgsets.len() as u64);
    flag_progress.set_style(
        ProgressStyle::default_bar()
            .template(
                "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
            )
            .progress_chars("=> "),
    );
    flag_progress.set_message("flagging b.lines");

    let baseline_flagmasks = baseline_imgsets
        .par_iter()
        .map(|imgset| {
            let flagmask = aoflagger
                .LoadStrategyFile(&strategy_filename.to_string())
                .Run(&imgset);
            flag_progress.inc(1);
            flagmask
        })
        .collect();
    flag_progress.finish();
    trace!("end flag_imgsets");
    baseline_flagmasks
}

/// Perform the binary or operation on the flag buffer of `this_flagmask` with
/// the buffer of `other_flagmask`, storing the result in `this_flagmask`
pub fn flagmask_or(
    this_flagmask: &mut UniquePtr<CxxFlagMask>,
    other_flagmask: &UniquePtr<CxxFlagMask>,
) {
    assert_eq!(
        this_flagmask.Width(),
        other_flagmask.Width(),
        "flagmasks should have the same width."
    );
    assert_eq!(
        this_flagmask.Height(),
        other_flagmask.Height(),
        "flagmasks should have the same height."
    );
    let this_buffer: &mut [bool] = this_flagmask.pin_mut().BufferMut();
    let other_buffer: &[bool] = other_flagmask.Buffer();

    this_buffer
        .iter_mut()
        .zip(other_buffer.iter())
        .for_each(|(this_flag, other_flag)| *this_flag |= other_flag);
}

/// Set the flag buffer of `this_flagmask` from the buffer of `other_flagmask`,
/// storing the result in `this_flagmask`
pub fn flagmask_set(
    this_flagmask: &mut UniquePtr<CxxFlagMask>,
    other_flagmask: &UniquePtr<CxxFlagMask>,
) {
    assert_eq!(
        this_flagmask.Width(),
        other_flagmask.Width(),
        "flagmasks should have the same width."
    );
    assert_eq!(
        this_flagmask.Height(),
        other_flagmask.Height(),
        "flagmasks should have the same height."
    );
    let this_buffer: &mut [bool] = this_flagmask.pin_mut().BufferMut();
    let other_buffer: &[bool] = other_flagmask.Buffer();

    this_buffer
        .iter_mut()
        .zip(other_buffer.iter())
        .for_each(|(this_flag, other_flag)| *this_flag = *other_flag);
}

/// Flag an observation's visibilities, given a [`CxxAOFlagger`] instance, a [`CxxStrategy`]
/// filename, and vectors of [`CxxImageSet`]s and [`CxxFlagMask`]s for each baseline in the
/// observation returning a vector of [`CxxFlagMask`]s.
///
/// # Examples
///
/// ```
/// use birli::{context_to_baseline_imgsets, flag_imgsets, write_flags,
///     cxx_aoflagger_new, flag_imgsets_existing, init_baseline_flagmasks,
///     get_antenna_flags, get_flaggable_timesteps};
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
/// // create a CxxAOFlagger object to perform AOFlagger operations
/// let aoflagger = unsafe { cxx_aoflagger_new() };
///
/// // specify which coarse_chan and timestep indices we want to load into an image.
/// let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
/// let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
///
/// let mut baseline_flagmasks = init_baseline_flagmasks(
///     &aoflagger,
///     &context,
///     &img_coarse_chan_idxs,
///     &img_timestep_idxs,
///     Some(get_antenna_flags(&context)),
/// );
///
/// // generate imagesets for each baseline in the format required by aoflagger
/// let baseline_imgsets = context_to_baseline_imgsets(
///     &aoflagger,
///     &context,
///     &img_coarse_chan_idxs,
///     &img_timestep_idxs,
///     None,
/// );
///
/// // use the default strategy file location for MWA
/// let strategy_filename = &aoflagger.FindStrategyFileMWA();
///
/// // run the strategy on the imagesets, and get the resulting flagmasks for each baseline
/// let baseline_flagmasks = flag_imgsets_existing(
///     &aoflagger,
///     &strategy_filename,
///     &baseline_imgsets,
///     &mut baseline_flagmasks,
///     true
/// );
/// ```
///
/// # TODO: just write back in to existing
pub fn flag_imgsets_existing(
    aoflagger: &CxxAOFlagger,
    strategy_filename: &str,
    baseline_imgsets: &[UniquePtr<CxxImageSet>],
    baseline_flagmasks: &mut Vec<UniquePtr<CxxFlagMask>>,
    re_apply_existing: bool,
) {
    trace!("start flag_imgsets_existing");

    let flag_progress = ProgressBar::new(baseline_imgsets.len() as u64);
    flag_progress.set_style(
        ProgressStyle::default_bar()
            .template(
                "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
            )
            .progress_chars("=> "),
    );
    flag_progress.set_message("flagging b.lines");

    assert_eq!(baseline_imgsets.len(), baseline_flagmasks.len());

    baseline_flagmasks
        .par_iter_mut()
        .zip(baseline_imgsets)
        .for_each(|(flagmask, imgset)| {
            let new_flagmask = aoflagger
                .LoadStrategyFile(&strategy_filename.to_string())
                .RunExisting(&imgset, &flagmask);
            if re_apply_existing {
                flagmask_or(flagmask, &new_flagmask);
            } else {
                flagmask_set(flagmask, &new_flagmask);
            }
            flag_progress.inc(1);
        });

    flag_progress.finish();
    trace!("end flag_imgsets_existing");
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
/// let baseline_imgsets = context_to_baseline_imgsets(
///     &aoflagger,
///     &context,
///     &context.common_coarse_chan_indices.clone(),
///     &context.common_timestep_indices.clone(),
///     None,
/// );
///
/// // use the default strategy file location for MWA
/// let strategy_filename = &aoflagger.FindStrategyFileMWA();
///
/// // run the strategy on the imagesets, and get the resulting flagmasks for each baseline
/// let baseline_flagmasks = flag_imgsets(&aoflagger, &strategy_filename, baseline_imgsets);
///
/// // Get a list of all gpubox IDs
/// let gpubox_ids: Vec<usize> = context
///     .common_coarse_chan_indices
///     .iter()
///     .map(|&chan| context.coarse_chans[chan].gpubox_number)
///     .collect();
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
    trace!("start write_flags");

    // TODO: error handling instead of unwrap.

    let mut flag_file_set =
        FlagFileSet::new(filename_template, &gpubox_ids, context.mwa_version).unwrap();
    flag_file_set
        .write_baseline_flagmasks(&context, baseline_flagmasks)
        .unwrap();

    trace!("end write_flags");
}

#[cfg(test)]
mod tests {
    // TODO: Why does clippy think CxxImageSet.ImageBuffer() is &[f64]?
    #![allow(clippy::float_cmp)]

    use super::{
        context_to_baseline_imgsets, flag_imgsets, flag_imgsets_existing, flagmask_or,
        flagmask_set, get_antenna_flags, get_flaggable_timesteps, init_baseline_flagmasks,
        write_flags, FlagFileSet,
    };
    use crate::{
        correct_cable_lengths,
        cxx_aoflagger::ffi::{cxx_aoflagger_new, CxxAOFlagger, CxxFlagMask, CxxImageSet},
        error::BirliError,
        SPEED_OF_LIGHT_IN_VACUUM,
    };
    use cxx::UniquePtr;
    use float_cmp::{approx_eq, F32Margin};
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

    /// Get a dummy MWA Ord context with no overlapping timesteps
    ///
    /// The gpubox (batch, hdu) tuples look like this:
    ///
    /// | gpubox \ timestep | 0      | 1      | 2      | 3      |
    /// | ----------------- | ------ | ------ | ------ | ------ |
    /// | 00                | (0, 0) | (0, 1) | .      | .      |
    /// | 01                | .      | .      | (0, 1) | (1, 0) |
    fn get_mwa_ord_no_overlap_context() -> CorrelatorContext {
        let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
        let gpufits_paths = vec![
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits",
        ];
        CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap()
    }

    /// Get a dummy MWA Ord context with no timesteps
    fn get_mwa_ord_no_timesteps_context() -> CorrelatorContext {
        // let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
        // let gpufits_paths =
        //     vec!["tests/data/1196175296_mwa_ord/empty/1196175296_20171201145440_gpubox01_00.fits"];
        let mut context = get_mwa_ord_no_overlap_context();
        context.provided_timestep_indices = vec![];
        context
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
    fn test_get_flaggable_timesteps_handles_no_overlap() {
        let context = get_mwa_ord_no_overlap_context();

        assert!(matches!(
            get_flaggable_timesteps(&context),
            Err(BirliError::NoCommonTimesteps { .. })
        ));
    }

    #[test]
    fn test_get_flaggable_timesteps_handles_no_provided() {
        let context = get_mwa_ord_no_timesteps_context();

        assert!(matches!(
            get_flaggable_timesteps(&context),
            Err(BirliError::NoProvidedTimesteps { .. })
        ));
    }

    #[test]
    fn test_get_flaggable_timesteps_handles_dodgy() {
        let context = get_mwa_ord_dodgy_context();

        let timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(timestep_idxs.len(), 4);
        assert_eq!(
            timestep_idxs.first(),
            context.common_timestep_indices.first()
        );
        assert_eq!(
            timestep_idxs.last(),
            context.provided_timestep_indices.last()
        );
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

    #[test]
    fn test_flag_imgsets_minimal() {
        let width = 64;
        let height = 64;
        let count = 8;
        let num_baselines = 2;

        let noise_bl = 1;
        let noise_x = 32;
        let noise_y = 32;
        let noise_z = 1;
        let signal_val = 0 as f32;
        let noise_val = 0xffffff as f32;

        let aoflagger = unsafe { cxx_aoflagger_new() };

        let mut baseline_imgsets: Vec<UniquePtr<CxxImageSet>> = (0..num_baselines)
            .into_iter()
            .map(|_| unsafe { aoflagger.MakeImageSet(width, height, count, signal_val, width) })
            .collect();

        let img_stride = baseline_imgsets[0].HorizontalStride();
        // modify the imageset to add some synthetic noise
        baseline_imgsets[noise_bl].pin_mut().ImageBufferMut(noise_z)
            [noise_y * img_stride + noise_x] = noise_val;

        let strategy_file_minimal = aoflagger.FindStrategyFileGeneric(&String::from("minimal"));

        let baseline_flagmasks = flag_imgsets(&aoflagger, &strategy_file_minimal, baseline_imgsets);

        let flag_stride = baseline_flagmasks[0].HorizontalStride();
        assert!(!baseline_flagmasks[0].Buffer()[0]);
        assert!(!baseline_flagmasks[0].Buffer()[noise_y * flag_stride + noise_x]);
        assert!(!baseline_flagmasks[noise_bl].Buffer()[0]);
        assert!(baseline_flagmasks[noise_bl].Buffer()[noise_y * flag_stride + noise_x]);
    }

    #[test]
    fn test_flag_imgsets_existing_minimal() {
        let width = 64;
        let height = 64;
        let count = 8;
        let num_baselines = 2;

        // parameters for simulated noise
        let noise_bl = 1;
        let noise_x = 32;
        let noise_y = 32;
        let noise_z = 1;
        let signal_val = 0 as f32;
        let noise_val = 0xffffff as f32;

        // parameters for simulated pre-existing flag
        let existing_bl = 0;
        let existing_x = 1;
        let existing_y = 2;

        let aoflagger = unsafe { cxx_aoflagger_new() };
        let mut baseline_imgsets: Vec<UniquePtr<CxxImageSet>> = (0..num_baselines)
            .into_iter()
            .map(|_| unsafe { aoflagger.MakeImageSet(width, height, count, signal_val, width) })
            .collect();

        // modify the imageset to add some synthetic noise
        let img_stride = baseline_imgsets[0].HorizontalStride();
        baseline_imgsets[noise_bl].pin_mut().ImageBufferMut(noise_z)
            [noise_y * img_stride + noise_x] = noise_val;

        let strategy_file_minimal = aoflagger.FindStrategyFileGeneric(&String::from("minimal"));

        let mut baseline_flagmasks: Vec<UniquePtr<CxxFlagMask>> = (0..num_baselines)
            .into_iter()
            .map(|_| unsafe { aoflagger.MakeFlagMask(width, height, false) })
            .collect();
        let flag_stride = baseline_flagmasks[0].HorizontalStride();
        baseline_flagmasks[existing_bl].pin_mut().BufferMut()
            [existing_y * flag_stride + existing_x] = true;

        flag_imgsets_existing(
            &aoflagger,
            &strategy_file_minimal,
            &baseline_imgsets,
            &mut baseline_flagmasks,
            true,
        );

        assert!(!baseline_flagmasks[0].Buffer()[0]);
        assert!(!baseline_flagmasks[0].Buffer()[noise_y * flag_stride + noise_x]);
        assert!(baseline_flagmasks[existing_bl].Buffer()[existing_y * flag_stride + existing_x]);
        assert!(!baseline_flagmasks[noise_bl].Buffer()[0]);
        assert!(baseline_flagmasks[noise_bl].Buffer()[noise_y * flag_stride + noise_x]);
    }

    #[test]
    fn test_flag_imgsets_existing_ord() {
        let aoflagger = unsafe { cxx_aoflagger_new() };
        let context = get_mwa_ord_context();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;

        let baseline_imgsets = context_to_baseline_imgsets(
            &aoflagger,
            &context,
            &img_coarse_chan_idxs,
            &img_timestep_idxs,
            None,
        );

        let strategy_filename = &aoflagger.FindStrategyFileMWA();
        let mut baseline_flagmasks = init_baseline_flagmasks(
            &aoflagger,
            &context,
            &img_coarse_chan_idxs,
            &img_timestep_idxs,
            Some(get_antenna_flags(&context)),
        );

        flag_imgsets_existing(
            &aoflagger,
            &strategy_filename,
            &baseline_imgsets,
            &mut baseline_flagmasks,
            true,
        );

        // test that flagged antennas are indeed flagged

        assert!(!baseline_flagmasks[0].Buffer()[7]);
        assert!(baseline_flagmasks[63].Buffer()[7]);
        assert!(baseline_flagmasks[62].Buffer()[7]);
        assert!(baseline_flagmasks[11].Buffer()[7]);
        assert!(baseline_flagmasks[111].Buffer()[7]);
        assert!(baseline_flagmasks[91].Buffer()[7]);
        assert!(baseline_flagmasks[95].Buffer()[7]);
        assert!(baseline_flagmasks[93].Buffer()[7]);
        assert!(baseline_flagmasks[80].Buffer()[7]);
        assert!(baseline_flagmasks[87].Buffer()[7]);
    }

    #[test]
    fn test_write_flags_mwax_minimal() {
        let context = get_mwax_context();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();

        let width = img_timestep_idxs.len();
        let height = context.num_common_coarse_chans
            * context.metafits_context.num_corr_fine_chans_per_coarse;

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
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| context.coarse_chans[chan].gpubox_number)
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
            filename_template.to_str().unwrap(),
            &selected_gpuboxes,
            context.mwa_version,
        )
        .unwrap();
        let chan_header_flags_raw = flag_file_set.read_chan_header_flags_raw().unwrap();
        let (chan1_header, chan1_flags_raw) =
            chan_header_flags_raw.get(&selected_gpuboxes[0]).unwrap();
        assert_eq!(chan1_header.gpubox_id, gpubox_ids[0]);
        let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;

        let num_baselines = chan1_header.num_ants * (chan1_header.num_ants + 1) / 2;
        assert_eq!(chan1_header.num_timesteps, context.num_timesteps);
        assert_eq!(num_baselines, context.metafits_context.num_baselines);
        assert_eq!(chan1_header.num_channels, fine_chans_per_coarse);
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
            let offset = row_idx * fine_chans_per_coarse + fine_chan_idx;
            assert_eq!(
                &chan1_flags_raw[offset], expected_flag,
                "with timestep {}, baseline {}, fine_chan {}, expected {} at row_idx {}, offset {}",
                timestep_idx, baseline_idx, fine_chan_idx, expected_flag, row_idx, offset
            );
        }
    }

    /// A vector of all fine channel frequencies for the provided coarse channels
    /// in order of provided `coarse_chan_idxs`, and ascending fine channel index
    fn _get_all_freqs_hz(context: &CorrelatorContext, coarse_chan_idxs: &[usize]) -> Vec<u32> {
        let coarse_chans = context.coarse_chans.clone();
        let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
        let fine_chan_width_hz = context.metafits_context.corr_fine_chan_width_hz;
        coarse_chan_idxs
            .iter()
            .flat_map(|&coarse_chan_idx| {
                let coarse_chan = &coarse_chans[coarse_chan_idx];
                let chan_start_hz = coarse_chan.chan_start_hz;
                (0..fine_chans_per_coarse).map(move |fine_chan_idx| {
                    chan_start_hz + (fine_chan_idx as u32 * fine_chan_width_hz)
                })
            })
            .collect()
    }

    #[test]
    fn test_cable_length_corrections_mwax() {
        let context = get_mwax_context();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(img_timestep_idxs.len(), 4);
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;

        let all_freqs_hz = _get_all_freqs_hz(&context, img_coarse_chan_idxs);

        let aoflagger = unsafe { cxx_aoflagger_new() };

        let width = img_timestep_idxs.len();
        let stride = (((width - 1) / 8) + 1) * 8;

        let mut baseline_imgsets = context_to_baseline_imgsets(
            &aoflagger,
            &context,
            &img_coarse_chan_idxs,
            &img_timestep_idxs,
            None,
        );

        let viz_0_xx_0_0_re = baseline_imgsets[0].ImageBuffer(0)[0];
        let viz_0_xx_0_0_im = baseline_imgsets[0].ImageBuffer(1)[0];
        let viz_1_xx_0_0_re = baseline_imgsets[1].ImageBuffer(0)[0];
        let viz_1_xx_0_0_im = baseline_imgsets[1].ImageBuffer(1)[0];
        let viz_1_xy_0_0_re = baseline_imgsets[1].ImageBuffer(2)[0];
        let viz_1_xy_0_0_im = baseline_imgsets[1].ImageBuffer(3)[0];
        let viz_1_yx_0_0_re = baseline_imgsets[1].ImageBuffer(4)[0];
        let viz_1_yx_0_0_im = baseline_imgsets[1].ImageBuffer(5)[0];
        let viz_1_yy_0_0_re = baseline_imgsets[1].ImageBuffer(6)[0];
        let viz_1_yy_0_0_im = baseline_imgsets[1].ImageBuffer(7)[0];
        let viz_1_yy_3_3_re = baseline_imgsets[1].ImageBuffer(6)[3 * stride + 3];
        let viz_1_yy_3_3_im = baseline_imgsets[1].ImageBuffer(7)[3 * stride + 3];

        assert_eq!(viz_0_xx_0_0_re as f32, 0x410000 as f32);
        assert_eq!(viz_0_xx_0_0_im as f32, 0x410001 as f32);
        assert_eq!(viz_1_xx_0_0_re as f32, 0x410010 as f32);
        assert_eq!(viz_1_xx_0_0_im as f32, 0x410011 as f32);
        assert_eq!(viz_1_xy_0_0_re as f32, 0x410012 as f32);
        assert_eq!(viz_1_xy_0_0_im as f32, 0x410013 as f32);
        assert_eq!(viz_1_yx_0_0_re as f32, 0x410014 as f32);
        assert_eq!(viz_1_yx_0_0_im as f32, 0x410015 as f32);
        assert_eq!(viz_1_yy_0_0_re as f32, 0x410016 as f32);
        assert_eq!(viz_1_yy_0_0_im as f32, 0x410017 as f32);
        assert_eq!(viz_1_yy_3_3_re as f32, 0x41071e as f32);
        assert_eq!(viz_1_yy_3_3_im as f32, 0x41071f as f32);

        // baseline 1, input 1, pol x
        let length_1_1_x = &context.metafits_context.antennas
            [context.metafits_context.baselines[1].ant1_index]
            .rfinput_x
            .electrical_length_m;
        // baseline 1, input 1, pol y
        let length_1_1_y = &context.metafits_context.antennas
            [context.metafits_context.baselines[1].ant1_index]
            .rfinput_y
            .electrical_length_m;
        // baseline 1, input 2, pol x
        let length_1_2_x = &context.metafits_context.antennas
            [context.metafits_context.baselines[1].ant2_index]
            .rfinput_x
            .electrical_length_m;
        // baseline 1, input 2, pol y
        let length_1_2_y = &context.metafits_context.antennas
            [context.metafits_context.baselines[1].ant2_index]
            .rfinput_y
            .electrical_length_m;

        // baseline 1, pol XX, cc 0, fc 0
        let length_m_1_xx = length_1_2_x - length_1_1_x;
        let angle_1_xx_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_1_xx * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM;
        let (sin_1_xx_0_f64, cos_1_xx_0_f64) = angle_1_xx_0.sin_cos();
        let (sin_1_xx_0, cos_1_xx_0) = (sin_1_xx_0_f64 as f32, cos_1_xx_0_f64 as f32);
        // baseline 1, pol XY, cc 0, fc 0
        let length_m_1_xy = length_1_2_y - length_1_1_x;
        let angle_1_xy_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_1_xy * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM;
        let (sin_1_xy_0_f64, cos_1_xy_0_f64) = angle_1_xy_0.sin_cos();
        let (sin_1_xy_0, cos_1_xy_0) = (sin_1_xy_0_f64 as f32, cos_1_xy_0_f64 as f32);
        // baseline 1, pol YX, cc 0, fc 0
        let length_m_1_yx = length_1_2_x - length_1_1_y;
        let angle_1_yx_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_1_yx * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM;
        let (sin_1_yx_0_f64, cos_1_yx_0_f64) = angle_1_yx_0.sin_cos();
        let (sin_1_yx_0, cos_1_yx_0) = (sin_1_yx_0_f64 as f32, cos_1_yx_0_f64 as f32);
        // baseline 1, pol YY, cc 0, fc 0
        let length_m_1_yy = length_1_2_y - length_1_1_y;
        let angle_1_yy_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_1_yy * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM;
        let (sin_1_yy_0_f64, cos_1_yy_0_f64) = angle_1_yy_0.sin_cos();
        let (sin_1_yy_0, cos_1_yy_0) = (sin_1_yy_0_f64 as f32, cos_1_yy_0_f64 as f32);
        // baseline 1, pol YY, cc 1, fc 1
        let angle_1_yy_3: f64 =
            -2.0 * std::f64::consts::PI * length_m_1_yy * (all_freqs_hz[3] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM;
        let (sin_1_yy_3_f64, cos_1_yy_3_f64) = angle_1_yy_3.sin_cos();
        let (sin_1_yy_3, cos_1_yy_3) = (sin_1_yy_3_f64 as f32, cos_1_yy_3_f64 as f32);

        correct_cable_lengths(&context, &mut baseline_imgsets);

        // there should be no difference in baseline 0
        test_imgset_val!(baseline_imgsets[0], 0, 0, 0, viz_0_xx_0_0_re);
        test_imgset_val!(baseline_imgsets[0], 1, 0, 0, viz_0_xx_0_0_im);

        ////
        // baseline 1 should be rotated
        ////
        // baseline 1, pol xx, cc 0, fc 0, ts 0
        let rot_1_xx_0_0_re = (cos_1_xx_0 * viz_1_xx_0_0_re - sin_1_xx_0 * viz_1_xx_0_0_im) as f32;
        let rot_1_xx_0_0_im = (sin_1_xx_0 * viz_1_xx_0_0_re + cos_1_xx_0 * viz_1_xx_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[1], 0, 0, 0, rot_1_xx_0_0_re);
        test_imgset_val!(baseline_imgsets[1], 1, 0, 0, rot_1_xx_0_0_im);
        // baseline 1, pol xy, cc 0, fc 0, ts 0
        let rot_1_xy_0_0_re = (cos_1_xy_0 * viz_1_xy_0_0_re - sin_1_xy_0 * viz_1_xy_0_0_im) as f32;
        let rot_1_xy_0_0_im = (sin_1_xy_0 * viz_1_xy_0_0_re + cos_1_xy_0 * viz_1_xy_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[1], 2, 0, 0, rot_1_xy_0_0_re);
        test_imgset_val!(baseline_imgsets[1], 3, 0, 0, rot_1_xy_0_0_im);
        // baseline 1, pol yx, cc 0, fc 0, ts 0
        let rot_1_yx_0_0_re = (cos_1_yx_0 * viz_1_yx_0_0_re - sin_1_yx_0 * viz_1_yx_0_0_im) as f32;
        let rot_1_yx_0_0_im = (sin_1_yx_0 * viz_1_yx_0_0_re + cos_1_yx_0 * viz_1_yx_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[1], 4, 0, 0, rot_1_yx_0_0_re);
        test_imgset_val!(baseline_imgsets[1], 5, 0, 0, rot_1_yx_0_0_im);
        // baseline 1, pol yy, cc 0, fc 0, ts 0
        let rot_1_yy_0_0_re = (cos_1_yy_0 * viz_1_yy_0_0_re - sin_1_yy_0 * viz_1_yy_0_0_im) as f32;
        let rot_1_yy_0_0_im = (sin_1_yy_0 * viz_1_yy_0_0_re + cos_1_yy_0 * viz_1_yy_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[1], 6, 0, 0, rot_1_yy_0_0_re);
        test_imgset_val!(baseline_imgsets[1], 7, 0, 0, rot_1_yy_0_0_im);
        // baseline 1, pol yy, cc 1, fc 1, ts 1
        let rot_1_yy_3_3_re = (cos_1_yy_3 * viz_1_yy_3_3_re - sin_1_yy_3 * viz_1_yy_3_3_im) as f32;
        let rot_1_yy_3_3_im = (sin_1_yy_3 * viz_1_yy_3_3_re + cos_1_yy_3 * viz_1_yy_3_3_im) as f32;
        test_imgset_val!(baseline_imgsets[1], 6, 3, 3, rot_1_yy_3_3_re);
        test_imgset_val!(baseline_imgsets[1], 7, 3, 3, rot_1_yy_3_3_im);
    }

    #[test]
    fn test_cable_length_corrections_ord() {
        let context = get_mwa_ord_context();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(img_timestep_idxs.len(), 4);
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;

        let all_freqs_hz = _get_all_freqs_hz(&context, img_coarse_chan_idxs);

        let aoflagger = unsafe { cxx_aoflagger_new() };

        let width = img_timestep_idxs.len();
        let stride = (((width - 1) / 8) + 1) * 8;

        let mut baseline_imgsets = context_to_baseline_imgsets(
            &aoflagger,
            &context,
            &img_coarse_chan_idxs,
            &img_timestep_idxs,
            None,
        );

        let viz_0_xx_0_0_re = baseline_imgsets[0].ImageBuffer(0)[0];
        let viz_0_xx_0_0_im = baseline_imgsets[0].ImageBuffer(1)[0];
        let viz_5_xx_0_0_re = baseline_imgsets[5].ImageBuffer(0)[0];
        let viz_5_xx_0_0_im = baseline_imgsets[5].ImageBuffer(1)[0];
        let viz_5_xy_0_0_re = baseline_imgsets[5].ImageBuffer(2)[0];
        let viz_5_xy_0_0_im = baseline_imgsets[5].ImageBuffer(3)[0];
        let viz_5_yx_0_0_re = baseline_imgsets[5].ImageBuffer(4)[0];
        let viz_5_yx_0_0_im = baseline_imgsets[5].ImageBuffer(5)[0];
        let viz_5_yy_0_0_re = baseline_imgsets[5].ImageBuffer(6)[0];
        let viz_5_yy_0_0_im = baseline_imgsets[5].ImageBuffer(7)[0];
        let viz_5_yy_3_3_re = baseline_imgsets[5].ImageBuffer(6)[3 * stride + 3];
        let viz_5_yy_3_3_im = baseline_imgsets[5].ImageBuffer(7)[3 * stride + 3];

        assert_eq!(viz_0_xx_0_0_re as f32, 0x10c5be as f32);
        assert_eq!(viz_0_xx_0_0_im as f32, 0x10c5bf as f32);
        assert_eq!(viz_5_xx_0_0_re as f32, 0x10f1ce as f32);
        assert_eq!(viz_5_xx_0_0_im as f32, -0x10f1cf as f32);
        assert_eq!(viz_5_xy_0_0_re as f32, 0x10ea26 as f32);
        assert_eq!(viz_5_xy_0_0_im as f32, -0x10ea27 as f32);
        assert_eq!(viz_5_yx_0_0_re as f32, 0x10f1be as f32);
        assert_eq!(viz_5_yx_0_0_im as f32, -0x10f1bf as f32);
        assert_eq!(viz_5_yy_0_0_re as f32, 0x10ea16 as f32);
        assert_eq!(viz_5_yy_0_0_im as f32, -0x10ea17 as f32);
        assert_eq!(viz_5_yy_3_3_re as f32, 0x0dec16 as f32);
        assert_eq!(viz_5_yy_3_3_im as f32, -0x0dec17 as f32);

        // baseline 5, input 1, pol x
        let length_5_1_x = &context.metafits_context.antennas
            [context.metafits_context.baselines[5].ant1_index]
            .rfinput_x
            .electrical_length_m;
        // baseline 5, input 1, pol y
        let length_5_1_y = &context.metafits_context.antennas
            [context.metafits_context.baselines[5].ant1_index]
            .rfinput_y
            .electrical_length_m;
        // baseline 5, input 2, pol x
        let length_5_2_x = &context.metafits_context.antennas
            [context.metafits_context.baselines[5].ant2_index]
            .rfinput_x
            .electrical_length_m;
        // baseline 5, input 2, pol y
        let length_5_2_y = &context.metafits_context.antennas
            [context.metafits_context.baselines[5].ant2_index]
            .rfinput_y
            .electrical_length_m;

        // baseline 1, pol XX, cc 0, fc 0
        let length_m_5_xx = length_5_2_x - length_5_1_x;
        let angle_5_xx_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_5_xx * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM;
        let (sin_5_xx_0_f64, cos_5_xx_0_f64) = angle_5_xx_0.sin_cos();
        let (sin_5_xx_0, cos_5_xx_0) = (sin_5_xx_0_f64 as f32, cos_5_xx_0_f64 as f32);
        // baseline 1, pol XY, cc 0, fc 0
        let length_m_5_xy = length_5_2_y - length_5_1_x;
        let angle_5_xy_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_5_xy * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM;
        let (sin_5_xy_0_f64, cos_5_xy_0_f64) = angle_5_xy_0.sin_cos();
        let (sin_5_xy_0, cos_5_xy_0) = (sin_5_xy_0_f64 as f32, cos_5_xy_0_f64 as f32);
        // baseline 1, pol YX, cc 0, fc 0
        let length_m_5_yx = length_5_2_x - length_5_1_y;
        let angle_5_yx_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_5_yx * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM;
        let (sin_5_yx_0_f64, cos_5_yx_0_f64) = angle_5_yx_0.sin_cos();
        let (sin_5_yx_0, cos_5_yx_0) = (sin_5_yx_0_f64 as f32, cos_5_yx_0_f64 as f32);
        // baseline 1, pol YY, cc 0, fc 0
        let length_m_5_yy = length_5_2_y - length_5_1_y;
        let angle_5_yy_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_5_yy * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM;
        let (sin_5_yy_0_f64, cos_5_yy_0_f64) = angle_5_yy_0.sin_cos();
        let (sin_5_yy_0, cos_5_yy_0) = (sin_5_yy_0_f64 as f32, cos_5_yy_0_f64 as f32);
        // baseline 1, pol YY, cc 1, fc 1
        let angle_5_yy_3: f64 =
            -2.0 * std::f64::consts::PI * length_m_5_yy * (all_freqs_hz[3] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM;
        let (sin_5_yy_3_f64, cos_5_yy_3_f64) = angle_5_yy_3.sin_cos();
        let (sin_5_yy_3, cos_5_yy_3) = (sin_5_yy_3_f64 as f32, cos_5_yy_3_f64 as f32);

        correct_cable_lengths(&context, &mut baseline_imgsets);

        // there should be no difference in baseline 0
        test_imgset_val!(baseline_imgsets[0], 0, 0, 0, viz_0_xx_0_0_re);
        test_imgset_val!(baseline_imgsets[0], 1, 0, 0, viz_0_xx_0_0_im);

        ////
        // baseline 1 should be rotated
        ////
        // baseline 1, pol xx, cc 0, fc 0, ts 0
        let rot_5_xx_0_0_re = (cos_5_xx_0 * viz_5_xx_0_0_re - sin_5_xx_0 * viz_5_xx_0_0_im) as f32;
        let rot_5_xx_0_0_im = (sin_5_xx_0 * viz_5_xx_0_0_re + cos_5_xx_0 * viz_5_xx_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[5], 0, 0, 0, rot_5_xx_0_0_re);
        test_imgset_val!(baseline_imgsets[5], 1, 0, 0, rot_5_xx_0_0_im);
        // baseline 1, pol xy, cc 0, fc 0, ts 0
        let rot_5_xy_0_0_re = (cos_5_xy_0 * viz_5_xy_0_0_re - sin_5_xy_0 * viz_5_xy_0_0_im) as f32;
        let rot_5_xy_0_0_im = (sin_5_xy_0 * viz_5_xy_0_0_re + cos_5_xy_0 * viz_5_xy_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[5], 2, 0, 0, rot_5_xy_0_0_re);
        test_imgset_val!(baseline_imgsets[5], 3, 0, 0, rot_5_xy_0_0_im);
        // baseline 1, pol yx, cc 0, fc 0, ts 0
        let rot_5_yx_0_0_re = (cos_5_yx_0 * viz_5_yx_0_0_re - sin_5_yx_0 * viz_5_yx_0_0_im) as f32;
        let rot_5_yx_0_0_im = (sin_5_yx_0 * viz_5_yx_0_0_re + cos_5_yx_0 * viz_5_yx_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[5], 4, 0, 0, rot_5_yx_0_0_re);
        test_imgset_val!(baseline_imgsets[5], 5, 0, 0, rot_5_yx_0_0_im);
        // baseline 1, pol yy, cc 0, fc 0, ts 0
        let rot_5_yy_0_0_re = (cos_5_yy_0 * viz_5_yy_0_0_re - sin_5_yy_0 * viz_5_yy_0_0_im) as f32;
        let rot_5_yy_0_0_im = (sin_5_yy_0 * viz_5_yy_0_0_re + cos_5_yy_0 * viz_5_yy_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[5], 6, 0, 0, rot_5_yy_0_0_re);
        test_imgset_val!(baseline_imgsets[5], 7, 0, 0, rot_5_yy_0_0_im);
        // baseline 1, pol yy, cc 1, fc 1, ts 1
        let rot_5_yy_3_3_re = (cos_5_yy_3 * viz_5_yy_3_3_re - sin_5_yy_3 * viz_5_yy_3_3_im) as f32;
        let rot_5_yy_3_3_im = (sin_5_yy_3 * viz_5_yy_3_3_re + cos_5_yy_3 * viz_5_yy_3_3_im) as f32;
        test_imgset_val!(baseline_imgsets[5], 6, 3, 3, rot_5_yy_3_3_re);
        test_imgset_val!(baseline_imgsets[5], 7, 3, 3, rot_5_yy_3_3_im);
    }

    fn _get_this_flagmask(aoflagger: &CxxAOFlagger) -> UniquePtr<CxxFlagMask> {
        let mut this_flagmask = unsafe { aoflagger.MakeFlagMask(2, 2, false) };
        let this_buffer = this_flagmask.pin_mut().BufferMut();
        this_buffer[0] = false;
        this_buffer[1] = false;
        this_buffer[2] = true;
        this_buffer[3] = true;
        this_flagmask
    }

    fn _get_other_flagmask(aoflagger: &CxxAOFlagger) -> UniquePtr<CxxFlagMask> {
        let mut other_flagmask = unsafe { aoflagger.MakeFlagMask(2, 2, false) };
        let other_buffer = other_flagmask.pin_mut().BufferMut();
        other_buffer[0] = false;
        other_buffer[1] = true;
        other_buffer[2] = false;
        other_buffer[3] = true;
        other_flagmask
    }

    #[test]
    fn test_flagmask_or() {
        let aoflagger = unsafe { cxx_aoflagger_new() };
        let mut this_flagmask = _get_this_flagmask(&aoflagger);
        let other_flagmask = _get_other_flagmask(&aoflagger);

        flagmask_or(&mut this_flagmask, &other_flagmask);

        let this_buffer = this_flagmask.Buffer();
        assert!(!this_buffer[0]);
        assert!(this_buffer[1]);
        assert!(this_buffer[2]);
        assert!(this_buffer[3]);
    }

    #[test]
    fn test_flagmask_set() {
        let aoflagger = unsafe { cxx_aoflagger_new() };
        let mut this_flagmask = _get_this_flagmask(&aoflagger);
        let other_flagmask = _get_other_flagmask(&aoflagger);

        flagmask_set(&mut this_flagmask, &other_flagmask);

        let this_buffer = this_flagmask.Buffer();
        assert!(!this_buffer[0]);
        assert!(this_buffer[1]);
        assert!(!this_buffer[2]);
        assert!(this_buffer[3]);
    }
}
