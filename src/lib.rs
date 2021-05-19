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
use error::BirliError;
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

use log::trace;

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

/// Read aan observation's visibilities into a vector containing a [`CxxImageSet`]s for each
/// baseline in the observation, given a [`CxxAOFlagger`] instance and that observation's
/// [`mwalib::CorrelatorContext`].
///
/// [`mwalib::CorrelatorContext`]: https://docs.rs/mwalib/0.7.0/mwalib/struct.CorrelatorContext.html
///
/// # Examples
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
    trace!("start context_to_baseline_imgsets");

    // TODO: although collect::Vec<_> is more readable, it uses the stack less
    // efficiently than Vec
    // TODO: error handling instead of unwrap.

    let num_coarse_chans = context.num_coarse_chans;
    let num_timesteps = context.num_timesteps;
    let num_baselines = context.metafits_context.num_baselines;
    let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
    let floats_per_finechan = context.metafits_context.num_visibility_pols * 2;
    let floats_per_baseline = fine_chans_per_coarse * floats_per_finechan;
    let height = context.num_coarse_chans * fine_chans_per_coarse;
    let width = context.num_timesteps;
    let img_stride = (((width - 1) / 8) + 1) * 8;

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
    let mut baseline_imgsets: Vec<UniquePtr<CxxImageSet>> = (0..num_baselines)
        .into_par_iter()
        .map(|_| {
            let imgset = unsafe { aoflagger.MakeImageSet(width, height, 8, 0 as f32, width) };
            allocation_progress.inc(1);
            imgset
        })
        .collect();

    let baseline_imgsets_arc = Arc::new(&mut baseline_imgsets);
    allocation_progress.finish();

    let num_producers = num_coarse_chans;
    let num_pols_complex = context.metafits_context.num_visibility_pols * 2;

    // A queue of coarse channel indices for the producers to work through
    let (tx_coarse_chan_idx, rx_coarse_chan_idx) = unbounded();
    // a queue of raw gpufits visibility image buffers for each complex polarization.
    let pol_img_queues: Vec<(Sender<_>, Receiver<_>)> = (0..num_pols_complex)
        .map(|_| bounded(num_producers))
        .collect();

    // a progress bar containing the progress bars associated with loading the
    // observation's HDUs
    let multi_progress = MultiProgress::new();
    // a vector of progress bars for the visibility reading progress of each
    // channel.
    let read_progress: Vec<ProgressBar> = (0..num_coarse_chans)
        .map(|coarse_chan_idx| {
            multi_progress.add(
                ProgressBar::new(num_timesteps as _)
                    .with_style(
                        ProgressStyle::default_bar()
                            .template("{msg:16}: [{wide_bar}] {pos:4}/{len:4}")
                            .progress_chars("=> "),
                    )
                    .with_position(0)
                    .with_message(format!("coarse chan {:3}", coarse_chan_idx)),
            )
        })
        .collect();
    // a vector of progress bars for the visibility writing progress of each
    // complex polarization,
    let write_progress: Vec<ProgressBar> = (0..num_pols_complex)
        .map(|pol_idx| {
            multi_progress.add(
                ProgressBar::new((num_timesteps * num_coarse_chans) as _)
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
        ProgressBar::new((num_timesteps * num_coarse_chans * num_pols_complex) as _)
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
    let total_progess_arc = Arc::new(total_progress);

    thread::scope(|scope| {
        // Spawn a thread to deal with the progress bar
        scope.spawn(|_| {
            multi_progress.join().unwrap();
        });

        // Queue up coarse channels to do
        (0..num_coarse_chans).for_each(|coarse_chan_idx| {
            tx_coarse_chan_idx.send(coarse_chan_idx).unwrap();
        });
        // This indicates the the producer threads that there are no more coarse
        // channels to process.
        drop(tx_coarse_chan_idx);

        // Create multiple producer threads
        for _ in 0..num_producers {
            let rx_coarse_chan_idx_worker = rx_coarse_chan_idx.clone();
            let pol_img_queues_worker = pol_img_queues.to_owned();
            let read_progress_worker = read_progress_arc.to_owned();
            scope.spawn(move |_| {
                // Each producer thread consumes the HDUs from one coarse
                // channel at a time.
                for coarse_chan_idx in rx_coarse_chan_idx_worker.iter() {
                    (0..num_timesteps).for_each(|timestep_idx| {
                        let img_buf = Arc::new(
                            context
                                .read_by_baseline(timestep_idx, coarse_chan_idx)
                                .unwrap(),
                        );
                        // Producer sends the visibility buffer to a separate
                        // queue for each complex polarization.
                        pol_img_queues_worker.iter().for_each(|(tx_img, _)| {
                            tx_img
                                .send((coarse_chan_idx, timestep_idx, img_buf.clone()))
                                .unwrap();
                        });
                        read_progress_worker[coarse_chan_idx].inc(1);
                    });
                    read_progress_worker[coarse_chan_idx].finish_and_clear();
                }
                pol_img_queues_worker.into_iter().for_each(|(tx_img, _)| {
                    drop(tx_img);
                });
            });
        }

        // create a consumer thread for each complex polarization
        for (pol_idx, (tx_img, rx_img)) in pol_img_queues.into_iter().enumerate() {
            // This ensures that
            drop(tx_img);
            let rx_img_worker = rx_img.to_owned();
            let pol_idx_worker = pol_idx.to_owned();
            let baseline_imgsets_worker = baseline_imgsets_arc.to_owned();
            let write_progress_worker = write_progress_arc.to_owned();
            let total_progess_worker = total_progess_arc.to_owned();
            scope.spawn(move |_| {
                // The thread consumes an image from it's img queue
                for (coarse_chan_idx, timestep_idx, img_buf) in rx_img_worker.iter() {
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
                            let imgset_buf = unsafe { imgset.ImageBufferMutUnsafe(pol_idx_worker) };
                            baseline_chunk
                                .chunks_exact(floats_per_finechan)
                                .enumerate()
                                .for_each(|(fine_chan_idx, fine_chan_chunk)| {
                                    let img_x = timestep_idx;
                                    let img_y =
                                        fine_chans_per_coarse * coarse_chan_idx + fine_chan_idx;
                                    unsafe {
                                        *imgset_buf.get_unchecked_mut(img_y * img_stride + img_x) =
                                            *fine_chan_chunk.get_unchecked(pol_idx_worker)
                                    };
                                });
                        });
                    write_progress_worker[pol_idx].inc(1);
                    total_progess_worker.inc(1);
                }
                write_progress_worker[pol_idx].finish_and_clear();
                if total_progess_worker.position() >= total_progess_worker.length() {
                    total_progess_worker.finish();
                }
            });
        }

        // consume the rx_img queue
        // total_progress.finish();
    })
    .unwrap();

    trace!("end context_to_baseline_imgsets");

    baseline_imgsets
}

const SPEED_OF_LIGHT_IN_VACUUM: f64 = 299792458.0; // speed of light in m/s

// void Cotter::correctCableLength(ImageSet& imageSet, size_t polarization, double cableDelay) const
// {
// 	float *reals = imageSet.ImageBuffer(polarization*2);
// 	float *imags = imageSet.ImageBuffer(polarization*2+1);

// 	for(size_t y=0; y!=imageSet.Height(); ++y)
// 	{
// 		double angle = -2.0 * M_PI * cableDelay * _channelFrequenciesHz[y] / SPEED_OF_LIGHT;
// 		double rotSinl, rotCosl;
// 		sincos(angle, &rotSinl, &rotCosl);
// 		float rotSin = rotSinl, rotCos = rotCosl;

// 		/// @todo This should use actual time step count in window
// 		float *realPtr = reals + y * imageSet.HorizontalStride();
// 		float *imagPtr = imags + y * imageSet.HorizontalStride();
// 		for(size_t x=0; x!=imageSet.Width(); ++x)
// 		{
// 			float r = *realPtr;
// 			*realPtr = rotCos * r - rotSin * (*imagPtr);
// 			*imagPtr = rotSin * r + rotCos * (*imagPtr);
// 			++realPtr;
// 			++imagPtr;
// 		}
// 	}
// }

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

// fn _correct_cable_length_buffers_precise(
//     freq_hz: &u32,
//     electrical_length_m: &f64,
//     buf_re: &mut [f32],
//     buf_im: &mut [f32],
// ) {
//     let angle: f64 = -2.0 * PI * electrical_length_m * (*freq_hz as f64) / SPEED_OF_LIGHT_IN_VACUUM;
//     let (sin_angle, cos_angle) = angle.sin_cos();

//     izip!(buf_re.iter_mut(), buf_im.iter_mut()).for_each(|(re, im)| {
//         let vis_re = *re as f64;
//         let vis_im = *im as f64;
//         *re = (cos_angle * vis_re - sin_angle * vis_im) as f32;
//         *im = (sin_angle * vis_re + cos_angle * vis_im) as f32;
//     })
// }

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
/// let baseline_imgsets = context_to_baseline_imgsets(&aoflagger, &context);
///
/// ccorrect_cable_lengths(&context, baseline_imgsets);
/// ```
pub fn correct_cable_lengths(
    context: &CorrelatorContext,
    mut baseline_imgsets: Vec<UniquePtr<CxxImageSet>>,
) -> Result<Vec<UniquePtr<CxxImageSet>>, BirliError> {
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
            let chan_start_hz = coarse_chan.chan_start_hz.clone();
            (0..fine_chans_per_coarse)
                .map(move |fine_chan_idx| {
                    chan_start_hz + (fine_chan_idx as u32 * fine_chan_width_hz)
                })
                .clone()
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
    Ok(baseline_imgsets)
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
    trace!("start write_flags");

    // TODO: error handling instead of unwrap.

    let mut flag_file_set = FlagFileSet::new(context, filename_template, &gpubox_ids).unwrap();
    flag_file_set
        .write_baseline_flagmasks(&context, baseline_flagmasks)
        .unwrap();

    trace!("end write_flags");
}

#[cfg(test)]
mod tests {
    // TODO: Why does clippy think CxxImageSet.ImageBuffer() is &[f64]?
    #![allow(clippy::float_cmp)]

    use super::{context_to_baseline_imgsets, flag_imgsets, write_flags, FlagFileSet};
    use crate::{
        correct_cable_lengths,
        cxx_aoflagger::ffi::{cxx_aoflagger_new, CxxFlagMask, CxxImageSet},
        SPEED_OF_LIGHT_IN_VACUUM,
    };
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
                $imgset.ImageBuffer($imgset_idx)[$x * $img_stride + $y].round(),
                ($val as f32).round()
            )
        };
    }

    #[test]
    fn test_context_to_baseline_imgsets_mwax() {
        let context = get_mwax_context();
        let width = context.num_timesteps;
        let stride = (((width - 1) / 8) + 1) * 8;

        let aoflagger = unsafe { cxx_aoflagger_new() };
        let baseline_imgsets = context_to_baseline_imgsets(&aoflagger, &context);

        test_imgset_val!(baseline_imgsets[0], 0, stride, 0, 0, 0x410000);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 0, 1, 0x410100);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 0, 2, 0x410200);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 0, 3, 0x410300);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 0, 4, 0);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 0, 5, 0);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 0, 6, 0);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 0, 7, 0);

        test_imgset_val!(baseline_imgsets[0], 0, stride, 1, 0, 0x410008);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 1, 1, 0x410108);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 1, 2, 0x410208);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 1, 3, 0x410308);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 1, 4, 0);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 1, 5, 0);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 1, 6, 0);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 1, 7, 0);

        test_imgset_val!(baseline_imgsets[0], 0, stride, 2, 0, 0x410400);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 2, 1, 0x410500);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 2, 2, 0x410600);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 2, 3, 0x410700);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 2, 4, 0);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 2, 5, 0);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 2, 6, 0);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 2, 7, 0);

        test_imgset_val!(baseline_imgsets[0], 0, stride, 3, 0, 0x410408);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 3, 1, 0x410508);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 3, 2, 0x410608);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 3, 3, 0x410708);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 3, 4, 0);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 3, 5, 0);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 3, 6, 0);
        test_imgset_val!(baseline_imgsets[0], 0, stride, 3, 7, 0);

        test_imgset_val!(baseline_imgsets[0], 1, stride, 0, 0, 0x410001);
        test_imgset_val!(baseline_imgsets[0], 1, stride, 0, 1, 0x410101);
        test_imgset_val!(baseline_imgsets[0], 1, stride, 0, 2, 0x410201);
        test_imgset_val!(baseline_imgsets[0], 1, stride, 0, 3, 0x410301);
        test_imgset_val!(baseline_imgsets[0], 1, stride, 0, 4, 0);
        test_imgset_val!(baseline_imgsets[0], 1, stride, 0, 5, 0);
        test_imgset_val!(baseline_imgsets[0], 1, stride, 0, 6, 0);
        test_imgset_val!(baseline_imgsets[0], 1, stride, 0, 7, 0);

        /* ... */
        test_imgset_val!(baseline_imgsets[0], 7, stride, 0, 0, 0x410007);
        test_imgset_val!(baseline_imgsets[0], 7, stride, 0, 1, 0x410107);
        test_imgset_val!(baseline_imgsets[0], 7, stride, 0, 2, 0x410207);
        test_imgset_val!(baseline_imgsets[0], 7, stride, 0, 3, 0x410307);
        test_imgset_val!(baseline_imgsets[0], 7, stride, 0, 4, 0);
        test_imgset_val!(baseline_imgsets[0], 7, stride, 0, 5, 0);
        test_imgset_val!(baseline_imgsets[0], 7, stride, 0, 6, 0);
        test_imgset_val!(baseline_imgsets[0], 7, stride, 0, 7, 0);

        /* ... */

        test_imgset_val!(baseline_imgsets[2], 0, stride, 0, 0, 0x410020);
        test_imgset_val!(baseline_imgsets[2], 0, stride, 0, 1, 0x410120);
        test_imgset_val!(baseline_imgsets[2], 0, stride, 0, 2, 0x410220);
        test_imgset_val!(baseline_imgsets[2], 0, stride, 0, 3, 0x410320);
        test_imgset_val!(baseline_imgsets[2], 0, stride, 0, 4, 0);
        test_imgset_val!(baseline_imgsets[2], 0, stride, 0, 5, 0);
        test_imgset_val!(baseline_imgsets[2], 0, stride, 0, 6, 0);
        test_imgset_val!(baseline_imgsets[2], 0, stride, 0, 7, 0);

        test_imgset_val!(baseline_imgsets[2], 0, stride, 1, 0, 0x410028);
        test_imgset_val!(baseline_imgsets[2], 0, stride, 1, 1, 0x410128);
        test_imgset_val!(baseline_imgsets[2], 0, stride, 1, 2, 0x410228);
        test_imgset_val!(baseline_imgsets[2], 0, stride, 1, 3, 0x410328);
        test_imgset_val!(baseline_imgsets[2], 0, stride, 1, 4, 0);
        test_imgset_val!(baseline_imgsets[2], 0, stride, 1, 5, 0);
        test_imgset_val!(baseline_imgsets[2], 0, stride, 1, 6, 0);
        test_imgset_val!(baseline_imgsets[2], 0, stride, 1, 7, 0);
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

        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 0, 0x10c5be);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 1, 0x14c5be);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 2, 0x18c5be);
        test_imgset_val!(baseline_imgsets[0], 0, img_stride, 0, 3, 0x1cc5be);

        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 0, 0x10c5bf);
        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 1, 0x14c5bf);
        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 2, 0x18c5bf);
        test_imgset_val!(baseline_imgsets[0], 1, img_stride, 0, 3, 0x1cc5bf);

        test_imgset_val!(baseline_imgsets[0], 2, img_stride, 0, 0, 0x10c5ae);
        test_imgset_val!(baseline_imgsets[0], 2, img_stride, 0, 1, 0x14c5ae);
        test_imgset_val!(baseline_imgsets[0], 2, img_stride, 0, 2, 0x18c5ae);
        test_imgset_val!(baseline_imgsets[0], 2, img_stride, 0, 3, 0x1cc5ae);

        test_imgset_val!(baseline_imgsets[0], 3, img_stride, 0, 0, -0x10c5af);
        test_imgset_val!(baseline_imgsets[0], 3, img_stride, 0, 1, -0x14c5af);
        test_imgset_val!(baseline_imgsets[0], 3, img_stride, 0, 2, -0x18c5af);
        test_imgset_val!(baseline_imgsets[0], 3, img_stride, 0, 3, -0x1cc5af);

        test_imgset_val!(baseline_imgsets[0], 4, img_stride, 0, 0, 0x10c5ae);
        test_imgset_val!(baseline_imgsets[0], 4, img_stride, 0, 1, 0x14c5ae);
        test_imgset_val!(baseline_imgsets[0], 4, img_stride, 0, 2, 0x18c5ae);
        test_imgset_val!(baseline_imgsets[0], 4, img_stride, 0, 3, 0x1cc5ae);

        test_imgset_val!(baseline_imgsets[0], 5, img_stride, 0, 0, 0x10c5af);
        test_imgset_val!(baseline_imgsets[0], 5, img_stride, 0, 1, 0x14c5af);
        test_imgset_val!(baseline_imgsets[0], 5, img_stride, 0, 2, 0x18c5af);
        test_imgset_val!(baseline_imgsets[0], 5, img_stride, 0, 3, 0x1cc5af);

        test_imgset_val!(baseline_imgsets[0], 6, img_stride, 0, 0, 0x10bec6);
        test_imgset_val!(baseline_imgsets[0], 6, img_stride, 0, 1, 0x14bec6);
        test_imgset_val!(baseline_imgsets[0], 6, img_stride, 0, 2, 0x18bec6);
        test_imgset_val!(baseline_imgsets[0], 6, img_stride, 0, 3, 0x1cbec6);

        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 0, 0x10bec7);
        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 1, 0x14bec7);
        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 2, 0x18bec7);
        test_imgset_val!(baseline_imgsets[0], 7, img_stride, 0, 3, 0x1cbec7);

        /* ... */

        test_imgset_val!(baseline_imgsets[5], 0, img_stride, 0, 0, 0x10f1ce);
        test_imgset_val!(baseline_imgsets[5], 0, img_stride, 0, 1, 0x14f1ce);
        test_imgset_val!(baseline_imgsets[5], 0, img_stride, 0, 2, 0x18f1ce);
        test_imgset_val!(baseline_imgsets[5], 0, img_stride, 0, 3, 0x1cf1ce);

        test_imgset_val!(baseline_imgsets[5], 1, img_stride, 0, 0, -0x10f1cf);
        test_imgset_val!(baseline_imgsets[5], 1, img_stride, 0, 1, -0x14f1cf);
        test_imgset_val!(baseline_imgsets[5], 1, img_stride, 0, 2, -0x18f1cf);
        test_imgset_val!(baseline_imgsets[5], 1, img_stride, 0, 3, -0x1cf1cf);

        test_imgset_val!(baseline_imgsets[5], 2, img_stride, 0, 0, 0x10ea26);
        test_imgset_val!(baseline_imgsets[5], 2, img_stride, 0, 1, 0x14ea26);
        test_imgset_val!(baseline_imgsets[5], 2, img_stride, 0, 2, 0x18ea26);
        test_imgset_val!(baseline_imgsets[5], 2, img_stride, 0, 3, 0x1cea26);

        test_imgset_val!(baseline_imgsets[5], 3, img_stride, 0, 0, -0x10ea27);
        test_imgset_val!(baseline_imgsets[5], 3, img_stride, 0, 1, -0x14ea27);
        test_imgset_val!(baseline_imgsets[5], 3, img_stride, 0, 2, -0x18ea27);
        test_imgset_val!(baseline_imgsets[5], 3, img_stride, 0, 3, -0x1cea27);

        test_imgset_val!(baseline_imgsets[5], 4, img_stride, 0, 0, 0x10f1be);
        test_imgset_val!(baseline_imgsets[5], 4, img_stride, 0, 1, 0x14f1be);
        test_imgset_val!(baseline_imgsets[5], 4, img_stride, 0, 2, 0x18f1be);
        test_imgset_val!(baseline_imgsets[5], 4, img_stride, 0, 3, 0x1cf1be);

        test_imgset_val!(baseline_imgsets[5], 5, img_stride, 0, 0, -0x10f1bf);
        test_imgset_val!(baseline_imgsets[5], 5, img_stride, 0, 1, -0x14f1bf);
        test_imgset_val!(baseline_imgsets[5], 5, img_stride, 0, 2, -0x18f1bf);
        test_imgset_val!(baseline_imgsets[5], 5, img_stride, 0, 3, -0x1cf1bf);

        test_imgset_val!(baseline_imgsets[5], 6, img_stride, 0, 0, 0x10ea16);
        test_imgset_val!(baseline_imgsets[5], 6, img_stride, 0, 1, 0x14ea16);
        test_imgset_val!(baseline_imgsets[5], 6, img_stride, 0, 2, 0x18ea16);
        test_imgset_val!(baseline_imgsets[5], 6, img_stride, 0, 3, 0x1cea16);

        test_imgset_val!(baseline_imgsets[5], 7, img_stride, 0, 0, -0x10ea17);
        test_imgset_val!(baseline_imgsets[5], 7, img_stride, 0, 1, -0x14ea17);
        test_imgset_val!(baseline_imgsets[5], 7, img_stride, 0, 2, -0x18ea17);
        test_imgset_val!(baseline_imgsets[5], 7, img_stride, 0, 3, -0x1cea17);
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

    #[test]
    fn test_cable_length_corrections_mwax() {
        let context = get_mwax_context();
        let width = context.num_timesteps;
        let stride = (((width - 1) / 8) + 1) * 8;

        let coarse_chans = context.coarse_chans.clone();
        let fine_chan_width_hz = context.metafits_context.corr_fine_chan_width_hz.clone();
        let fine_chans_per_coarse = context
            .metafits_context
            .num_corr_fine_chans_per_coarse
            .clone();

        // A vector of all fine channel frequencies for all coarse channels in ascending order.
        let all_freqs_hz: Vec<u32> = coarse_chans
            .iter()
            .flat_map(|coarse_chan| {
                let chan_start_hz = coarse_chan.chan_start_hz.clone();
                (0..fine_chans_per_coarse)
                    .map(move |fine_chan_idx| {
                        chan_start_hz + (fine_chan_idx as u32 * fine_chan_width_hz)
                    })
                    .clone()
            })
            .collect();

        let baseline_imgsets = unsafe {
            let aoflagger = cxx_aoflagger_new();
            context_to_baseline_imgsets(&aoflagger, &context)
        };

        let viz_0_xx_0_0_re = baseline_imgsets[0].ImageBuffer(0)[0 * stride + 0];
        let viz_0_xx_0_0_im = baseline_imgsets[0].ImageBuffer(1)[0 * stride + 0];
        let viz_1_xx_0_0_re = baseline_imgsets[1].ImageBuffer(0)[0 * stride + 0];
        let viz_1_xx_0_0_im = baseline_imgsets[1].ImageBuffer(1)[0 * stride + 0];
        let viz_1_xy_0_0_re = baseline_imgsets[1].ImageBuffer(2)[0 * stride + 0];
        let viz_1_xy_0_0_im = baseline_imgsets[1].ImageBuffer(3)[0 * stride + 0];
        let viz_1_yx_0_0_re = baseline_imgsets[1].ImageBuffer(4)[0 * stride + 0];
        let viz_1_yx_0_0_im = baseline_imgsets[1].ImageBuffer(5)[0 * stride + 0];
        let viz_1_yy_0_0_re = baseline_imgsets[1].ImageBuffer(6)[0 * stride + 0];
        let viz_1_yy_0_0_im = baseline_imgsets[1].ImageBuffer(7)[0 * stride + 0];
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

        let baseline_imgsets = correct_cable_lengths(&context, baseline_imgsets).unwrap();

        // there should be no difference in baseline 0
        test_imgset_val!(baseline_imgsets[0], 0, stride, 0, 0, viz_0_xx_0_0_re);
        test_imgset_val!(baseline_imgsets[0], 1, stride, 0, 0, viz_0_xx_0_0_im);

        ////
        // baseline 1 should be rotated
        ////
        // baseline 1, pol xx, cc 0, fc 0, ts 0
        let rot_1_xx_0_0_re = (cos_1_xx_0 * viz_1_xx_0_0_re - sin_1_xx_0 * viz_1_xx_0_0_im) as f32;
        let rot_1_xx_0_0_im = (sin_1_xx_0 * viz_1_xx_0_0_re + cos_1_xx_0 * viz_1_xx_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[1], 0, stride, 0, 0, rot_1_xx_0_0_re);
        test_imgset_val!(baseline_imgsets[1], 1, stride, 0, 0, rot_1_xx_0_0_im);
        // baseline 1, pol xy, cc 0, fc 0, ts 0
        let rot_1_xy_0_0_re = (cos_1_xy_0 * viz_1_xy_0_0_re - sin_1_xy_0 * viz_1_xy_0_0_im) as f32;
        let rot_1_xy_0_0_im = (sin_1_xy_0 * viz_1_xy_0_0_re + cos_1_xy_0 * viz_1_xy_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[1], 2, stride, 0, 0, rot_1_xy_0_0_re);
        test_imgset_val!(baseline_imgsets[1], 3, stride, 0, 0, rot_1_xy_0_0_im);
        // baseline 1, pol yx, cc 0, fc 0, ts 0
        let rot_1_yx_0_0_re = (cos_1_yx_0 * viz_1_yx_0_0_re - sin_1_yx_0 * viz_1_yx_0_0_im) as f32;
        let rot_1_yx_0_0_im = (sin_1_yx_0 * viz_1_yx_0_0_re + cos_1_yx_0 * viz_1_yx_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[1], 4, stride, 0, 0, rot_1_yx_0_0_re);
        test_imgset_val!(baseline_imgsets[1], 5, stride, 0, 0, rot_1_yx_0_0_im);
        // baseline 1, pol yy, cc 0, fc 0, ts 0
        let rot_1_yy_0_0_re = (cos_1_yy_0 * viz_1_yy_0_0_re - sin_1_yy_0 * viz_1_yy_0_0_im) as f32;
        let rot_1_yy_0_0_im = (sin_1_yy_0 * viz_1_yy_0_0_re + cos_1_yy_0 * viz_1_yy_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[1], 6, stride, 0, 0, rot_1_yy_0_0_re);
        test_imgset_val!(baseline_imgsets[1], 7, stride, 0, 0, rot_1_yy_0_0_im);
        // baseline 1, pol yy, cc 1, fc 1, ts 1
        let rot_1_yy_3_3_re = (cos_1_yy_3 * viz_1_yy_3_3_re - sin_1_yy_3 * viz_1_yy_3_3_im) as f32;
        let rot_1_yy_3_3_im = (sin_1_yy_3 * viz_1_yy_3_3_re + cos_1_yy_3 * viz_1_yy_3_3_im) as f32;
        test_imgset_val!(baseline_imgsets[1], 6, stride, 3, 3, rot_1_yy_3_3_re);
        test_imgset_val!(baseline_imgsets[1], 7, stride, 3, 3, rot_1_yy_3_3_im);
    }

    #[test]
    fn test_cable_length_corrections_ord() {
        let context = get_mwa_ord_context();

        let width = context.num_timesteps;
        let stride = (((width - 1) / 8) + 1) * 8;

        let coarse_chans = context.coarse_chans.clone();
        let fine_chan_width_hz = context.metafits_context.corr_fine_chan_width_hz.clone();
        let fine_chans_per_coarse = context
            .metafits_context
            .num_corr_fine_chans_per_coarse
            .clone();

        // A vector of all fine channel frequencies for all coarse channels in ascending order.
        let all_freqs_hz: Vec<u32> = coarse_chans
            .iter()
            .flat_map(|coarse_chan| {
                let chan_start_hz = coarse_chan.chan_start_hz.clone();
                (0..fine_chans_per_coarse)
                    .map(move |fine_chan_idx| {
                        chan_start_hz + (fine_chan_idx as u32 * fine_chan_width_hz)
                    })
                    .clone()
            })
            .collect();

        let baseline_imgsets = unsafe {
            let aoflagger = cxx_aoflagger_new();
            context_to_baseline_imgsets(&aoflagger, &context)
        };

        let viz_0_xx_0_0_re = baseline_imgsets[0].ImageBuffer(0)[0 * stride + 0];
        let viz_0_xx_0_0_im = baseline_imgsets[0].ImageBuffer(1)[0 * stride + 0];
        let viz_5_xx_0_0_re = baseline_imgsets[5].ImageBuffer(0)[0 * stride + 0];
        let viz_5_xx_0_0_im = baseline_imgsets[5].ImageBuffer(1)[0 * stride + 0];
        let viz_5_xy_0_0_re = baseline_imgsets[5].ImageBuffer(2)[0 * stride + 0];
        let viz_5_xy_0_0_im = baseline_imgsets[5].ImageBuffer(3)[0 * stride + 0];
        let viz_5_yx_0_0_re = baseline_imgsets[5].ImageBuffer(4)[0 * stride + 0];
        let viz_5_yx_0_0_im = baseline_imgsets[5].ImageBuffer(5)[0 * stride + 0];
        let viz_5_yy_0_0_re = baseline_imgsets[5].ImageBuffer(6)[0 * stride + 0];
        let viz_5_yy_0_0_im = baseline_imgsets[5].ImageBuffer(7)[0 * stride + 0];
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

        let baseline_imgsets = correct_cable_lengths(&context, baseline_imgsets).unwrap();

        // there should be no difference in baseline 0
        test_imgset_val!(baseline_imgsets[0], 0, stride, 0, 0, viz_0_xx_0_0_re);
        test_imgset_val!(baseline_imgsets[0], 1, stride, 0, 0, viz_0_xx_0_0_im);

        ////
        // baseline 1 should be rotated
        ////
        // baseline 1, pol xx, cc 0, fc 0, ts 0
        let rot_5_xx_0_0_re = (cos_5_xx_0 * viz_5_xx_0_0_re - sin_5_xx_0 * viz_5_xx_0_0_im) as f32;
        let rot_5_xx_0_0_im = (sin_5_xx_0 * viz_5_xx_0_0_re + cos_5_xx_0 * viz_5_xx_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[5], 0, stride, 0, 0, rot_5_xx_0_0_re);
        test_imgset_val!(baseline_imgsets[5], 1, stride, 0, 0, rot_5_xx_0_0_im);
        // baseline 1, pol xy, cc 0, fc 0, ts 0
        let rot_5_xy_0_0_re = (cos_5_xy_0 * viz_5_xy_0_0_re - sin_5_xy_0 * viz_5_xy_0_0_im) as f32;
        let rot_5_xy_0_0_im = (sin_5_xy_0 * viz_5_xy_0_0_re + cos_5_xy_0 * viz_5_xy_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[5], 2, stride, 0, 0, rot_5_xy_0_0_re);
        test_imgset_val!(baseline_imgsets[5], 3, stride, 0, 0, rot_5_xy_0_0_im);
        // baseline 1, pol yx, cc 0, fc 0, ts 0
        let rot_5_yx_0_0_re = (cos_5_yx_0 * viz_5_yx_0_0_re - sin_5_yx_0 * viz_5_yx_0_0_im) as f32;
        let rot_5_yx_0_0_im = (sin_5_yx_0 * viz_5_yx_0_0_re + cos_5_yx_0 * viz_5_yx_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[5], 4, stride, 0, 0, rot_5_yx_0_0_re);
        test_imgset_val!(baseline_imgsets[5], 5, stride, 0, 0, rot_5_yx_0_0_im);
        // baseline 1, pol yy, cc 0, fc 0, ts 0
        let rot_5_yy_0_0_re = (cos_5_yy_0 * viz_5_yy_0_0_re - sin_5_yy_0 * viz_5_yy_0_0_im) as f32;
        let rot_5_yy_0_0_im = (sin_5_yy_0 * viz_5_yy_0_0_re + cos_5_yy_0 * viz_5_yy_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[5], 6, stride, 0, 0, rot_5_yy_0_0_re);
        test_imgset_val!(baseline_imgsets[5], 7, stride, 0, 0, rot_5_yy_0_0_im);
        // baseline 1, pol yy, cc 1, fc 1, ts 1
        let rot_5_yy_3_3_re = (cos_5_yy_3 * viz_5_yy_3_3_re - sin_5_yy_3 * viz_5_yy_3_3_im) as f32;
        let rot_5_yy_3_3_im = (sin_5_yy_3 * viz_5_yy_3_3_re + cos_5_yy_3 * viz_5_yy_3_3_im) as f32;
        test_imgset_val!(baseline_imgsets[5], 6, stride, 3, 3, rot_5_yy_3_3_re);
        test_imgset_val!(baseline_imgsets[5], 7, stride, 3, 3, rot_5_yy_3_3_im);
    }
}
