//! Methods for manipulating flagmasks and flagging imagesets

use std::ops::Range;

use crate::{
    error::{
        BirliError,
        BirliError::{NoCommonCoarseChans, NoCommonTimesteps, NoProvidedTimesteps},
    },
    io::error::IOError,
    ndarray::Axis,
    ndarray::{Array, Array3, Array4, ArrayView, ArrayView3, Dimension},
    FlagFileSet,
};
use cfg_if::cfg_if;
use log::trace;
use marlu::mwalib::CorrelatorContext;

cfg_if! {
    if #[cfg(feature = "aoflagger")] {
        use crate::{
            flag_baseline_view_to_flagmask, jones_baseline_view_to_imageset,
        };
        use aoflagger_sys::{CxxAOFlagger, flagmask_or,
            flagmask_set};
        use indicatif::{ProgressBar, ProgressStyle};
        use marlu::{Jones, rayon::prelude::*};
    }
}

/// Produce a vector of timesteps which can be used for creating imagesets and
/// flagmasks for aoflagger flagging, given an [`mwalib::CorrelatorContext`].
///
/// These timesteps:
/// - are contiguous, and are each separated by an integration time
/// - start at the latest start time of visibilities in all provided gpubox files
///   (AKA the first common timestep)
/// - end at the latest provided visibility
///
/// Start times are determined from the TIME and MILLITIM headers of individual
/// gpubox visibility HDUs
///
/// [`mwalib::CorrelatorContext`]: https://docs.rs/mwalib/latest/mwalib/struct.CorrelatorContext.html
/// mwalib core concepts: https://github.com/MWATelescope/mwalib/wiki/Key-Concepts#timesteps-and-coarse-channels
///
/// # Examples
///
/// ```rust
/// use birli::{get_flaggable_timesteps, mwalib::CorrelatorContext};
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
        (Some(&first), Some(&last)) if first <= last => Ok((first..last + 1).collect()),
        (.., None) => Err(NoProvidedTimesteps {
            hdu_info: format!("{:?}", &context.gpubox_time_map),
        }),
        _ => Err(NoCommonTimesteps {
            hdu_info: format!("{:?}", &context.gpubox_time_map),
        }),
    }
}

/// Produce a vector of flags for antennas whose imputs are flagged in `context`
///
/// # Examples
///
/// ```rust
/// use birli::{get_antenna_flags, mwalib::CorrelatorContext};
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

/// Produce a contiguous range of coarse channel indices from [`MWALib::CorrelatorContext.common_coarse_chan_indices`]
///
/// # Errors
///
/// - Can throw NoCommonCoarseChannels if the common coarse channel range is empty.
///
pub fn get_coarse_chan_range(context: &CorrelatorContext) -> Result<Range<usize>, BirliError> {
    let common_indices = &context.common_coarse_chan_indices;

    match (common_indices.first(), common_indices.last()) {
        (Some(first), Some(last)) if first <= last => Ok((*first)..(*last + 1)),
        _ => Err(NoCommonCoarseChans {
            hdu_info: format!("{:?}", &context.gpubox_time_map),
        }),
    }
}
/// Produce a vector of flags for all known coarse channels.
///
/// A channel is flagged if it appears in the metafits, but is not provided.
///
pub fn get_coarse_chan_flags(context: &CorrelatorContext) -> Vec<bool> {
    let provided_indices = &context.provided_coarse_chan_indices;

    let mut flags = vec![false; context.num_coarse_chans];
    for (idx, flag) in flags.iter_mut().enumerate() {
        if !provided_indices.contains(&idx) {
            *flag = true;
        }
    }

    flags
}

/// Produce a contiguous range of timestep indices from the first [`MWALib::CorrelatorContext.common_timestep_indices`]
/// to the last [`MWALib::CorrelatorContext.provided_timestep_indices`]
///
/// # Errors
///
/// - can throw NoProvidedTimesteps if no timesteps have been provided.
/// - can throw NoCommonTimesteps if there are no common timesteps.
pub fn get_timestep_range(context: &CorrelatorContext) -> Result<Range<usize>, BirliError> {
    let common_indices = &context.common_timestep_indices;
    let provided_indices = &context.provided_timestep_indices;
    match (common_indices.first(), provided_indices.last()) {
        (Some(first), Some(last)) if first <= last => Ok((*first)..(*last + 1)),
        (.., None) => Err(NoProvidedTimesteps {
            hdu_info: format!("{:?}", &context.gpubox_time_map),
        }),
        _ => Err(NoCommonTimesteps {
            hdu_info: format!("{:?}", &context.gpubox_time_map),
        }),
    }
}

/// Produce a vector of flags for all known timesteps.
///
/// A timestep is flagged if it appears in the metafits, but is not provided.
pub fn get_timestep_flags(context: &CorrelatorContext, flag_idxs: Option<Vec<usize>>) -> Vec<bool> {
    let provided_indices = &context.provided_timestep_indices;

    let mut flags = vec![false; context.num_timesteps];
    for (idx, flag) in flags.iter_mut().enumerate() {
        if !provided_indices.contains(&idx) {
            *flag = true;
        }
    }

    for idx in flag_idxs.unwrap_or_default() {
        flags[idx] = true;
    }

    flags
}

/// Produce a vector of flags for baslines whose antennae are flagged in `antenna_flags`
pub fn get_baseline_flags(context: &CorrelatorContext, antenna_flags: &[bool]) -> Vec<bool> {
    assert_eq!(antenna_flags.len(), context.metafits_context.num_ants);
    context
        .metafits_context
        .baselines
        .clone()
        .into_iter()
        .map(|baseline| antenna_flags[baseline.ant1_index] || antenna_flags[baseline.ant2_index])
        .collect()
}

/// Produce a vector of flags for fine each channel in the selected range of channels.
// pub fn get_chan_flags(
//     context: &CorrelatorContext,
//     coarse_chan_flags: Option<&[bool]>,
//     fine_chan_flags: Option<&[bool]>,
// ) -> Vec<bool> {

//     let num_coarse_chans = context.num_coarse_chans;
//     let num_fine_chans = context.metafits_context.num_corr_fine_chans_per_coarse;

//     let coarse_chan_flags: Vec<_> = match coarse_chan_flags {
//         Some(coarse_chan_flags) => {
//             assert_eq!(coarse_chan_flags.len(), num_coarse_chans);
//             coarse_chan_flags.to_vec()
//         },
//         None => (0..num_coarse_chans).into_iter()
//             .map(|_| false)
//             .collect(),
//     };
//     let fine_chan_flags: Vec<_> = match fine_chan_flags {
//         Some(fine_chan_flags) => {
//             assert_eq!(fine_chan_flags.len(), num_fine_chans);
//             fine_chan_flags.to_vec()
//         },
//         None => (0..num_fine_chans)
//             .map(|_| false)
//             .collect(),
//     };

//     coarse_chan_flags.iter().flat_map(|coarse_chan_flag| {
//         if *coarse_chan_flag {
//             fine_chan_flags.to_vec()
//         } else {
//             vec![false; fine_chan_flags.len()]
//         }
//     }).collect()
// }

/// Initialize a 3D array of `bool` flags for each timestep, channel, baseline provided.
///
/// Optionally provide flags for the timestep, channel and baseline axes
///
/// TODO: fix clippy::too_many_arguments
#[allow(clippy::too_many_arguments)]
pub fn init_flag_array(
    context: &CorrelatorContext,
    mwalib_timestep_range: &Range<usize>,
    mwalib_coarse_chan_range: &Range<usize>,
    timestep_flags: Option<&[bool]>,
    coarse_chan_flags: Option<&[bool]>,
    fine_chan_flags: Option<&[bool]>,
    baseline_flags: Option<&[bool]>,
) -> Array3<bool> {
    let mwalib_timestep_range = mwalib_timestep_range.clone();
    let mwalib_coarse_chan_range = mwalib_coarse_chan_range.clone();
    let num_timesteps_total = context.num_timesteps;
    let num_timesteps = mwalib_timestep_range.len();
    let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
    let num_coarse_chans_total = context.num_coarse_chans;
    let num_coarse_chans = mwalib_coarse_chan_range.len();
    let num_chans = num_coarse_chans * fine_chans_per_coarse;
    let num_baselines = context.metafits_context.num_baselines;
    let shape = (num_timesteps, num_chans, num_baselines);

    let timestep_flags: Vec<_> = match timestep_flags {
        Some(timestep_flags) => {
            assert_eq!(timestep_flags.len(), num_timesteps_total);
            timestep_flags[mwalib_timestep_range].to_vec()
        }
        None => (0..num_timesteps).map(|_| false).collect(),
    };
    let coarse_chan_flags: Vec<_> = match coarse_chan_flags {
        Some(coarse_chan_flags) => {
            assert_eq!(coarse_chan_flags.len(), num_coarse_chans_total);
            coarse_chan_flags[mwalib_coarse_chan_range].to_vec()
        }
        None => (0..num_coarse_chans).map(|_| false).collect(),
    };
    let fine_chan_flags: Vec<_> = match fine_chan_flags {
        Some(fine_chan_flags) => {
            assert_eq!(fine_chan_flags.len(), fine_chans_per_coarse);
            fine_chan_flags.to_vec()
        }
        None => (0..fine_chans_per_coarse).map(|_| false).collect(),
    };
    let baseline_flags: Vec<_> = match baseline_flags {
        Some(baseline_flags) => {
            assert_eq!(baseline_flags.len(), num_baselines);
            baseline_flags.to_vec()
        }
        None => (0..num_baselines).map(|_| false).collect(),
    };

    let chan_flags: Vec<_> = coarse_chan_flags
        .iter()
        .flat_map(|coarse_chan_flag| {
            if *coarse_chan_flag {
                vec![true; fine_chan_flags.len()]
            } else {
                fine_chan_flags.to_vec()
            }
        })
        .collect();

    Array3::from_shape_fn(shape, |(ts_idx, ch_idx, bl_idx)| {
        timestep_flags[ts_idx] || chan_flags[ch_idx] || baseline_flags[bl_idx]
    })
}

/// Expand an array into a new axis by repeating each element `size` times
pub fn add_dimension<T>(array: ArrayView3<T>, size: usize) -> Array4<T>
where
    T: std::clone::Clone,
{
    let shape = array.dim();
    array
        .insert_axis(Axis(3))
        .broadcast((shape.0, shape.1, shape.2, size))
        .unwrap()
        .to_owned()
}

/// Expand a flag array into a new axis by repeating each element `size` times
#[deprecated(since = "0.4.1", note = "please use `add_dimension` instead")]
pub fn expand_flag_array(array: ArrayView3<bool>, size: usize) -> Array4<bool> {
    add_dimension(array, size)
}

/// Flag an ndarray of [`Jones`] visibilities, given a [`CxxAOFlagger`] instance,
/// a [`CxxStrategy`] filename, returning an [`ndarray::Array3`] of boolean flags.
///
/// Providing some existing flags is optional, however these flags must be the same
/// dimension as the provided Jones array. If these are not provided, an empty flag
/// array is created instead
///
/// if [`re_apply_existing`] is true, then the new flags are binary or'd with
/// the existing flags, otherwise they overwrite them.
///
/// # Performance
///
/// Because of all the memory juggling required to use aoflagger flagmasks,
/// providing existing flagmasks is slower.
///
///
/// # Examples
///
/// ```
/// use birli::{context_to_jones_array, init_flag_array, flag_jones_array_existing, write_flags,
///     get_antenna_flags, get_flaggable_timesteps, mwalib::CorrelatorContext, cxx_aoflagger_new};
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
/// let img_timestep_range =
///     *img_timestep_idxs.first().unwrap()..(*img_timestep_idxs.last().unwrap() + 1);
/// let img_coarse_chan_range =
///     *img_coarse_chan_idxs.first().unwrap()..(*img_coarse_chan_idxs.last().unwrap() + 1);
///
/// // read visibilities out of the gpubox files
/// let (jones_array, mut flag_array) = context_to_jones_array(
///     &context,
///     &img_timestep_range,
///     &img_coarse_chan_range,
///     None,
///     false,
/// ).unwrap();
///
/// // use the default strategy file location for MWA
/// let strategy_filename = &aoflagger.FindStrategyFileMWA();
///
/// // run the strategy on the imagesets, and get the resulting flagmasks for each baseline
/// flag_jones_array_existing(
///    &aoflagger,
///    &strategy_filename,
///    &jones_array,
///    &mut flag_array,
///    true,
///    false,
/// );
/// ```
#[cfg(feature = "aoflagger")]
pub fn flag_jones_array_existing(
    aoflagger: &CxxAOFlagger,
    strategy_filename: &str,
    jones_array: &Array3<Jones<f32>>,
    flag_array: &mut Array3<bool>,
    re_apply_existing: bool,
    draw_progress: bool,
) {
    use indicatif::ProgressDrawTarget;

    trace!("start flag_jones_array");

    let jones_shape = jones_array.dim();

    let draw_target = if draw_progress {
        ProgressDrawTarget::stderr()
    } else {
        ProgressDrawTarget::hidden()
    };

    // The total reading progress.
    let flag_progress = ProgressBar::with_draw_target(jones_shape.2 as _, draw_target)
        .with_style(
            ProgressStyle::default_bar()
                .template(
                    "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
                )
                .progress_chars("=> "),
        )
        .with_position(0)
        .with_message("flagging b'lines");

    jones_array
        .axis_iter(Axis(2))
        .into_par_iter()
        .zip(flag_array.axis_iter_mut(Axis(2)))
        .for_each(|(jones_baseline_view, mut flag_baseine_view)| {
            let imgset = jones_baseline_view_to_imageset(aoflagger, &jones_baseline_view).unwrap();
            let flag_strategy = aoflagger.LoadStrategyFile(&strategy_filename.to_string());
            let flag_baseline_view_immutable = flag_baseine_view.view();
            // This lets us pass in our mutable flag array view to something not expecting a mutable.
            let mut flagmask =
                flag_baseline_view_to_flagmask(aoflagger, &flag_baseline_view_immutable).unwrap();
            let new_flagmask = flag_strategy.RunExisting(&imgset, &flagmask);

            if re_apply_existing {
                flagmask_or(&mut flagmask, &new_flagmask);
            } else {
                flagmask_set(&mut flagmask, &new_flagmask);
            }
            let flag_buf = flagmask.Buffer();
            let stride = flagmask.HorizontalStride();

            // TODO: assign by slice
            for (img_timestep_idx, mut flag_timestep_view) in
                flag_baseine_view.outer_iter_mut().enumerate()
            {
                for (img_chan_idx, mut flag_singular_view) in
                    flag_timestep_view.outer_iter_mut().enumerate()
                {
                    flag_singular_view.fill(flag_buf[img_chan_idx * stride + img_timestep_idx]);
                }
            }
            flag_progress.inc(1);
        });

    flag_progress.finish();
    trace!("end flag_jones_array");
}

/// Shorthand for [`flag_jones_array_existing`] with `flag_array` as None.
#[cfg(feature = "aoflagger")]
pub fn flag_jones_array(
    aoflagger: &CxxAOFlagger,
    strategy_filename: &str,
    jones_array: &Array3<Jones<f32>>,
) -> Array3<bool> {
    let mut flag_array = Array3::from_elem(jones_array.dim(), false);
    flag_jones_array_existing(
        aoflagger,
        strategy_filename,
        jones_array,
        &mut flag_array,
        false,
        false,
    );
    flag_array
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
/// use birli::{
///     context_to_jones_array, write_flags, mwalib::CorrelatorContext,
///     init_flag_array, get_antenna_flags, get_baseline_flags};
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
/// // Determine which timesteps and coarse channels we want to use
/// let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
/// let img_timestep_idxs = &context.common_timestep_indices;
///
/// let img_timestep_range =
///     *img_timestep_idxs.first().unwrap()..(*img_timestep_idxs.last().unwrap() + 1);
/// let img_coarse_chan_range =
///     *img_coarse_chan_idxs.first().unwrap()..(*img_coarse_chan_idxs.last().unwrap() + 1);
///
/// // Prepare our flagmasks with known bad antennae
/// let flag_array = init_flag_array(
///     &context,
///     &img_timestep_range,
///     &img_coarse_chan_range,
///     None,
///     None,
///     None,
///     Some(&get_baseline_flags(&context, &get_antenna_flags(&context))),
/// );
///
/// // read visibilities out of the gpubox files
/// let (jones_array, flag_array) = context_to_jones_array(
///     &context,
///     &img_timestep_range,
///     &img_coarse_chan_range,
///     Some(flag_array),
///     false,
/// ).unwrap();
///
/// // write the flags to disk as .mwaf
/// write_flags(&context, &flag_array, flag_template.to_str().unwrap(), &img_coarse_chan_range).unwrap();
/// ```
///
/// # Errors
///
/// - Will error with [IOError::FitsOpen] if there are files already present at the paths specified in filename template.
/// - Will error with [IOError::InvalidFlagFilenameTemplate] if an invalid flag filename template is provided (wrong number of percents).

pub fn write_flags(
    context: &CorrelatorContext,
    flag_array: &Array3<bool>,
    filename_template: &str,
    mwalib_coarse_chan_range: &Range<usize>,
) -> Result<(), IOError> {
    trace!("start write_flags");

    let gpubox_ids = context.coarse_chans[mwalib_coarse_chan_range.clone()]
        .iter()
        .map(|chan| chan.gpubox_number)
        .collect::<Vec<_>>();

    trace!(
        "writing flags to template: {}, gpubox ids: {:?}",
        filename_template,
        gpubox_ids
    );

    let mut flag_file_set = FlagFileSet::new(filename_template, &gpubox_ids, context.mwa_version)?;
    flag_file_set.write_flag_array(context, flag_array, &gpubox_ids)?;

    trace!("end write_flags");
    Ok(())
}

#[cfg(test)]
mod tests {
    #![allow(clippy::float_cmp)]
    use super::{get_flaggable_timesteps, init_flag_array, write_flags};
    use glob::glob;
    use tempfile::tempdir;

    use crate::{
        error::BirliError::{NoCommonTimesteps, NoProvidedTimesteps},
        marlu::mwalib::CorrelatorContext,
        FlagFileSet,
    };

    // TODO: Why does clippy think CxxImageSet.ImageBuffer() is &[f64]?
    // TODO: deduplicate from lib.rs
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

    #[test]
    fn test_get_flaggable_timesteps_handles_no_overlap() {
        let context = get_mwa_ord_no_overlap_context();

        assert!(matches!(
            get_flaggable_timesteps(&context),
            Err(NoCommonTimesteps { .. })
        ));
    }

    #[test]
    fn test_get_flaggable_timesteps_handles_no_provided() {
        let context = get_mwa_ord_no_timesteps_context();

        assert!(matches!(
            get_flaggable_timesteps(&context),
            Err(NoProvidedTimesteps { .. })
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
    fn test_write_flags_mwax_minimal() {
        let flag_timestep = 1;
        let flag_channel = 1;
        let flag_baseline = 1;

        let context = get_mwax_context();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(img_timestep_idxs.len(), 4);
        let img_timestep_range =
            *img_timestep_idxs.first().unwrap()..(*img_timestep_idxs.last().unwrap() + 1);

        let img_coarse_chan_idxs = &context.common_coarse_chan_indices[..1].to_vec();
        let img_coarse_chan_range =
            *img_coarse_chan_idxs.first().unwrap()..(*img_coarse_chan_idxs.last().unwrap() + 1);

        let mut flag_array = init_flag_array(
            &context,
            &img_timestep_range,
            &img_coarse_chan_range,
            None,
            None,
            None,
            None,
        );

        flag_array[[flag_timestep, flag_channel, flag_baseline]] = true;

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
            &flag_array,
            filename_template.to_str().unwrap(),
            &img_coarse_chan_range,
        )
        .unwrap();

        let flag_files = glob(tmp_dir.path().join("Flagfile*.mwaf").to_str().unwrap()).unwrap();

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
}

/// Get the weight factor of an observation's context.
///
/// This is a conceptfrom Cotter, and the legacy MWA correlator where the value
/// is a multiple of the frequency averaging factor (relative to 10kHz), and the
/// time averaging factor (relative to 1s).
pub fn get_weight_factor(context: &CorrelatorContext) -> f64 {
    let integration_time_s = context.metafits_context.corr_int_time_ms as f64 / 1000.0;
    let fine_chan_width_hz = context.metafits_context.corr_fine_chan_width_hz as f64;
    fine_chan_width_hz * integration_time_s / 10000.0
}

/// Convert the given ndarray of boolean flags to an ndarray of float weights
/// TODO: deal with weight factor when doing averaging.
/// TODO: deal with passband gains
pub fn flag_to_weight_array<D>(flag_array: ArrayView<bool, D>, weight_factor: f64) -> Array<f32, D>
where
    D: Dimension,
{
    Array::from_elem(flag_array.dim(), weight_factor as f32)
}

#[cfg(test)]
#[cfg(feature = "aoflagger")]
/// Tests which require the use of the aoflagger feature
mod tests_aoflagger {
    use super::init_flag_array;
    use marlu::{mwalib::CorrelatorContext, Complex, Jones};
    use ndarray::Array3;

    use crate::{
        context_to_jones_array,
        flags::{
            flag_jones_array, flag_jones_array_existing, get_antenna_flags, get_baseline_flags,
            get_flaggable_timesteps,
        },
    };
    use aoflagger_sys::cxx_aoflagger_new;

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

    /// Test that flags are sane when not using existing flags.
    #[test]
    fn test_flag_jones_array_minimal() {
        let width = 64; // number of timesteps
        let height = 64; // number of channels
        let num_baselines = 2;

        let noise_bl = 1;
        let noise_x = 32;
        let noise_y = 32;
        let noise_z = 1;
        let signal_val = 0 as f32;
        let noise_val = 0xffffff as f32;

        let aoflagger = unsafe { cxx_aoflagger_new() };

        let signal_jones = Jones::from([
            Complex::new(signal_val, signal_val),
            Complex::new(signal_val, signal_val),
            Complex::new(signal_val, signal_val),
            Complex::new(signal_val, signal_val),
        ]);

        let mut jones_array = Array3::from_elem((width, height, num_baselines), signal_jones);
        let noise_jones = jones_array.get_mut((noise_x, noise_y, noise_bl)).unwrap();
        if noise_z % 2 == 0 {
            noise_jones[noise_z / 2].re = noise_val;
        } else {
            noise_jones[noise_z / 2].im = noise_val;
        }

        let strategy_filename = aoflagger.FindStrategyFileGeneric(&String::from("minimal"));

        let flag_array = flag_jones_array(&aoflagger, &strategy_filename, &jones_array);

        assert!(!flag_array.get((0, 0, 0)).unwrap());
        assert!(!flag_array.get((noise_x, noise_y, 0)).unwrap());
        assert!(!flag_array.get((0, 0, noise_bl)).unwrap());
        assert!(flag_array.get((noise_x, noise_y, noise_bl)).unwrap());
    }

    /// Test that a single noise value is preserved when using existing flags
    #[test]
    fn test_flag_jones_array_existing_minimal() {
        let width = 64;
        let height = 64;
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

        let signal_jones = Jones::from([
            Complex::new(signal_val, signal_val),
            Complex::new(signal_val, signal_val),
            Complex::new(signal_val, signal_val),
            Complex::new(signal_val, signal_val),
        ]);

        let mut jones_array = Array3::from_elem((width, height, num_baselines), signal_jones);
        let noise_jones = jones_array.get_mut((noise_x, noise_y, noise_bl)).unwrap();
        if noise_z % 2 == 0 {
            noise_jones[noise_z / 2].re = noise_val;
        } else {
            noise_jones[noise_z / 2].im = noise_val;
        }

        let strategy_filename = aoflagger.FindStrategyFileGeneric(&String::from("minimal"));

        let mut existing_flag_array = Array3::from_elem((width, height, num_baselines), false);
        existing_flag_array[[existing_x, existing_y, existing_bl]] = true;

        flag_jones_array_existing(
            &aoflagger,
            &strategy_filename,
            &jones_array,
            &mut existing_flag_array,
            true,
            false,
        );

        assert!(!existing_flag_array.get((0, 0, 0)).unwrap());
        assert!(!existing_flag_array.get((noise_x, noise_y, 0)).unwrap());
        assert!(existing_flag_array
            .get((existing_x, existing_y, existing_bl))
            .unwrap());
        assert!(!existing_flag_array.get((0, 0, noise_bl)).unwrap());
        assert!(existing_flag_array
            .get((noise_x, noise_y, noise_bl))
            .unwrap());
    }

    /// Test that flagged antennas remain flagged when used as pre-existing flags.
    #[test]
    fn test_flag_jones_array_existing_ord() {
        let aoflagger = unsafe { cxx_aoflagger_new() };

        let context = get_mwa_ord_context();
        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(img_timestep_idxs.len(), 4);
        let img_timestep_range =
            *img_timestep_idxs.first().unwrap()..(*img_timestep_idxs.last().unwrap() + 1);
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let img_coarse_chan_range =
            *img_coarse_chan_idxs.first().unwrap()..(*img_coarse_chan_idxs.last().unwrap() + 1);

        let (jones_array, _) = context_to_jones_array(
            &context,
            &img_timestep_range,
            &img_coarse_chan_range,
            None,
            false,
        )
        .unwrap();

        let mut existing_flag_array = init_flag_array(
            &context,
            &img_timestep_range,
            &img_coarse_chan_range,
            None,
            None,
            None,
            Some(&get_baseline_flags(&context, &get_antenna_flags(&context))),
        );

        let strategy_filename = &aoflagger.FindStrategyFileMWA();

        flag_jones_array_existing(
            &aoflagger,
            strategy_filename,
            &jones_array,
            &mut existing_flag_array,
            true,
            false,
        );

        // test that flagged antennas are indeed flagged

        assert!(!existing_flag_array.get((1, 0, 0)).unwrap());
        assert!(existing_flag_array.get((1, 0, 11)).unwrap());
        assert!(existing_flag_array.get((1, 0, 11)).unwrap());
        assert!(existing_flag_array.get((1, 0, 62)).unwrap());
        assert!(existing_flag_array.get((1, 0, 63)).unwrap());
        assert!(existing_flag_array.get((1, 0, 80)).unwrap());
        assert!(existing_flag_array.get((1, 0, 87)).unwrap());
        assert!(existing_flag_array.get((1, 0, 91)).unwrap());
        assert!(existing_flag_array.get((1, 0, 93)).unwrap());
        assert!(existing_flag_array.get((1, 0, 95)).unwrap());
        assert!(existing_flag_array.get((1, 0, 111)).unwrap());
        assert!(!existing_flag_array.get((1, 0, 113)).unwrap());
    }
}
