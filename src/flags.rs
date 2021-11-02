//! Methods for manipulating flagmasks and flagging imagesets

use std::ops::Range;

use crate::{
    error::{
        BirliError,
        BirliError::{NoCommonTimesteps, NoProvidedTimesteps},
    },
    io::error::IOError,
    FlagFileSet,
};
use cfg_if::cfg_if;
use log::trace;
use mwa_rust_core::mwalib::CorrelatorContext;
use ndarray::{Array2, Array3, ArrayView3, Axis, Zip};
use rayon::prelude::*;

cfg_if! {
    if #[cfg(feature = "aoflagger")] {
        use crate::{
            flag_baseline_view_to_flagmask, jones_baseline_view_to_imageset
        };
        use aoflagger_sys::{CxxAOFlagger, CxxFlagMask, UniquePtr, flagmask_or,
            flagmask_set};
        use indicatif::{ProgressBar, ProgressStyle};
        use mwa_rust_core::Jones;
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

/// Initialize an array of `bool` flags for each timestep, channel provided.
///
/// Assumptions:
/// - you want all baselines
///
/// TODO: use:
/// - [x] antenna_flags
/// - [ ] coarse_chan_flags
/// - [ ] fine_chan_flags
/// - [ ] timestep_flags
/// to set correlator mask
pub fn init_flag_array(
    context: &CorrelatorContext,
    mwalib_timestep_range: &Range<usize>,
    mwalib_coarse_chan_range: &Range<usize>,
    antenna_flags: Option<Vec<bool>>,
) -> Array3<bool> {
    let num_timesteps = mwalib_timestep_range.len();
    let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
    let num_coarse_chans = mwalib_coarse_chan_range.len();
    let num_chans = num_coarse_chans * fine_chans_per_coarse;
    let num_baselines = context.metafits_context.num_baselines;

    let shape = (num_timesteps, num_chans, num_baselines);

    let correlator_baseline_flags: Array2<bool> =
        Array2::from_elem((num_timesteps, num_chans), false);
    // TODO: there are other things that can affect correlator flags
    let flagged_baseline_flags: Array2<bool> = Array2::from_elem((num_timesteps, num_chans), true);

    let antenna_flags = match antenna_flags {
        Some(antenna_flags) => antenna_flags,
        None => (0..context.metafits_context.num_ants)
            .map(|_| false)
            .collect(),
    };

    let mut flag_array = Array3::from_elem(shape, false);

    flag_array
        .axis_iter_mut(Axis(2))
        .into_par_iter()
        .zip(get_baseline_flags(context, antenna_flags))
        .for_each(|(mut baseline_flag_view, baseline_flag_val)| {
            baseline_flag_view.assign(if baseline_flag_val {
                &flagged_baseline_flags
            } else {
                &correlator_baseline_flags
            });
        });

    flag_array
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
/// let (jones_array, flag_array) = context_to_jones_array(
///     &context,
///     &img_timestep_range,
///     &img_coarse_chan_range,
///     None
/// );
///
/// // use the default strategy file location for MWA
/// let strategy_filename = &aoflagger.FindStrategyFileMWA();
///
/// // run the strategy on the imagesets, and get the resulting flagmasks for each baseline
/// let flag_array = flag_jones_array_existing(
///    &aoflagger,
///    &strategy_filename,
///    &jones_array,
///    Some(flag_array),
///    true,
/// );
/// ```
#[cfg(feature = "aoflagger")]
pub fn flag_jones_array_existing(
    aoflagger: &CxxAOFlagger,
    strategy_filename: &str,
    jones_array: &Array3<Jones<f32>>,
    flag_array: Option<Array3<bool>>,
    re_apply_existing: bool,
) -> Array3<bool> {
    trace!("start flag_jones_array");

    let jones_shape = jones_array.dim();

    let mut flags_provided = false;
    let mut flag_array = match flag_array {
        Some(flag_array) => {
            let flag_shape = flag_array.dim();
            assert_eq!(jones_shape, flag_shape);
            flags_provided = true;
            flag_array
        }
        None => Array3::from_elem(jones_shape, false),
    };

    // The total reading progress.
    let flag_progress = ProgressBar::new(jones_shape.2 as _)
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
            let mut flagmask: UniquePtr<CxxFlagMask>;
            let flag_strategy = aoflagger.LoadStrategyFile(&strategy_filename.to_string());
            let flag_baseline_view_immutable = flag_baseine_view.view();
            if flags_provided {
                // This lets us pass in our mutable flag array view to something not expecting a mutable.
                flagmask = flag_baseline_view_to_flagmask(aoflagger, &flag_baseline_view_immutable)
                    .unwrap();
                let new_flagmask = flag_strategy.RunExisting(&imgset, &flagmask);

                if re_apply_existing {
                    flagmask_or(&mut flagmask, &new_flagmask);
                } else {
                    flagmask_set(&mut flagmask, &new_flagmask);
                }
            } else {
                flagmask = flag_strategy.Run(&imgset);
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

    flag_array
}

/// Shorthand for [`flag_jones_array_existing`] with `flag_array` as None.
#[cfg(feature = "aoflagger")]
pub fn flag_jones_array(
    aoflagger: &CxxAOFlagger,
    strategy_filename: &str,
    jones_array: &Array3<Jones<f32>>,
) -> Array3<bool> {
    flag_jones_array_existing(aoflagger, strategy_filename, jones_array, None, false)
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
/// use birli::{context_to_jones_array, write_flags, mwalib::CorrelatorContext, init_flag_array, get_antenna_flags};
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
///     Some(get_antenna_flags(&context)),
/// );
///
/// // read visibilities out of the gpubox files
/// let (jones_array, flag_array) = context_to_jones_array(
///     &context,
///     &img_timestep_range,
///     &img_coarse_chan_range,
///     Some(flag_array)
/// );
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
        mwa_rust_core::mwalib::CorrelatorContext,
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

        let mut flag_array =
            init_flag_array(&context, &img_timestep_range, &img_coarse_chan_range, None);

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
pub fn get_weight_factor(
    context: &CorrelatorContext,
) -> f64 {
    let integration_time_s = context.metafits_context.corr_int_time_ms as f64 / 1000.0;
    let fine_chan_width_hz = context.metafits_context.corr_fine_chan_width_hz as f64;
    fine_chan_width_hz * integration_time_s / 10000.0
}

/// Convert the given ndarray of boolean flags to an ndarray of float weights
/// TODO: deal with weight factor when doing averaging.
/// TODO: deal with passband gains
pub fn flag_to_weight_array(
    flag_array: ArrayView3<bool>,
    weight_factor: f64,
) -> Array3<f32> {
    let mut weight_array = Array3::from_elem(flag_array.dim(), 0.0_f32);

    Zip::from(&mut weight_array)
        .and(&flag_array)
        .par_apply(|weight, &flag| {
            *weight = if flag { -weight_factor } else { weight_factor } as f32;
        });

    weight_array
}

#[cfg(test)]
#[cfg(feature = "aoflagger")]
/// Tests which require the use of the aoflagger feature
mod tests_aoflagger {
    use super::{get_flaggable_timesteps, init_flag_array};
    use mwa_rust_core::{mwalib::CorrelatorContext, Complex, Jones};
    use ndarray::Array3;

    use crate::{
        context_to_jones_array,
        flags::{flag_jones_array, flag_jones_array_existing, get_antenna_flags},
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

        let existing_flag_array = flag_jones_array_existing(
            &aoflagger,
            &strategy_filename,
            &jones_array,
            Some(existing_flag_array),
            true,
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

        let (jones_array, _) =
            context_to_jones_array(&context, &img_timestep_range, &img_coarse_chan_range, None);

        let existing_flag_array = init_flag_array(
            &context,
            &img_timestep_range,
            &img_coarse_chan_range,
            Some(get_antenna_flags(&context)),
        );

        let strategy_filename = &aoflagger.FindStrategyFileMWA();

        let existing_flag_array = flag_jones_array_existing(
            &aoflagger,
            strategy_filename,
            &jones_array,
            Some(existing_flag_array),
            true,
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
