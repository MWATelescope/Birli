//! Methods for manipulating flagmasks and flagging imagesets

use crate::{
    error::{
        BirliError,
        BirliError::{NoCommonTimesteps, NoProvidedTimesteps},
    },
    io::error::IOError,
    CxxAOFlagger, CxxFlagMask, CxxImageSet, FlagFileSet,
};
use cxx::UniquePtr;
use indicatif::{ProgressBar, ProgressStyle};
use log::trace;
use mwa_rust_core::mwalib::CorrelatorContext;

use rayon::prelude::*;

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
/// use birli::{get_flaggable_timesteps, mwalib};
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
/// use birli::{get_antenna_flags, cxx_aoflagger_new, mwalib};
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

    get_baseline_flags(context, antenna_flags)
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
/// use birli::{context_to_baseline_imgsets, flag_imgsets, write_flags, cxx_aoflagger_new, mwalib};
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
/// let baseline_flagmasks = flag_imgsets(&aoflagger, strategy_filename, baseline_imgsets);
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
                .Run(imgset);
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
///     get_antenna_flags, get_flaggable_timesteps, mwalib};
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
///     img_coarse_chan_idxs,
///     &img_timestep_idxs,
///     Some(get_antenna_flags(&context)),
/// );
///
/// // generate imagesets for each baseline in the format required by aoflagger
/// let baseline_imgsets = context_to_baseline_imgsets(
///     &aoflagger,
///     &context,
///     img_coarse_chan_idxs,
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
///     strategy_filename,
///     &baseline_imgsets,
///     &mut baseline_flagmasks,
///     true
/// );
/// ```
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
                .RunExisting(imgset, flagmask);
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
/// use birli::{context_to_baseline_imgsets, flag_imgsets, write_flags, cxx_aoflagger_new, mwalib};
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
/// // Determine which timesteps and coarse channels we want to use
/// let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
/// let img_timestep_idxs = &context.common_timestep_indices;
///
/// // generate imagesets for each baseline in the format required by aoflagger
/// let baseline_imgsets = context_to_baseline_imgsets(
///     &aoflagger,
///     &context,
///     img_coarse_chan_idxs,
///     &img_timestep_idxs,
///     None,
/// );
///
/// // use the default strategy file location for MWA
/// let strategy_filename = &aoflagger.FindStrategyFileMWA();
///
/// // run the strategy on the imagesets, and get the resulting flagmasks for each baseline
/// let baseline_flagmasks = flag_imgsets(&aoflagger, strategy_filename, baseline_imgsets);
///
/// // write the flags to disk as .mwaf
/// write_flags(
///     &context,
///     &baseline_flagmasks,
///     flag_template.to_str().unwrap(),
///     img_coarse_chan_idxs
/// );
/// ```
///
/// # Errors
///
/// - Will error with [IOError::FitsOpen] if there are files already present at the paths specified in filename template.
/// - Will error with [IOError::InvalidFlagFilenameTemplate] if an invalid flag filename template is provided (wrong number of percents).

pub fn write_flags(
    context: &CorrelatorContext,
    baseline_flagmasks: &[UniquePtr<CxxFlagMask>],
    filename_template: &str,
    img_coarse_chan_idxs: &[usize],
) -> Result<(), IOError> {
    trace!("start write_flags");

    let gpubox_ids: Vec<usize> = img_coarse_chan_idxs
        .iter()
        .map(|&chan| context.coarse_chans[chan].gpubox_number)
        .collect();

    trace!(
        "writing flags to template: {}, gpubox ids: {:?}",
        filename_template,
        gpubox_ids
    );

    let mut flag_file_set = FlagFileSet::new(filename_template, &gpubox_ids, context.mwa_version)?;
    flag_file_set.write_baseline_flagmasks(context, baseline_flagmasks, img_coarse_chan_idxs)?;

    trace!("end write_flags");
    Ok(())
}

#[cfg(test)]
mod tests {
    #![allow(clippy::float_cmp)]

    use super::{
        flag_imgsets, flag_imgsets_existing, flagmask_or, flagmask_set, get_antenna_flags,
        get_flaggable_timesteps, init_baseline_flagmasks, write_flags,
    };
    use cxx::UniquePtr;
    use glob::glob;
    use mwa_rust_core::mwalib::CorrelatorContext;
    use tempfile::tempdir;

    use crate::{
        context_to_baseline_imgsets, cxx_aoflagger_new,
        error::BirliError::{NoCommonTimesteps, NoProvidedTimesteps},
        CxxAOFlagger, CxxFlagMask, CxxImageSet, FlagFileSet,
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
            img_coarse_chan_idxs,
            &img_timestep_idxs,
            None,
        );

        let strategy_filename = &aoflagger.FindStrategyFileMWA();
        let mut baseline_flagmasks = init_baseline_flagmasks(
            &aoflagger,
            &context,
            img_coarse_chan_idxs,
            &img_timestep_idxs,
            Some(get_antenna_flags(&context)),
        );

        flag_imgsets_existing(
            &aoflagger,
            strategy_filename,
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
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices[..1].to_vec();

        let width = img_timestep_idxs.len();
        let height =
            img_coarse_chan_idxs.len() * context.metafits_context.num_corr_fine_chans_per_coarse;

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
            &baseline_flagmasks,
            filename_template.to_str().unwrap(),
            img_coarse_chan_idxs,
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
