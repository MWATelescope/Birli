//! Methods for manipulating flagmasks and flagging imagesets

use std::ops::Range;

use crate::{
    io::error::IOError,
    marlu::{
        mwalib::{CorrelatorContext, MWAVersion},
        ndarray::prelude::*,
    },
    BirliError, FlagFileSet,
};
use cfg_if::cfg_if;
use derive_builder::Builder;
use itertools::izip;
use log::trace;
use marlu::{io::error::BadArrayShape, VisSelection};

cfg_if! {
    if #[cfg(feature = "aoflagger")] {
        use aoflagger_sys::{CxxAOFlagger, flagmask_or,
            flagmask_set, CxxFlagMask, UniquePtr, CxxImageSet};
        use indicatif::{ProgressBar, ProgressStyle};
        use marlu::{Jones, rayon::prelude::*, ndarray::{ArrayView2, Axis}};
    }
}

/// Which timesteps, channels and baselines are flagged in a given observation
#[derive(Builder, Debug, Default)]
pub struct FlagContext {
    // TODO: remove _flags suffix
    /// Which mwalib timestep indices are flagged
    pub timestep_flags: Vec<bool>,
    /// Which mwalib coarse channel indices are flagged
    pub coarse_chan_flags: Vec<bool>,
    /// Which fine channel indices are flagged in every coarse channel
    pub fine_chan_flags: Vec<bool>,
    /// Which mwalib fine channel indices are flagged
    pub raw_chan_flags: Vec<bool>,
    /// Which mwalib antenna indices are flagged
    pub antenna_flags: Vec<bool>,
    /// Whether auto-correlations are flagged
    #[builder(default = "false")]
    pub autos: bool,
    /// Whether DC (centre) fine channel indices are flagged
    #[builder(default = "false")]
    pub flag_dc: bool,
    /// How many seconds to flag from the start of the observation
    pub flag_init: f32,
    /// How many seconds to flag from the end of the observation
    pub flag_end: f32,
}

impl FlagContext {
    /// Create a new [`FlagContext`] with all flags false, in the given dimensions
    pub fn blank_from_dimensions(
        num_timesteps: usize,
        num_coarse_chans: usize,
        num_fine_chans_per_coarse: usize,
        num_ants: usize,
    ) -> Self {
        Self {
            timestep_flags: vec![false; num_timesteps],
            coarse_chan_flags: vec![false; num_coarse_chans],
            fine_chan_flags: vec![false; num_fine_chans_per_coarse],
            raw_chan_flags: vec![false; num_coarse_chans * num_fine_chans_per_coarse],
            antenna_flags: vec![false; num_ants],
            ..Self::default()
        }
    }

    /// Create a new [`FlagContext`] from a [`marlu::mwalib::CorrelatorContext`], flagging anything not provided.
    ///
    /// - Timesteps are flagged in they are not provided in any of the gpubox files.
    /// - Coarse channels are flagged if they appear in the metafits, but are not provided.
    /// - No fine channel flags are set by default.
    /// -
    ///
    ///
    ///
    /// # Examples
    ///
    /// ```rust
    /// use birli::{FlagContext, mwalib::CorrelatorContext};
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
    /// let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();
    /// let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
    ///
    /// assert!(!flag_ctx.antenna_flags[0]);
    /// assert!(flag_ctx.antenna_flags[63]);
    /// assert!(flag_ctx.coarse_chan_flags[0]);
    /// assert!(!flag_ctx.coarse_chan_flags[23]);
    /// ```
    ///
    /// TODO:
    /// - quack time
    /// - flag side bandwidth
    /// - auto flag centre
    pub fn from_mwalib(corr_ctx: &CorrelatorContext) -> Self {
        let mut result = Self::blank_from_dimensions(
            corr_ctx.num_timesteps,
            corr_ctx.num_coarse_chans,
            corr_ctx.metafits_context.num_corr_fine_chans_per_coarse,
            corr_ctx.metafits_context.num_ants,
        );
        for (i, flag) in result.timestep_flags.iter_mut().enumerate() {
            *flag = !corr_ctx.provided_timestep_indices.contains(&i);
        }

        for (i, flag) in result.coarse_chan_flags.iter_mut().enumerate() {
            *flag = !corr_ctx.provided_coarse_chan_indices.contains(&i);
        }

        for (antenna, flag) in izip!(
            corr_ctx.metafits_context.antennas.iter(),
            result.antenna_flags.iter_mut()
        ) {
            *flag = antenna.rfinput_x.flagged || antenna.rfinput_y.flagged;
        }

        result.flag_dc = matches!(
            corr_ctx.mwa_version,
            MWAVersion::CorrOldLegacy | MWAVersion::CorrLegacy
        );

        result.flag_init = corr_ctx.metafits_context.quack_time_duration_ms as f32 / 1000.0;
        result
    }

    /// Get the mwalib antenna indices which are flagged
    pub fn get_flagged_antenna_idxs(&self) -> Vec<usize> {
        self.antenna_flags
            .iter()
            .enumerate()
            .filter_map(|(i, &flag)| if flag { Some(i) } else { None })
            .collect()
    }

    /// Produce a vector of flags for baslines where either antenna is flagged in `antenna_flags`
    /// or if `autos` is true and it is an autocorrelation.
    pub fn get_baseline_flags(&self, ant_pairs: &[(usize, usize)]) -> Vec<bool> {
        ant_pairs
            .iter()
            .map(|&(ant1, ant2)| {
                self.antenna_flags[ant1] || self.antenna_flags[ant2] || (self.autos && ant1 == ant2)
            })
            .collect()
    }

    /// Apply timestep flags from `flag_init` and `flag_end`
    ///
    /// Default values for these are inferred from metadata, but may be
    /// overridden by command line arguments, so we apply the flags only after
    /// both have been considered.
    ///
    /// TODO: move DC flagging and other flags that combine contextual defaults
    /// and command line arguments here.
    pub fn finalise_flag_settings(&mut self, corr_ctx: &CorrelatorContext) {
        let flag_before = corr_ctx.common_start_unix_time_ms + (self.flag_init * 1000.0) as u64;
        let flag_after = corr_ctx.common_end_unix_time_ms - (self.flag_end * 1000.0) as u64;
        for (flag, timestep) in self.timestep_flags.iter_mut().zip(&corr_ctx.timesteps) {
            let time = timestep.unix_time_ms;
            *flag |= !(time >= flag_before && time < flag_after);
        }
    }

    /// Get the flags for the fine channels in the given coarse channel range.
    ///
    /// This combines the flags from the fine channels, the coarse channels, and
    /// the raw channels.
    pub fn get_raw_chan_flags(&self, coarse_chan_range: &Range<usize>) -> Vec<bool> {
        let coarse_chan_flags = &self.coarse_chan_flags[coarse_chan_range.clone()];
        let fine_chan_count = self.fine_chan_flags.len();
        let mut fine_chan_flags = self.fine_chan_flags.clone();
        if self.flag_dc {
            fine_chan_flags[fine_chan_count / 2] = true;
        }
        let raw_chan_flags = &self.raw_chan_flags
            [coarse_chan_range.start * fine_chan_count..coarse_chan_range.end * fine_chan_count];

        let chan_flags: Vec<_> = coarse_chan_flags
            .iter()
            .enumerate()
            .flat_map(|(local_coarse_idx, coarse_chan_flag)| {
                if *coarse_chan_flag {
                    vec![true; fine_chan_count]
                } else {
                    let raw_start = local_coarse_idx * fine_chan_count;
                    let raw_end = raw_start + fine_chan_count;
                    let coarse_raw_flags = &raw_chan_flags[raw_start..raw_end];

                    // Combine fine_chan_flags with raw_chan_flags for this coarse channel
                    fine_chan_flags
                        .iter()
                        .zip(coarse_raw_flags.iter())
                        .map(|(&fine_flag, &raw_flag)| fine_flag || raw_flag)
                        .collect::<Vec<bool>>()
                }
            })
            .collect();
        chan_flags
    }

    /// Set flags from this context in an existing array.
    ///
    /// # Errors
    ///
    /// Can throw error if array is not the correct shape.
    pub fn set_flags(
        &self,
        mut flag_array: ArrayViewMut3<bool>,
        timestep_range: &Range<usize>,
        coarse_chan_range: &Range<usize>,
        ant_pairs: &[(usize, usize)],
    ) -> Result<(), BirliError> {
        let timestep_flags = &self.timestep_flags[timestep_range.clone()];
        let baseline_flags = self.get_baseline_flags(ant_pairs);

        let chan_flags = self.get_raw_chan_flags(coarse_chan_range);
        let shape = (timestep_range.len(), chan_flags.len(), ant_pairs.len());

        let flag_shape = flag_array.dim();
        if flag_shape.0 > shape.0 || flag_shape.1 > shape.1 || flag_shape.2 > shape.2 {
            return Err(BirliError::BadArrayShape(BadArrayShape {
                argument: "flag_array",
                function: "FlagContext::set_flags",
                expected: format!("dims less than {shape:?}"),
                received: format!("{flag_shape:?}"),
            }));
        };

        flag_array
            .indexed_iter_mut()
            .for_each(|((ts_idx, ch_idx, bl_idx), flag)| {
                *flag = timestep_flags[ts_idx] || chan_flags[ch_idx] || baseline_flags[bl_idx];
            });

        Ok(())
    }
}

/// Create an aoflagger [`CxxImageSet`] for a particular baseline from the given jones array
///
/// # Assumptions
///
/// - `baseline_jones_view` is `[timestep][channel]` for one baseline
/// - imageset is timesteps wide, and channels high
/// - jones matrics are always XX, YY, XY, YX
///
#[cfg(feature = "aoflagger")]
pub fn jones_baseline_view_to_imageset(
    aoflagger: &CxxAOFlagger,
    baseline_jones_view: ArrayView2<Jones<f32>>,
) -> UniquePtr<CxxImageSet> {
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

    imgset
}

/// Create an aoflagger [`CxxImageSet`] for array of amplitudes
///
/// since most strategies assume jones matrix, the imageset will have a count
/// of 2x the polarization count. imaginary values are set to 0.
///
/// # Assumptions
///
/// - `amps_tfp` is `[timestep][channel][pol]` for one baseline
/// - imageset is timesteps wide, and channels high
///
#[cfg(feature = "aoflagger")]
pub fn amps_tfp_to_imageset(
    aoflagger: &CxxAOFlagger,
    amps_tfp: ArrayView3<f32>,
) -> UniquePtr<CxxImageSet> {
    let array_dims = amps_tfp.dim();
    let img_count = array_dims.2 * 2;
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

    for ((t, f, q), &amp) in amps_tfp.indexed_iter() {
        img_bufs[q * 2][f * img_stride + t] = amp;
        img_bufs[q * 2 + 1][f * img_stride + t] = 0.0;
    }

    imgset
}

/// Create an aoflagger [`CxxFlagMask`] for a from the given flag array view
///
/// # Assumptions
///
/// - flag array view is `[timestep][channel]` for one baseline
/// - flagmask is timesteps wide, and channels high
///
#[cfg(feature = "aoflagger")]
pub fn flag_baseline_view_to_flagmask(
    aoflagger: &CxxAOFlagger,
    baseline_flag_view: ArrayView2<bool>,
) -> UniquePtr<CxxFlagMask> {
    let array_dims = baseline_flag_view.dim();
    let mut flag_mask = unsafe { aoflagger.MakeFlagMask(array_dims.0, array_dims.1, false) };
    let stride = flag_mask.HorizontalStride();
    let flag_buf = flag_mask.pin_mut().BufferMut();

    // TODO: assign by slice
    // flag_buf.copy_from_slice(baseline_flag_view.as_slice().unwrap());
    for (img_timestep_idx, timestep_flag_view) in baseline_flag_view.outer_iter().enumerate() {
        for (img_chan_idx, singular_flag_view) in timestep_flag_view.outer_iter().enumerate() {
            flag_buf[img_chan_idx * stride + img_timestep_idx] =
                *singular_flag_view.get(()).unwrap();
        }
    }
    flag_mask
}

/// Flag an ndarray of [`Jones`] visibilities, given a [`CxxAOFlagger`]
/// instance, a [`aoflagger_sys::CxxAOFlagger`] filename, returning an
/// [`ndarray::Array3`](crate::ndarray::Array3) of boolean flags.
///
/// Providing some existing flags is optional, however these flags must be the
/// same dimension as the provided Jones array. If these are not provided, an
/// empty flag array is created instead
///
/// if `re_apply_existing` is true, then the new flags are binary or'd with
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
/// use birli::{FlagContext, flag_jones_array_existing, write_flags,
///     mwalib::CorrelatorContext, cxx_aoflagger_new, VisSelection, io::read_mwalib};
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
/// let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();
///
/// // create a CxxAOFlagger object to perform AOFlagger operations
/// let aoflagger = unsafe { cxx_aoflagger_new() };
///
/// // specify which coarse_chan and timestep indices we want to load into an image.
/// let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
///
/// // Create a blank array to store flags and visibilities
/// let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
/// let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
/// let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
///
/// // read visibilities out of the gpubox files
/// read_mwalib(&vis_sel, &corr_ctx, jones_array.view_mut(), flag_array.view_mut(), false)
///     .unwrap();
///
/// // use the default strategy file location for MWA
/// let strategy_filename = &aoflagger.FindStrategyFileMWA();
///
/// // run the strategy on the imagesets, and get the resulting flagmasks for each baseline
/// flag_jones_array_existing(
///    &aoflagger,
///    &strategy_filename,
///    jones_array.view(),
///    flag_array.view_mut(),
///    true,
///    false,
/// );
/// ```
#[cfg(feature = "aoflagger")]
pub fn flag_jones_array_existing(
    aoflagger: &CxxAOFlagger,
    strategy_filename: &str,
    jones_array: ArrayView3<Jones<f32>>,
    mut flag_array: ArrayViewMut3<bool>,
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
    let flag_progress = ProgressBar::with_draw_target(Some(jones_shape.2 as _), draw_target)
        .with_style(
            ProgressStyle::default_bar()
                .template(
                    "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
                )
                .unwrap()
                .progress_chars("=> "),
        )
        .with_position(0)
        .with_message("flagging b'lines");

    jones_array
        .axis_iter(Axis(2))
        .into_par_iter()
        .zip(flag_array.axis_iter_mut(Axis(2)))
        .for_each(|(jones_baseline_view, mut flag_baseine_view)| {
            let imgset = jones_baseline_view_to_imageset(aoflagger, jones_baseline_view.view());
            let flag_strategy = aoflagger.LoadStrategyFile(&strategy_filename.to_string());
            let flag_baseline_view_immutable = flag_baseine_view.view();
            // This lets us pass in our mutable flag array view to something not expecting a mutable.
            let mut flagmask =
                flag_baseline_view_to_flagmask(aoflagger, flag_baseline_view_immutable.view());
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
    jones_array: ArrayView3<Jones<f32>>,
) -> Array3<bool> {
    let mut flag_array = Array3::from_elem(jones_array.dim(), false);
    flag_jones_array_existing(
        aoflagger,
        strategy_filename,
        jones_array,
        flag_array.view_mut(),
        false,
        false,
    );
    flag_array
}

/// Write flags to disk, given an observation's [`marlu::mwalib::CorrelatorContext`], a vector of
/// [`CxxFlagMask`]s for each baseline in the observation, a filename template and a vector of
/// gpubox IDs.
///
/// The filename template should contain two or 3 percentage (`%`) characters which will be replaced
/// by the gpubox id or channel number (depending on correlator type) provided in `gpubox_ids`. See
/// [`crate::io::mwaf::FlagFileSet`] for more details.
///
/// # Examples
///
/// Here's an example of how to flag some visibility files
///
/// ```rust
/// use birli::{FlagContext, write_flags, mwalib::CorrelatorContext, VisSelection, io::read_mwalib};
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
/// let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();
///
/// // Determine which timesteps and coarse channels we want to use
/// let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
///
/// // Get the GPS time of the first timestep
/// let gps_start = corr_ctx.timesteps[vis_sel.timestep_range.start].gps_time_ms as f64 / 1e3;
///
/// // Prepare our flagmasks with known bad antennae
/// let mut flag_ctx = FlagContext::from_mwalib(&corr_ctx);
///
/// // Create a blank array to store flags and visibilities
/// let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
/// let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
/// flag_ctx.set_flags(
///     flag_array.view_mut(),
///     &vis_sel.timestep_range,
///     &vis_sel.coarse_chan_range,
///     &vis_sel.get_ant_pairs(&corr_ctx.metafits_context)
/// );
/// let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
///
/// // read visibilities out of the gpubox files
/// read_mwalib(&vis_sel, &corr_ctx, jones_array.view_mut(), flag_array.view_mut(), false)
///     .unwrap();
///
/// // write the flags to disk as .mwaf
/// write_flags(flag_template.to_str().unwrap(),
///             &corr_ctx,
///             &vis_sel,
///             flag_array.view(),
///             true,
///             None,
///             None,
/// ).unwrap();
/// ```
///
/// # Errors
///
/// - Will error with [`IOError::FitsOpen`] if there are files already present at the paths specified in filename template.
/// - Will error with [`IOError::InvalidFlagFilenameTemplate`] if an invalid flag filename template is provided (wrong number of percents).
#[allow(clippy::too_many_arguments)]
pub fn write_flags(
    filename_template: &str,
    corr_ctx: &CorrelatorContext,
    vis_sel: &VisSelection,
    flag_array: ArrayView3<bool>,
    draw_progress: bool,
    aoflagger_version: Option<String>,
    aoflagger_strategy: Option<String>,
) -> Result<(), IOError> {
    trace!("start write_flags");

    let gpubox_ids = corr_ctx.coarse_chans[vis_sel.coarse_chan_range.clone()]
        .iter()
        .map(|chan| chan.gpubox_number)
        .collect::<Vec<_>>();

    trace!("writing flags to template: {filename_template}, gpubox ids: {gpubox_ids:?}");

    let mut flag_file_set = FlagFileSet::new(
        filename_template,
        corr_ctx,
        vis_sel,
        aoflagger_version,
        aoflagger_strategy,
    )?;
    flag_file_set.write_flag_array(flag_array, draw_progress)?;
    flag_file_set.finalise()?;

    trace!("end write_flags");
    Ok(())
}

/// Get the weight factor of an observation's `corr_ctx`.
///
/// This is a concept from Cotter, and the legacy MWA correlator where the value
/// is a multiple of the frequency averaging factor (relative to 10kHz), and the
/// time averaging factor (relative to 1s). These factors have been moved to
/// Marlu for better visibility.
pub fn get_weight_factor(corr_ctx: &CorrelatorContext) -> f64 {
    let integration_time_s = corr_ctx.metafits_context.corr_int_time_ms as f64 / 1000.0;
    let fine_chan_width_hz = corr_ctx.metafits_context.corr_fine_chan_width_hz as f64;
    (fine_chan_width_hz / marlu::constants::FREQ_WEIGHT_FACTOR)
        * (integration_time_s / marlu::constants::TIME_WEIGHT_FACTOR)
}

/// Convert the given ndarray of boolean flags to an ndarray of float weights
#[allow(clippy::needless_pass_by_value)]
pub fn flag_to_weight_array<D>(flag_array: ArrayView<bool, D>, weight_factor: f64) -> Array<f32, D>
where
    D: Dimension,
{
    flag_array.map(|f| if *f { -weight_factor } else { weight_factor } as f32)
}

pub(crate) fn get_unflagged_timestep_ranges(flag_chunk: ArrayView3<bool>) -> Vec<Range<usize>> {
    flag_chunk
        .axis_iter(Axis(0))
        .enumerate()
        .filter_map(|(timestep_idx, flag_chunk)| {
            if flag_chunk.iter().any(|&x| !x) {
                Some(timestep_idx)
            } else {
                None
            }
        })
        // now convert this iterator of unflagged timestep ranges to an iterator of contiguous ranges
        .fold(Vec::new(), |mut acc, timestep_idx| {
            if acc.is_empty() {
                acc.push(timestep_idx..(timestep_idx + 1));
            } else {
                let mut last_acc = acc.pop().unwrap();
                if last_acc.end == timestep_idx {
                    last_acc.end = timestep_idx + 1;
                    acc.push(last_acc);
                } else {
                    acc.push(last_acc);
                    acc.push(timestep_idx..(timestep_idx + 1));
                }
            }
            acc
        })
}

#[cfg(test)]
mod tests {
    use super::write_flags;
    use glob::glob;
    use ndarray::Array3;
    use std::ffi::c_char;
    use tempfile::tempdir;

    use crate::{
        flags::get_unflagged_timestep_ranges,
        marlu::selection::SelectionError::{NoCommonTimesteps, NoProvidedTimesteps},
        test_common::{
            get_mwa_ord_dodgy_context, get_mwa_ord_no_overlap_context,
            get_mwa_ord_no_timesteps_context, get_mwax_context,
        },
        FlagFileSet, VisSelection,
    };

    #[test]
    fn test_get_flaggable_timesteps_handles_no_overlap() {
        let corr_ctx = get_mwa_ord_no_overlap_context();
        let vis_sel = VisSelection::from_mwalib(&corr_ctx);

        assert!(matches!(vis_sel, Err(NoCommonTimesteps { .. })));
    }

    #[test]
    fn test_get_flaggable_timesteps_handles_no_provided() {
        let corr_ctx = get_mwa_ord_no_timesteps_context();
        let vis_sel = VisSelection::from_mwalib(&corr_ctx);

        assert!(matches!(vis_sel, Err(NoProvidedTimesteps { .. })));
    }

    #[test]
    fn test_get_flaggable_timesteps_handles_dodgy() {
        let corr_ctx = get_mwa_ord_dodgy_context();
        let vis_sel = VisSelection::from_mwalib(&corr_ctx);

        let timestep_idxs: Vec<_> = vis_sel.unwrap().timestep_range.collect();
        assert_eq!(timestep_idxs.len(), 4);
        assert_eq!(
            timestep_idxs.first(),
            corr_ctx.common_timestep_indices.first()
        );
        assert_eq!(
            timestep_idxs.last(),
            corr_ctx.provided_timestep_indices.last()
        );
    }

    #[test]
    fn test_write_flags_mwax_minimal() {
        let flag_timestep = 1;
        let flag_channel = 1;
        let flag_baseline = 1;

        let corr_ctx = get_mwax_context();
        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        vis_sel.coarse_chan_range =
            vis_sel.coarse_chan_range.start..vis_sel.coarse_chan_range.start + 1;

        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();

        flag_array[[flag_timestep, flag_channel, flag_baseline]] = true;

        let tmp_dir = tempdir().unwrap();

        let gpubox_ids: Vec<usize> = corr_ctx
            .common_coarse_chan_indices
            .iter()
            .map(|&chan| corr_ctx.coarse_chans[chan].gpubox_number)
            .collect();

        let filename_template = tmp_dir.path().join("Flagfile%%%.mwaf");
        let selected_gpuboxes = gpubox_ids[..1].to_vec();

        write_flags(
            filename_template.to_str().unwrap(),
            &corr_ctx,
            &vis_sel,
            flag_array.view(),
            false,
            None,
            None,
        )
        .unwrap();

        let flag_files = glob(tmp_dir.path().join("Flagfile*.mwaf").to_str().unwrap()).unwrap();

        assert_eq!(flag_files.count(), selected_gpuboxes.len());

        let flag_file_set = FlagFileSet::open(
            filename_template.to_str().unwrap(),
            &selected_gpuboxes,
            corr_ctx.mwa_version,
        )
        .unwrap();
        assert_eq!(flag_file_set.gpuboxes.len(), 1);
        assert_eq!(flag_file_set.gpuboxes[0].id, gpubox_ids[0]);
        let flags = flag_file_set.read_flags().unwrap();

        let num_baselines =
            (flag_file_set.header.num_ants * (flag_file_set.header.num_ants + 1)) / 2;
        assert_eq!(
            flag_file_set.header.num_timesteps,
            corr_ctx.num_timesteps as u32
        );
        assert_eq!(
            num_baselines,
            corr_ctx.metafits_context.num_baselines as u32
        );
        assert_eq!(
            flag_file_set.header.num_channels,
            fine_chans_per_coarse as u32
        );
        assert_eq!(
            flags.dim(),
            (
                flag_file_set.header.num_timesteps as usize,
                num_baselines as usize,
                flag_file_set.header.num_channels as usize,
            )
        );

        let tests = [
            (0, 0, 0, c_char::from(false)),
            (0, 0, 1, c_char::from(false)),
            (0, 1, 0, c_char::from(false)),
            (0, 1, 1, c_char::from(false)),
            (0, 2, 0, c_char::from(false)),
            (0, 2, 1, c_char::from(false)),
            (1, 0, 0, c_char::from(false)),
            (1, 0, 1, c_char::from(false)),
            (1, 1, 0, c_char::from(false)),
            (1, 1, 1, c_char::from(true)),
            (1, 2, 0, c_char::from(false)),
            (1, 2, 1, c_char::from(false)),
        ];
        for (timestep_idx, baseline_idx, fine_chan_idx, expected_flag) in tests {
            let row_idx = timestep_idx * num_baselines + baseline_idx;
            let offset = row_idx as usize * fine_chans_per_coarse + fine_chan_idx as usize;
            assert_eq!(
                *flags.iter().nth(offset).unwrap(), expected_flag,
                "with timestep {timestep_idx}, baseline {baseline_idx}, fine_chan {fine_chan_idx}, expected {expected_flag} at row_idx {row_idx}, offset {offset}"
            );
        }
    }

    #[test]
    fn test_get_unflagged_timestep_ranges() {
        let mut flag_array = Array3::from_elem((3, 3, 3), false);
        flag_array[[1, 0, 0]] = true;
        flag_array[[1, 0, 1]] = true;
        flag_array[[1, 0, 2]] = true;
        flag_array[[1, 1, 0]] = true;
        flag_array[[1, 1, 1]] = true;
        flag_array[[1, 1, 2]] = true;
        flag_array[[1, 2, 0]] = true;
        flag_array[[1, 2, 1]] = true;
        flag_array[[1, 2, 2]] = true;
        let unflagged_timestep_ranges = get_unflagged_timestep_ranges(flag_array.view());
        assert_eq!(unflagged_timestep_ranges, vec![0..1, 2..3]);
    }
}

#[cfg(test)]
#[cfg(feature = "aoflagger")]
/// Tests which require the use of the aoflagger feature
mod tests_aoflagger {
    use marlu::{mwalib::CorrelatorContext, Complex, Jones};
    use ndarray::Array3;

    use crate::{
        flags::{flag_jones_array, flag_jones_array_existing, FlagContext},
        io::read_mwalib,
        BirliError, VisSelection,
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
        CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap()
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

        let flag_array = flag_jones_array(&aoflagger, &strategy_filename, jones_array.view());

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
            jones_array.view(),
            existing_flag_array.view_mut(),
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

        let corr_ctx = get_mwa_ord_context();
        let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

        let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        flag_ctx
            .set_flags(
                flag_array.view_mut(),
                &vis_sel.timestep_range,
                &vis_sel.coarse_chan_range,
                &vis_sel.get_ant_pairs(&corr_ctx.metafits_context),
            )
            .unwrap();
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        read_mwalib(
            &vis_sel,
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();

        let strategy_filename = &aoflagger.FindStrategyFileMWA();

        flag_jones_array_existing(
            &aoflagger,
            strategy_filename,
            jones_array.view(),
            flag_array.view_mut(),
            true,
            false,
        );

        // test that flagged antennas are indeed flagged

        assert!(!flag_array.get((1, 0, 0)).unwrap());
        assert!(flag_array.get((1, 0, 11)).unwrap());
        assert!(flag_array.get((1, 0, 11)).unwrap());
        assert!(flag_array.get((1, 0, 62)).unwrap());
        assert!(flag_array.get((1, 0, 63)).unwrap());
        assert!(flag_array.get((1, 0, 80)).unwrap());
        assert!(flag_array.get((1, 0, 87)).unwrap());
        assert!(flag_array.get((1, 0, 91)).unwrap());
        assert!(flag_array.get((1, 0, 93)).unwrap());
        assert!(flag_array.get((1, 0, 95)).unwrap());
        assert!(flag_array.get((1, 0, 111)).unwrap());
        assert!(!flag_array.get((1, 0, 113)).unwrap());
    }

    #[test]
    fn test_set_flags_checks_array_shape() {
        let corr_ctx = get_mwa_ord_context();
        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

        let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();

        // mess with flag_ctx to make it think it's a different shape
        vis_sel.timestep_range = 0..1;
        assert!(matches!(
            flag_ctx.set_flags(
                flag_array.view_mut(),
                &vis_sel.timestep_range,
                &vis_sel.coarse_chan_range,
                &vis_sel.get_ant_pairs(&corr_ctx.metafits_context),
            ),
            Err(BirliError::BadArrayShape(_))
        ));
    }

    #[test]
    fn test_set_flags_with_raw_chan_flags() {
        // Create a simple FlagContext with known dimensions
        let num_timesteps = 2;
        let num_coarse_chans = 3;
        let num_fine_chans_per_coarse = 4;
        let num_ants = 2;

        let mut flag_ctx = FlagContext::blank_from_dimensions(
            num_timesteps,
            num_coarse_chans,
            num_fine_chans_per_coarse,
            num_ants,
        );

        // Set up some test flags
        flag_ctx.timestep_flags[0] = true; // Flag first timestep
        flag_ctx.coarse_chan_flags[1] = true; // Flag second coarse channel
        flag_ctx.fine_chan_flags[2] = true; // Flag third fine channel in all coarse channels

        // Set up raw_chan_flags - flag specific fine channels in specific coarse channels
        // For coarse channel 0, flag fine channel 1
        flag_ctx.raw_chan_flags[1] = true;
        // For coarse channel 2, flag fine channel 3
        flag_ctx.raw_chan_flags[2 * num_fine_chans_per_coarse + 3] = true;

        // Create test ranges
        let timestep_range = 0..num_timesteps;
        let coarse_chan_range = 0..num_coarse_chans;
        let ant_pairs = vec![(0, 1)]; // One baseline

        // Create flag array
        let mut flag_array = Array3::from_elem(
            (
                num_timesteps,
                num_coarse_chans * num_fine_chans_per_coarse,
                1,
            ),
            false,
        );

        // Apply flags
        flag_ctx
            .set_flags(
                flag_array.view_mut(),
                &timestep_range,
                &coarse_chan_range,
                &ant_pairs,
            )
            .unwrap();

        // Test that timestep flags work (first timestep should be flagged)
        assert!(flag_array[[0, 0, 0]]); // First timestep, first channel, first baseline
        assert!(!flag_array[[1, 0, 0]]); // Second timestep, first channel, first baseline

        // Test that coarse channel flags work (second coarse channel should be flagged)
        let second_coarse_start = num_fine_chans_per_coarse;
        let second_coarse_end = second_coarse_start + num_fine_chans_per_coarse;
        for ch in second_coarse_start..second_coarse_end {
            assert!(flag_array[[1, ch, 0]]); // Second timestep, all channels in second coarse channel
        }

        // Test that fine channel flags work (third fine channel should be flagged in all coarse channels)
        for coarse in 0..num_coarse_chans {
            if coarse != 1 {
                // Skip the already flagged coarse channel
                let fine_ch_2_idx = coarse * num_fine_chans_per_coarse + 2;
                assert!(flag_array[[1, fine_ch_2_idx, 0]]);
            }
        }

        // Test that raw_chan_flags work
        // Fine channel 1 in coarse channel 0 should be flagged
        assert!(flag_array[[1, 1, 0]]);
        // Fine channel 3 in coarse channel 2 should be flagged
        assert!(flag_array[[1, 2 * num_fine_chans_per_coarse + 3, 0]]);

        // Test that non-flagged channels remain unflagged
        assert!(!flag_array[[1, 0, 0]]); // Fine ch 0 in coarse ch 0
        assert!(!flag_array[[1, 2 * num_fine_chans_per_coarse + 0, 0]]); // Fine ch 0 in coarse ch 2
    }
}
