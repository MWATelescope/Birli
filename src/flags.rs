//! Methods for manipulating flagmasks and flagging imagesets

use std::ops::Range;

use crate::{
    io::error::IOError,
    marlu::{
        mwalib::CorrelatorContext,
        ndarray::{Array, Array3, ArrayView, Dimension},
    },
    BirliError, FlagFileSet,
};
use cfg_if::cfg_if;
use derive_builder::Builder;
use itertools::izip;
use log::trace;

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
    /// Which mwalib antenna indices are flagged
    pub antenna_flags: Vec<bool>,
    /// Whether auto-correlations are flagged
    #[builder(default = "false")]
    pub autos: bool,
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
    /// let corr_ctx = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
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

        result
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

    /// Set flags from this context in an existing array.
    ///
    /// # Errors
    ///
    /// Can throw error if array is not the correct shape.
    pub fn set_flags(
        &self,
        flag_array: &mut Array3<bool>,
        timestep_range: &Range<usize>,
        coarse_chan_range: &Range<usize>,
        ant_pairs: &[(usize, usize)],
    ) -> Result<(), BirliError> {
        let timestep_flags = &self.timestep_flags[timestep_range.clone()];
        let coarse_chan_flags = &self.coarse_chan_flags[coarse_chan_range.clone()];
        let baseline_flags = self.get_baseline_flags(ant_pairs);

        let chan_flags: Vec<_> = coarse_chan_flags
            .iter()
            .flat_map(|coarse_chan_flag| {
                if *coarse_chan_flag {
                    vec![true; self.fine_chan_flags.len()]
                } else {
                    self.fine_chan_flags.clone()
                }
            })
            .collect();
        let shape = (timestep_range.len(), chan_flags.len(), ant_pairs.len());

        let flag_shape = flag_array.dim();
        if flag_shape.0 > shape.0 || flag_shape.1 > shape.1 || flag_shape.2 > shape.2 {
            return Err(BirliError::BadArrayShape {
                argument: "flag_array".to_string(),
                function: "FlagContext::set_flags".to_string(),
                expected: format!("dims less than {:?}", shape),
                received: format!("{:?}", flag_shape),
            });
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
/// - `baseline_jones_view` is [timestep][channel] for one baseline
/// - imageset is timesteps wide, and channels high
/// - jones matrics are always XX, YY, XY, YX
///
/// # Errors
///
/// TODO: this doesn't actually throw any errors?
#[cfg(feature = "aoflagger")]
pub fn jones_baseline_view_to_imageset(
    aoflagger: &CxxAOFlagger,
    baseline_jones_view: &ArrayView2<Jones<f32>>,
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
/// TODO: this doesn't actually throw any errors?
#[cfg(feature = "aoflagger")]
pub fn flag_baseline_view_to_flagmask(
    aoflagger: &CxxAOFlagger,
    baseline_flag_view: &ArrayView2<bool>,
) -> Result<UniquePtr<CxxFlagMask>, BirliError> {
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
    Ok(flag_mask)
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
/// use birli::{FlagContext, flag_jones_array_existing, write_flags,
///     mwalib::CorrelatorContext, cxx_aoflagger_new, VisSelection};
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
/// let corr_ctx = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
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
/// vis_sel
///     .read_mwalib(&corr_ctx, jones_array.view_mut(), flag_array.view_mut(), false)
///     .unwrap();
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

/// Write flags to disk, given an observation's [`marlu::mwalib::CorrelatorContext`], a vector of
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
/// use birli::{FlagContext, write_flags, mwalib::CorrelatorContext, VisSelection};
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
/// let corr_ctx = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
///
/// // Determine which timesteps and coarse channels we want to use
/// let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
///
/// // Prepare our flagmasks with known bad antennae
/// let mut flag_ctx = FlagContext::from_mwalib(&corr_ctx);
///
/// // Create a blank array to store flags and visibilities
/// let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
/// let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
/// flag_ctx.set_flags(
///     &mut flag_array,
///     &vis_sel.timestep_range,
///     &vis_sel.coarse_chan_range,
///     &vis_sel.get_ant_pairs(&corr_ctx.metafits_context)
/// );
/// let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
///
/// // read visibilities out of the gpubox files
/// vis_sel
///     .read_mwalib(&corr_ctx, jones_array.view_mut(), flag_array.view_mut(), false)
///     .unwrap();
///
/// // write the flags to disk as .mwaf
/// write_flags(&corr_ctx, &flag_array, flag_template.to_str().unwrap(), &vis_sel.coarse_chan_range).unwrap();
/// ```
///
/// # Errors
///
/// - Will error with [`IOError::FitsOpen`] if there are files already present at the paths specified in filename template.
/// - Will error with [`IOError::InvalidFlagFilenameTemplate`] if an invalid flag filename template is provided (wrong number of percents).

pub fn write_flags(
    corr_ctx: &CorrelatorContext,
    flag_array: &Array3<bool>,
    filename_template: &str,
    coarse_chan_range: &Range<usize>,
) -> Result<(), IOError> {
    trace!("start write_flags");

    let gpubox_ids = corr_ctx.coarse_chans[coarse_chan_range.clone()]
        .iter()
        .map(|chan| chan.gpubox_number)
        .collect::<Vec<_>>();

    trace!(
        "writing flags to template: {}, gpubox ids: {:?}",
        filename_template,
        gpubox_ids
    );

    let mut flag_file_set = FlagFileSet::new(filename_template, &gpubox_ids, corr_ctx.mwa_version)?;
    flag_file_set.write_flag_array(corr_ctx, flag_array, &gpubox_ids)?;

    trace!("end write_flags");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::write_flags;
    use glob::glob;
    use tempfile::tempdir;

    use crate::{
        marlu::selection::SelectionError::{NoCommonTimesteps, NoProvidedTimesteps},
        test_common::get_mwax_context,
        test_common::{
            get_mwa_ord_dodgy_context, get_mwa_ord_no_overlap_context,
            get_mwa_ord_no_timesteps_context,
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
            &corr_ctx,
            &flag_array,
            filename_template.to_str().unwrap(),
            &vis_sel.coarse_chan_range,
        )
        .unwrap();

        let flag_files = glob(tmp_dir.path().join("Flagfile*.mwaf").to_str().unwrap()).unwrap();

        assert_eq!(flag_files.count(), selected_gpuboxes.len());

        let mut flag_file_set = FlagFileSet::open(
            filename_template.to_str().unwrap(),
            &selected_gpuboxes,
            corr_ctx.mwa_version,
        )
        .unwrap();
        let chan_header_flags_raw = flag_file_set.read_chan_header_flags_raw().unwrap();
        let (chan1_header, chan1_flags_raw) =
            chan_header_flags_raw.get(&selected_gpuboxes[0]).unwrap();
        assert_eq!(chan1_header.gpubox_id, gpubox_ids[0]);
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;

        let num_baselines = chan1_header.num_ants * (chan1_header.num_ants + 1) / 2;
        assert_eq!(chan1_header.num_timesteps, corr_ctx.num_timesteps);
        assert_eq!(num_baselines, corr_ctx.metafits_context.num_baselines);
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
        for (timestep_idx, baseline_idx, fine_chan_idx, expected_flag) in tests {
            let row_idx = timestep_idx * num_baselines + baseline_idx;
            let offset = row_idx * fine_chans_per_coarse + fine_chan_idx;
            assert_eq!(
                chan1_flags_raw[offset], expected_flag,
                "with timestep {}, baseline {}, fine_chan {}, expected {} at row_idx {}, offset {}",
                timestep_idx, baseline_idx, fine_chan_idx, expected_flag, row_idx, offset
            );
        }
    }
}

/// Get the weight factor of an observation's `corr_ctx`.
///
/// This is a conceptfrom Cotter, and the legacy MWA correlator where the value
/// is a multiple of the frequency averaging factor (relative to 10kHz), and the
/// time averaging factor (relative to 1s).
pub fn get_weight_factor(corr_ctx: &CorrelatorContext) -> f64 {
    let integration_time_s = corr_ctx.metafits_context.corr_int_time_ms as f64 / 1000.0;
    let fine_chan_width_hz = corr_ctx.metafits_context.corr_fine_chan_width_hz as f64;
    fine_chan_width_hz * integration_time_s / 10000.0
}

/// Convert the given ndarray of boolean flags to an ndarray of float weights
pub fn flag_to_weight_array<D>(flag_array: &ArrayView<bool, D>, weight_factor: f64) -> Array<f32, D>
where
    D: Dimension,
{
    flag_array.map(|f| if *f { -weight_factor } else { weight_factor } as f32)
}

#[cfg(test)]
#[cfg(feature = "aoflagger")]
/// Tests which require the use of the aoflagger feature
mod tests_aoflagger {
    use marlu::{mwalib::CorrelatorContext, Complex, Jones};
    use ndarray::Array3;

    use crate::{
        flags::{flag_jones_array, flag_jones_array_existing, FlagContext},
        VisSelection,
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

        let corr_ctx = get_mwa_ord_context();
        let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

        let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        flag_ctx
            .set_flags(
                &mut flag_array,
                &vis_sel.timestep_range,
                &vis_sel.coarse_chan_range,
                &vis_sel.get_ant_pairs(&corr_ctx.metafits_context),
            )
            .unwrap();
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        vis_sel
            .read_mwalib(
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
            &jones_array,
            &mut flag_array,
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
}
