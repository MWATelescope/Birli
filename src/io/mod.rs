//! Input and Ouput data file format modules
//!
//! Most of this was blatently stolen (with permission) from [Chris Jordan](https://github.com/cjordan)

pub mod aocal;
pub mod error;
pub mod mwaf;

use std::{
    ops::Range,
    path::{Path, PathBuf},
};

use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};
use itertools::izip;
use log::{trace, warn};

use crate::{
    marlu::{
        constants::MWA_LAT_RAD,
        fitsio,
        hifitime::Duration,
        io::{ms::MeasurementSetWriter, uvfits::UvfitsWriter, VisWrite},
        mwalib,
        mwalib::{
            CorrelatorContext, MwalibError, _get_optional_fits_key, _open_hdu, fits_open_hdu,
            get_optional_fits_key,
        },
        rayon::prelude::*,
        Jones, LatLngHeight, MwaObsContext, ObsContext, RADec, SelectionError, VisContext,
        VisSelection, ENH,
    },
    ndarray::prelude::*,
};

use self::error::IOError;

/// Groups together parameters related to I/O
#[derive(Debug, Default, Clone)]
pub struct IOContext {
    // in
    /// The path to the .metafits input file
    pub metafits_in: PathBuf,
    /// A vector of gpufits .fits input paths
    pub gpufits_in: Vec<PathBuf>,
    /// Optional path to a .bin ao calibration solutions input file
    pub aocalsols_in: Option<PathBuf>,

    // out
    /// Optional .uvfits output path
    pub uvfits_out: Option<PathBuf>,
    /// Optional .ms measurement set output path
    pub ms_out: Option<PathBuf>,
    /// Optional .mwaf flag file path template (see `io::mwaf::FlagFileSet`)
    pub flag_template: Option<String>,
}

impl IOContext {
    // get the scale factor for the raw files
    // will be deprecated after <https://github.com/MWATelescope/mwalib/issues/85> is resolved
    fn get_raw_scale_factor(&self) -> f32 {
        let mut fptr = fitsio::FitsFile::open(&self.metafits_in).unwrap();
        let hdu0 = fits_open_hdu!(&mut fptr, 1).unwrap();
        let mut scale_factor: Option<f32> =
            get_optional_fits_key!(&mut fptr, &hdu0, "RAWSCALE").unwrap();
        drop(fptr);
        for raw_path in &self.gpufits_in {
            let mut fptr = fitsio::FitsFile::open(raw_path).unwrap();
            let hdu1 = fits_open_hdu!(&mut fptr, 1).unwrap();
            for key in ["SCALEFAC", "BSCALE"] {
                let this_scale_factor: Option<f32> =
                    get_optional_fits_key!(&mut fptr, &hdu1, key).unwrap();
                match (scale_factor, this_scale_factor) {
                    (Some(sf), Some(this_sf)) => {
                        assert!(
                            ((sf - this_sf).abs() < f32::EPSILON),
                            "Different scale factors found in raw files: {} and {}",
                            sf,
                            this_sf
                        );
                    }
                    (None, Some(this_sf)) => {
                        scale_factor = Some(this_sf);
                    }
                    _ => {}
                }
            }
        }
        scale_factor.unwrap_or(1.0)
    }

    /// Get the `mwalib::CorrelatorContext` from metafits and gpufits
    ///
    /// # Errors
    ///
    /// see `mwalib::CorrelatorContext::new`
    pub fn get_corr_ctx(&self) -> Result<CorrelatorContext, MwalibError> {
        CorrelatorContext::new(&self.metafits_in, &self.gpufits_in).map(|mut corr_ctx| {
            corr_ctx.metafits_context.corr_raw_scale_factor = self.get_raw_scale_factor();
            corr_ctx
        })
    }

    // TODO: pub fn validate_params(&self), checks permissions
}

/// The container has visibilities which can be read by passing in a mwalib
/// [`CorrelatorContext`] and the range of values to read.
pub trait ReadableVis: Sync + Send {
    /// Read the visibilities and weights for the selected timesteps, coarse
    /// channels and baselines into the provided arrays.
    ///
    /// # Errors
    ///
    /// Can throw [`IOError`] if there is an issue reading.
    ///
    /// TODO: reduce number of arguments.
    #[allow(clippy::too_many_arguments)]
    fn read_vis_mwalib(
        &self,
        jones_array: ArrayViewMut3<Jones<f32>>,
        weight_array: ArrayViewMut3<f32>,
        corr_ctx: &CorrelatorContext,
        timestep_range: &Range<usize>,
        coarse_chan_range: &Range<usize>,
        baseline_idxs: &[usize],
    ) -> Result<(), IOError>;
}

/// Read the visibilities for this selection into the jones array using mwalib,
/// flag visiblities if they are not provided.
///
/// # Errors
///
/// Can raise [`SelectionError::BadArrayShape`] if `jones_array` or `flag_array` does not match the
/// expected shape of this selection.
///
/// # Examples
///
/// ```rust
/// use marlu::{mwalib::CorrelatorContext, VisSelection};
/// use birli::io::read_mwalib;
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
/// // Determine which timesteps and coarse channels we want to use
/// let img_timestep_idxs = &corr_ctx.common_timestep_indices;
/// let good_timestep_idxs = &corr_ctx.common_good_timestep_indices;
///
/// let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
/// vis_sel.timestep_range =
///     *img_timestep_idxs.first().unwrap()..(*img_timestep_idxs.last().unwrap() + 1);
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
/// let dims_common = jones_array.dim();
///
/// // now try only with good timesteps
/// vis_sel.timestep_range =
///     *good_timestep_idxs.first().unwrap()..(*good_timestep_idxs.last().unwrap() + 1);
///
/// // read visibilities out of the gpubox files
/// let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
/// let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
/// read_mwalib(&vis_sel, &corr_ctx, jones_array.view_mut(), flag_array.view_mut(), false)
///     .unwrap();
///
/// let dims_good = jones_array.dim();
///
/// // different selections have different sized arrays.
/// assert_ne!(dims_common, dims_good);
/// ```
pub fn read_mwalib(
    vis_sel: &VisSelection,
    corr_ctx: &CorrelatorContext,
    mut jones_array: ArrayViewMut3<Jones<f32>>,
    mut flag_array: ArrayViewMut3<bool>,
    draw_progress: bool,
) -> Result<(), SelectionError> {
    let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
    let shape = vis_sel.get_shape(fine_chans_per_coarse);
    let (num_timesteps, _, _) = shape;
    let num_coarse_chans = vis_sel.coarse_chan_range.len();
    let max_bl_idx = corr_ctx.metafits_context.baselines.len();

    // check output array dimensions
    if jones_array.dim() != shape {
        return Err(SelectionError::BadArrayShape {
            argument: "jones_array".to_string(),
            function: "VisSelection::read_mwalib".to_string(),
            expected: format!("{shape:?}"),
            received: format!("{:?}", jones_array.dim()),
        });
    };

    if flag_array.dim() != shape {
        return Err(SelectionError::BadArrayShape {
            argument: "flag_array".to_string(),
            function: "VisSelection::read_mwalib".to_string(),
            expected: format!("{shape:?}"),
            received: format!("{:?}", flag_array.dim()),
        });
    };

    // check all selected baseline idxs are < max_bl_idx
    if vis_sel.baseline_idxs.iter().any(|&idx| idx >= max_bl_idx) {
        return Err(SelectionError::BadBaselineIdx {
            function: "VisSelection::read_mwalib".to_string(),
            expected: format!(" < {max_bl_idx}"),
            received: format!("{:?}", vis_sel.baseline_idxs.clone()),
        });
    }

    // since we are using read_by_baseline_into_buffer, the visibilities are read in order:
    // baseline,frequency,pol,r,i

    // compiler optimization
    let floats_per_chan = 8;
    assert_eq!(
        corr_ctx.metafits_context.num_visibility_pols * 2,
        floats_per_chan
    );

    let floats_per_baseline = floats_per_chan * fine_chans_per_coarse;
    let floats_per_hdu = floats_per_baseline * corr_ctx.metafits_context.num_baselines;

    // Progress bar draw target
    let draw_target = if draw_progress {
        ProgressDrawTarget::stderr()
    } else {
        ProgressDrawTarget::hidden()
    };
    // a progress bar containing the progress bars associated with this method
    let multi_progress = MultiProgress::with_draw_target(draw_target);
    // a vector of progress bars for the visibility reading progress of each channel.
    let read_progress: Vec<ProgressBar> = vis_sel
        .coarse_chan_range
        .clone()
        .map(|mwalib_coarse_chan_idx| {
            let channel_progress = multi_progress.add(
                ProgressBar::new(num_timesteps as _)
                    .with_style(
                        ProgressStyle::default_bar()
                            .template("{msg:16}: [{wide_bar:.blue}] {pos:4}/{len:4}")
                            .unwrap()
                            .progress_chars("=> "),
                    )
                    .with_position(0)
                    .with_message(format!("coarse_chan {mwalib_coarse_chan_idx:03}")),
            );
            channel_progress.set_position(0);
            channel_progress
        })
        .collect();
    // The total reading progress bar.
    let total_progress = multi_progress.add(
            ProgressBar::new((num_timesteps * num_coarse_chans) as _)
                .with_style(
                    ProgressStyle::default_bar()
                        .template(
                            "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
                        )
                        .unwrap()
                        .progress_chars("=> "),
                )
                .with_position(0)
                .with_message("loading hdus"),
        );

    // Load HDUs from each coarse channel. arrays: [timestep][chan][baseline]
    jones_array
        .axis_chunks_iter_mut(Axis(1), fine_chans_per_coarse)
        .into_par_iter()
        .zip(flag_array.axis_chunks_iter_mut(Axis(1), fine_chans_per_coarse))
        .zip(vis_sel.coarse_chan_range.clone())
        .zip(read_progress)
        .try_for_each(
            |(((mut jones_array, mut flag_array), coarse_chan_idx), progress)| {
                progress.set_position(0);

                // buffer: [baseline][chan][pol][complex]
                let mut hdu_buffer: Vec<f32> = vec![0.0; floats_per_hdu];

                // arrays: [chan][baseline]
                for (mut jones_array, mut flag_array, timestep_idx) in izip!(
                    jones_array.outer_iter_mut(),
                    flag_array.outer_iter_mut(),
                    vis_sel.timestep_range.clone(),
                ) {
                    match corr_ctx.read_by_baseline_into_buffer(
                        timestep_idx,
                        coarse_chan_idx,
                        hdu_buffer.as_mut_slice(),
                    ) {
                        Ok(()) => {
                            // arrays: [chan]
                            for (mut jones_array, baseline_idx) in izip!(
                                jones_array.axis_iter_mut(Axis(1)),
                                vis_sel.baseline_idxs.iter()
                            ) {
                                // buffer: [chan][pol][complex]
                                let hdu_baseline_chunk = &hdu_buffer
                                    [baseline_idx * floats_per_baseline..][..floats_per_baseline];
                                for (jones, hdu_chan_chunk) in izip!(
                                    jones_array.iter_mut(),
                                    hdu_baseline_chunk.chunks_exact(floats_per_chan)
                                ) {
                                    *jones = Jones::from([
                                        hdu_chan_chunk[0],
                                        hdu_chan_chunk[1],
                                        hdu_chan_chunk[2],
                                        hdu_chan_chunk[3],
                                        hdu_chan_chunk[4],
                                        hdu_chan_chunk[5],
                                        hdu_chan_chunk[6],
                                        hdu_chan_chunk[7],
                                    ]);
                                }
                            }
                        }
                        Err(mwalib::GpuboxError::NoDataForTimeStepCoarseChannel { .. }) => {
                            warn!(
                                "Flagging missing HDU @ ts={}, cc={}",
                                timestep_idx, coarse_chan_idx
                            );
                            flag_array.fill(true);
                        }
                        Err(e) => return Err(e),
                    }

                    progress.inc(1);
                    total_progress.inc(1);
                }
                progress.finish();
                Ok(())
            },
        )?;

    // We're done!
    total_progress.finish();

    Ok(())
}

/// Write the given ndarrays of flags and [`Jones`] matrix visibilities to a
/// uvfits file.
///
/// mwalib timestep, coarse channel and baseline indices are needed to map between
/// indices in the arrays and indices according to mwalib, which are not the same.
///
/// # Examples
///
/// ```rust
/// use tempfile::tempdir;
/// use birli::{
///     FlagContext,
///     marlu::mwalib::CorrelatorContext,
///     get_weight_factor,
///     flag_to_weight_array,
///     VisSelection, io::{read_mwalib, write_uvfits}
/// };
///
/// // define our input files
/// let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
/// let gpufits_paths = vec![
///     "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits",
/// ];
///
/// // define a temporary directory for output files
/// let tmp_dir = tempdir().unwrap();
///
/// // define our output flag file template
/// let uvfits_out = tmp_dir.path().join("synthetic.uvfits");
///
/// // Create an mwalib::CorrelatorContext for accessing visibilities.
/// let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();
///
/// // Determine which timesteps and coarse channels we want to use
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
/// // write the visibilities to disk as .uvfits
/// let num_pols = corr_ctx.metafits_context.num_visibility_pols;
/// let weight_factor = get_weight_factor(&corr_ctx);
/// let weight_array = flag_to_weight_array(flag_array.view(), weight_factor);
///
/// write_uvfits(
///     uvfits_out.as_path(),
///     &corr_ctx,
///     jones_array.view(),
///     weight_array.view(),
///     &vis_sel.timestep_range,
///     &vis_sel.coarse_chan_range,
///     &vis_sel.baseline_idxs,
///     None,
///     None,
///     1,
///     1,
/// )
/// .unwrap();
/// ```
/// # Errors
///
/// See: [`UvfitsWriter`]
///
/// TODO: reduce number of arguments.
#[allow(clippy::too_many_arguments)]
pub fn write_uvfits<T: AsRef<Path>>(
    path: T,
    corr_ctx: &CorrelatorContext,
    jones_array: ArrayView3<Jones<f32>>,
    weight_array: ArrayView3<f32>,
    timestep_range: &Range<usize>,
    coarse_chan_range: &Range<usize>,
    baseline_idxs: &[usize],
    array_pos: Option<LatLngHeight>,
    phase_centre: Option<RADec>,
    avg_time: usize,
    avg_freq: usize,
) -> Result<(), IOError> {
    trace!("start write_uvfits to {:?}", path.as_ref());

    let vis_ctx = VisContext::from_mwalib(
        corr_ctx,
        timestep_range,
        coarse_chan_range,
        baseline_idxs,
        avg_time,
        avg_freq,
    );

    let mut obs_ctx = ObsContext::from_mwalib(&corr_ctx.metafits_context);
    if let Some(phase_centre) = phase_centre {
        obs_ctx.phase_centre = phase_centre;
    }

    let array_latitude = match array_pos {
        Some(array_pos) => {
            obs_ctx.array_pos = array_pos;
            array_pos.latitude_rad
        }
        None => {
            // The writer will assume MWA, so we do here too.
            MWA_LAT_RAD
        }
    };
    let (antenna_names, antenna_positions) = corr_ctx
        .metafits_context
        .antennas
        .iter()
        .map(|a| {
            let enh = ENH {
                e: a.east_m,
                n: a.north_m,
                h: a.height_m,
            };
            let (s_lat, c_lat) = array_latitude.sin_cos();
            let xyz = enh.to_xyz_inner(s_lat, c_lat);
            (a.tile_name.clone(), xyz)
        })
        .unzip();

    let mut uvfits_writer = UvfitsWriter::from_marlu(
        path,
        &vis_ctx,
        obs_ctx.array_pos,
        obs_ctx.phase_centre,
        Duration::from_seconds(corr_ctx.metafits_context.dut1.unwrap_or(0.0)),
        obs_ctx.name.as_deref(),
        antenna_names,
        antenna_positions,
        true,
        None,
    )?;

    uvfits_writer.write_vis(jones_array, weight_array, &vis_ctx)?;

    uvfits_writer.finalise()?;

    trace!("end write_uvfits");

    Ok(())
}

/// Write the given ndarrays of flags and [`Jones`] matrix visibilities to a
/// measurement set.
///
/// mwalib timestep, coarse channel and baseline indices are needed to map between
/// indices in the arrays and indices according to mwalib, which are not the same.
///
/// # Examples
///
/// ```rust
/// use tempfile::tempdir;
/// use birli::{
///     VisSelection,
///     marlu::mwalib::CorrelatorContext,
///     get_weight_factor,
///     flag_to_weight_array,
///     FlagContext, io::{read_mwalib, write_ms}
/// };
///
/// // define our input files
/// let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
/// let gpufits_paths = vec![
///     "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits",
/// ];
///
/// // define a temporary directory for output files
/// let tmp_dir = tempdir().unwrap();
///
/// // define our output flag file template
/// let ms_out = tmp_dir.path().join("synthetic.ms");
///
/// // Create an mwalib::CorrelatorContext for accessing visibilities.
/// let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();
///
/// // Determine which timesteps and coarse channels we want to use
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
/// // write the visibilities to disk as .ms
///
/// let num_pols = corr_ctx.metafits_context.num_visibility_pols;
/// let weight_factor = get_weight_factor(&corr_ctx);
/// let weight_array = flag_to_weight_array(flag_array.view(), weight_factor);
/// // time and frequency averaging
/// let (avg_time, avg_freq) = (1, 1);
/// write_ms(
///     ms_out.as_path(),
///     &corr_ctx,
///     jones_array.view(),
///     weight_array.view(),
///     &vis_sel.timestep_range,
///     &vis_sel.coarse_chan_range,
///     &vis_sel.baseline_idxs,
///     None,
///     None,
///     avg_time,
///     avg_freq,
/// )
/// .unwrap();
/// ```
/// # Errors
///
/// See: [`UvfitsWriter`]
///
/// TODO: reduce number of arguments.
#[allow(clippy::too_many_arguments)]
pub fn write_ms<T: AsRef<Path>>(
    path: T,
    corr_ctx: &CorrelatorContext,
    jones_array: ArrayView3<Jones<f32>>,
    weight_array: ArrayView3<f32>,
    timestep_range: &Range<usize>,
    coarse_chan_range: &Range<usize>,
    baseline_idxs: &[usize],
    array_pos: Option<LatLngHeight>,
    phase_centre: Option<RADec>,
    avg_time: usize,
    avg_freq: usize,
) -> Result<(), IOError> {
    trace!("start write_ms to {:?}", path.as_ref());

    let vis_ctx = VisContext::from_mwalib(
        corr_ctx,
        timestep_range,
        coarse_chan_range,
        baseline_idxs,
        avg_time,
        avg_freq,
    );

    let mut obs_ctx = ObsContext::from_mwalib(&corr_ctx.metafits_context);
    if let Some(phase_centre) = phase_centre {
        obs_ctx.phase_centre = phase_centre;
    }
    if let Some(array_pos) = array_pos {
        obs_ctx.array_pos = array_pos;
    }
    let mwa_ctx = MwaObsContext::from_mwalib(&corr_ctx.metafits_context);

    let mut ms_writer = MeasurementSetWriter::new(
        path,
        obs_ctx.phase_centre,
        obs_ctx.array_pos,
        obs_ctx.ant_positions_geodetic().collect(),
        Duration::from_seconds(corr_ctx.metafits_context.dut1.unwrap_or(0.0)),
        true,
    );

    ms_writer
        .initialize_mwa(&vis_ctx, &obs_ctx, &mwa_ctx, None, coarse_chan_range)
        .unwrap();

    ms_writer
        .write_vis(jones_array.view(), weight_array.view(), &vis_ctx)
        .unwrap();

    trace!("end write_ms");

    Ok(())
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use marlu::{Complex, Jones};

    use crate::{compare_jones, test_common::get_mwax_context, VisSelection};

    use super::read_mwalib;

    // test read_mwalib with bad vis_sel.baseline_idxs
    #[test]
    fn test_read_bad_baseline_sel() {
        let corr_ctx = get_mwax_context();

        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        vis_sel.baseline_idxs = vec![5];

        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        assert!(read_mwalib(
            &vis_sel,
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .is_err());
    }

    // test read_mwalib with custom vis_sel.baseline_idxs
    #[test]
    fn test_read_baseline_sel() {
        let corr_ctx = get_mwax_context();

        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        vis_sel.baseline_idxs = vec![1];

        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        assert_eq!(flag_array.shape(), &[4, 4, 1]);
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        assert_eq!(jones_array.shape(), &[4, 4, 1]);
        let weight_array = vis_sel.allocate_weights(fine_chans_per_coarse).unwrap();
        assert_eq!(weight_array.shape(), &[4, 4, 1]);
        read_mwalib(
            &vis_sel,
            &corr_ctx,
            jones_array.view_mut(),
            flag_array.view_mut(),
            false,
        )
        .unwrap();

        // ts 0, chan 0 (cc 0, fc 0), baseline 1
        let viz_0_0_1 = jones_array[(0, 0, 0)];
        compare_jones!(
            viz_0_0_1,
            Jones::from([
                Complex::new(0x410010 as f32, 0x410011 as f32),
                Complex::new(0x410012 as f32, 0x410013 as f32),
                Complex::new(0x410014 as f32, 0x410015 as f32),
                Complex::new(0x410016 as f32, 0x410017 as f32),
            ])
        );

        // ts 3, chan 3 (cc 1, fc 1), baseline 1
        let viz_3_3_1 = jones_array[(3, 3, 0)];
        compare_jones!(
            viz_3_3_1,
            Jones::from([
                Complex::new(0x410718 as f32, 0x410719 as f32),
                Complex::new(0x41071a as f32, 0x41071b as f32),
                Complex::new(0x41071c as f32, 0x41071d as f32),
                Complex::new(0x41071e as f32, 0x41071f as f32),
            ])
        );
    }
}

#[cfg(test)]
#[cfg(feature = "aoflagger")]
/// Tests which require the use of the aoflagger feature
mod tests_aoflagger {
    use crate::{
        flags::{flag_jones_array_existing, flag_to_weight_array, get_weight_factor},
        io::{read_mwalib, write_uvfits},
        FlagContext, VisSelection,
    };
    use aoflagger_sys::cxx_aoflagger_new;
    use fitsio::errors::check_status as fits_check_status;
    use float_cmp::{approx_eq, F32Margin};
    use itertools::izip;
    use marlu::{
        fitsio, fitsio_sys,
        mwalib::{
            _get_required_fits_key, _open_fits, _open_hdu, fits_open, fits_open_hdu,
            get_required_fits_key, CorrelatorContext,
        },
    };
    use tempfile::tempdir;

    #[test]
    fn write_uvfits_flags() {
        // define our input files
        let metafits_path = "tests/data/1196175296_mwa_ord/1196175296.metafits";
        let gpufits_paths = vec![
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits",
            "tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits",
        ];

        // define a temporary directory for output files
        let tmp_dir = tempdir().unwrap();

        // define our output uvfits path
        let uvfits_out = tmp_dir.path().join("1297526432.uvfits");

        // Create an mwalib::CorrelatorContext for accessing visibilities.
        let corr_ctx = CorrelatorContext::new(metafits_path, &gpufits_paths).unwrap();

        // create a CxxAOFlagger object to perform AOFlagger operations
        let aoflagger = unsafe { cxx_aoflagger_new() };

        // Determine which timesteps and coarse channels we want to use
        let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

        // Prepare our flagmasks with known bad antennae
        let mut flag_ctx = FlagContext::from_mwalib(&corr_ctx);
        flag_ctx.flag_dc = false;
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

        // use the default strategy file location for MWA
        let strategy_filename = &aoflagger.FindStrategyFileMWA();

        // run the strategy on the imagesets, and get the resulting flagmasks for each baseline
        flag_jones_array_existing(
            &aoflagger,
            strategy_filename,
            jones_array.view(),
            flag_array.view_mut(),
            true,
            false,
        );

        let weight_factor = get_weight_factor(&corr_ctx);
        let weight_array = flag_to_weight_array(flag_array.view(), weight_factor);

        // write the visibilities to disk as .uvfits

        write_uvfits(
            uvfits_out.as_path(),
            &corr_ctx,
            jones_array.view(),
            weight_array.view(),
            &vis_sel.timestep_range,
            &vis_sel.coarse_chan_range,
            &vis_sel.baseline_idxs,
            None,
            None,
            1,
            1,
        )
        .unwrap();

        // Test the values have been correctly populated.
        let mut fptr = fits_open!(&uvfits_out).unwrap();
        let vis_hdu = fits_open_hdu!(&mut fptr, 0).unwrap();

        let pcount: usize = get_required_fits_key!(&mut fptr, &vis_hdu, "PCOUNT").unwrap();
        let floats_per_pol: usize = get_required_fits_key!(&mut fptr, &vis_hdu, "NAXIS2").unwrap();
        let num_pols: usize = get_required_fits_key!(&mut fptr, &vis_hdu, "NAXIS3").unwrap();
        let num_fine_freq_chans: usize =
            get_required_fits_key!(&mut fptr, &vis_hdu, "NAXIS4").unwrap();

        let vis_len = num_fine_freq_chans * num_pols * floats_per_pol;
        assert_eq!(vis_len, 48);

        let expected = [
            (
                1, // first non-autocorrelated baseline
                [
                    -0.00000077589283,
                    0.0000005552067,
                    -0.00000009305131,
                    258.0,
                    -0.37870082,
                ],
                [
                    0x10c59e, 0x10c59f, 32, // XX 0
                    0x10bea6, 0x10bea7, 32, // YY 0
                    0x10c58e, 0x10c58f, 32, // XY 0
                    0x10beb6, 0x10beb7, 32, // YX 0
                    0x11c79e, 0x11c79f, 32, // XX 1
                    0x11c0a6, 0x11c0a7, 32, // YY 1
                    0x11c78e, 0x11c78f, 32, // XY 1
                    0x11c0b6, 0x11c0b7, 32, // YX 1
                    0x00c59e, 0x00c59f, 32, // XX 2
                    0x00bea6, 0x00bea7, 32, // YY 2
                    0x00c58e, 0x00c58f, 32, // XY 2
                    0x00beb6, 0x00beb7, 32, // YX 2
                    0x01c79e, 0x01c79f, 32, // XX 3
                    0x01c0a6, 0x01c0a7, 32, // YY 3
                    0x01c78e, 0x01c78f, 32, // XY 3
                    0x01c0b6, 0x01c0b7, 32, // YX 3
                ],
            ),
            (
                11, // where we start to see flagged antennae
                [
                    0.0000011107809,
                    0.00000055093767,
                    -0.000000102761746,
                    268.0,
                    -0.37870082,
                ],
                [
                    0x10c25e, 0x10c25f, -32, // XX 0
                    0x10bb66, 0x10bb67, -32, // YY 0
                    0x10c24e, 0x10c24f, -32, // XY 0
                    0x10bb76, 0x10bb77, -32, // YX 0
                    0x11c45e, 0x11c45f, -32, // XX 1
                    0x11bd66, 0x11bd67, -32, // YY 1
                    0x11c44e, 0x11c44f, -32, // XY 1
                    0x11bd76, 0x11bd77, -32, // YX 1
                    0x00c25e, 0x00c25f, -32, // XX 2
                    0x00bb66, 0x00bb67, -32, // YY 2
                    0x00c24e, 0x00c24f, -32, // XY 2
                    0x00bb76, 0x00bb77, -32, // YX 2
                    0x01c45e, 0x01c45f, -32, // XX 3
                    0x01bd66, 0x01bd67, -32, // YY 3
                    0x01c44e, 0x01c44f, -32, // XY 3
                    0x01bd76, 0x01bd77, -32, // YX 3
                ],
            ),
        ];

        let mut status = 0;
        let mut obs_vis: Vec<f32> = vec![0.0; vis_len];
        let mut obs_group_params: Vec<f32> = vec![0.0; pcount];

        for (row_idx, exp_group_params, exp_vis) in &expected {
            unsafe {
                // ffggpe = fits_read_grppar_flt
                fitsio_sys::ffggpe(
                    fptr.as_raw(),                 /* I - FITS file pointer                       */
                    1 + *row_idx as i64,           /* I - group to read (1 = 1st group)           */
                    1,                             /* I - first vector element to read (1 = 1st)  */
                    pcount as i64,                 /* I - number of values to read                */
                    obs_group_params.as_mut_ptr(), /* O - array of values that are returned       */
                    &mut status,                   /* IO - error status                           */
                );
            }
            fits_check_status(status).unwrap();

            for (param_idx, (obs_group_param, exp_group_param)) in
                izip!(&obs_group_params, exp_group_params).enumerate()
            {
                assert!(
                    approx_eq!(
                        f32,
                        *obs_group_param,
                        *exp_group_param,
                        F32Margin::default()
                    ),
                    "cells don't match in param {param_idx}, row {row_idx}. {obs_group_params:?} != {exp_group_params:?}"
                );
            }

            unsafe {
                // ffgpve = fits_read_img_flt
                fitsio_sys::ffgpve(
                    fptr.as_raw(),        /* I - FITS file pointer                       */
                    1 + *row_idx as i64,  /* I - group to read (1 = 1st group)           */
                    1,                    /* I - first vector element to read (1 = 1st)  */
                    obs_vis.len() as i64, /* I - number of values to read                */
                    0.0,                  /* I - value for undefined pixels              */
                    obs_vis.as_mut_ptr(), /* O - array of values that are returned       */
                    &mut 0,               /* O - set to 1 if any values are null; else 0 */
                    &mut status,          /* IO - error status                           */
                );
            };
            fits_check_status(status).unwrap();

            for (vis_idx, (obs_val, exp_val)) in izip!(&obs_vis, exp_vis).enumerate() {
                assert!(
                    approx_eq!(f32, *obs_val, *exp_val as f32, F32Margin::default()),
                    "cells don't match in row {}, vis index {}. observed: {:?} != expected: {:?}",
                    row_idx,
                    vis_idx,
                    &obs_vis,
                    &exp_vis
                );
            }
        }
    }
}
