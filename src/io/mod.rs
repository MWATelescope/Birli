//! Input and Ouput data file format modules
//!
//! Most of this was blatently stolen (with permission) from [Chris Jordan](https://github.com/cjordan)

pub mod aocal;
pub mod error;
pub mod mwaf;

use std::{ops::Range, path::Path};

use log::trace;
use marlu::{MwaObsContext, ObsContext, VisContext};

use crate::{
    marlu::{
        io::{ms::MeasurementSetWriter, uvfits::UvfitsWriter, VisWritable},
        mwalib::{CorrelatorContext, MwalibError},
        Jones, LatLngHeight, RADec,
    },
    ndarray::{ArrayView3, ArrayViewMut3},
};

use self::error::IOError;

/// Groups together parameters related to I/O
#[derive(Debug, Default)]
pub struct IOContext {
    // in
    /// The path to the .metafits input file
    pub metafits_in: String,
    /// A vector of gpufits .fits input paths
    pub gpufits_in: Vec<String>,
    /// Optional path to a .bin ao calibration solutions input file
    pub aocalsols_in: Option<String>,

    // out
    /// Optional .uvfits output path
    pub uvfits_out: Option<String>,
    /// Optional .ms measurement set output path
    pub ms_out: Option<String>,
    /// Optional .mwaf flag file path template (see `io::mwaf::FlagFileSet`)
    pub flag_template: Option<String>,
}

impl IOContext {
    /// Get the `mwalib::CorrelatorContext` from metafits and gpufits
    ///
    /// # Errors
    ///
    /// see `mwalib::CorrelatorContext::new`
    pub fn get_corr_ctx(&self) -> Result<CorrelatorContext, MwalibError> {
        CorrelatorContext::new(&self.metafits_in, &self.gpufits_in)
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
///     write_uvfits,
///     marlu::mwalib::CorrelatorContext,
///     get_weight_factor,
///     flag_to_weight_array,
///     VisSelection
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
/// let corr_ctx = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
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
/// vis_sel
///     .read_mwalib(&corr_ctx, jones_array.view_mut(), flag_array.view_mut(), false)
///     .unwrap();
///
/// // write the visibilities to disk as .uvfits
/// let num_pols = corr_ctx.metafits_context.num_visibility_pols;
/// let weight_factor = get_weight_factor(&corr_ctx);
/// let weight_array = flag_to_weight_array(&flag_array.view(), weight_factor);
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
///     false,
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
    draw_progress: bool,
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
    if let Some(array_pos) = array_pos {
        obs_ctx.array_pos = array_pos;
    }

    let mut uvfits_writer = UvfitsWriter::from_marlu(
        path,
        &vis_ctx,
        Some(obs_ctx.array_pos),
        obs_ctx.phase_centre,
        obs_ctx.name.clone(),
    )?;

    let ant_positions_geodetic: Vec<_> = obs_ctx.ant_positions_geodetic().collect();

    uvfits_writer.write_vis_marlu(
        jones_array,
        weight_array,
        &vis_ctx,
        &ant_positions_geodetic,
        draw_progress,
    )?;

    uvfits_writer.write_ants_from_mwalib(&corr_ctx.metafits_context)?;

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
///     write_ms,
///     marlu::mwalib::CorrelatorContext,
///     get_weight_factor,
///     flag_to_weight_array,
///     FlagContext,
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
/// let corr_ctx = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
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
/// vis_sel
///     .read_mwalib(&corr_ctx, jones_array.view_mut(), flag_array.view_mut(), false)
///     .unwrap();
///
/// // write the visibilities to disk as .ms
///
/// let num_pols = corr_ctx.metafits_context.num_visibility_pols;
/// let weight_factor = get_weight_factor(&corr_ctx);
/// let weight_array = flag_to_weight_array(&flag_array.view(), weight_factor);
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
///     false,
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
    draw_progress: bool,
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

    let mut ms_writer =
        MeasurementSetWriter::new(path, obs_ctx.phase_centre, Some(obs_ctx.array_pos));

    ms_writer
        .initialize_mwa(&vis_ctx, &obs_ctx, &mwa_ctx, coarse_chan_range)
        .unwrap();

    let ant_positions_geodetic: Vec<_> = obs_ctx.ant_positions_geodetic().collect();

    ms_writer
        .write_vis_marlu(
            jones_array.view(),
            weight_array.view(),
            &vis_ctx,
            &ant_positions_geodetic,
            draw_progress,
        )
        .unwrap();

    trace!("end write_ms");

    Ok(())
}

#[cfg(test)]
#[cfg(feature = "aoflagger")]
/// Tests which require the use of the aoflagger feature
mod tests_aoflagger {
    use crate::{
        flags::{flag_jones_array_existing, flag_to_weight_array, get_weight_factor},
        write_uvfits, FlagContext, VisSelection,
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
        let corr_ctx = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();

        // create a CxxAOFlagger object to perform AOFlagger operations
        let aoflagger = unsafe { cxx_aoflagger_new() };

        // Determine which timesteps and coarse channels we want to use
        let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

        // Prepare our flagmasks with known bad antennae
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

        // use the default strategy file location for MWA
        let strategy_filename = &aoflagger.FindStrategyFileMWA();

        // run the strategy on the imagesets, and get the resulting flagmasks for each baseline
        flag_jones_array_existing(
            &aoflagger,
            strategy_filename,
            &jones_array,
            &mut flag_array,
            true,
            false,
        );

        let weight_factor = get_weight_factor(&corr_ctx);
        let weight_array = flag_to_weight_array(&flag_array.view(), weight_factor);

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
            false,
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
                    "cells don't match in param {}, row {}. {:?} != {:?}",
                    param_idx,
                    row_idx,
                    obs_group_params,
                    exp_group_params
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
