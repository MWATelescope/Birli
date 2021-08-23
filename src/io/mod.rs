//! Input and Ouput data file format modules
//!
//! Most of this was blatently stolen (with permission) from [Chris Jordan](https://github.com/cjordan)

pub mod error;
pub mod mwaf;
pub mod uvfits;

use std::path::Path;

use cxx::UniquePtr;
use log::trace;
use mwa_rust_core::{mwalib, LatLngHeight};
use mwalib::CorrelatorContext;
use uvfits::UvfitsWriter;

use crate::{CxxFlagMask, CxxImageSet};

use self::error::UvfitsWriteError;

/// Write the given [`BirliContext`] to a uvfits file.
///
/// # Errors
///
/// See: [`UvfitsWriter`]
///
/// TODO: replace all these args with birli_context
#[allow(clippy::too_many_arguments)]
pub fn write_uvfits<'a>(
    filename: &'a Path,
    context: &CorrelatorContext,
    baseline_idxs: &[usize],
    baseline_imgsets: &[UniquePtr<CxxImageSet>],
    baseline_flagmasks: &[UniquePtr<CxxFlagMask>],
    img_timestep_idxs: &[usize],
    img_coarse_chan_idxs: &[usize],
    array_pos: Option<LatLngHeight>,
) -> Result<(), UvfitsWriteError> {
    trace!("start write_uvfits");

    let mut uvfits_writer = UvfitsWriter::from_mwalib(
        filename,
        context,
        img_timestep_idxs,
        img_coarse_chan_idxs,
        baseline_idxs,
        array_pos,
    )?;

    let mut fits_file = uvfits_writer.open().unwrap();

    uvfits_writer.write_baseline_imgset_flagmasks(
        &mut fits_file,
        context,
        baseline_idxs,
        baseline_imgsets,
        baseline_flagmasks,
        img_timestep_idxs,
        img_coarse_chan_idxs,
    )?;

    uvfits_writer.write_ants_from_mwalib(&context.metafits_context)?;

    trace!("end write_uvfits");

    Ok(())
}

#[cfg(test)]
mod tests {

    use fitsio::errors::check_status as fits_check_status;
    use float_cmp::{approx_eq, F32Margin};
    use itertools::izip;
    use mwa_rust_core::{fitsio, fitsio_sys, mwalib};
    use mwalib::{
        CorrelatorContext, _get_required_fits_key, _open_fits, _open_hdu, fits_open, fits_open_hdu,
        get_required_fits_key,
    };
    use tempfile::tempdir;

    use crate::{
        context_to_baseline_imgsets, cxx_aoflagger_new, flag_imgsets_existing, get_antenna_flags,
        get_flaggable_timesteps, init_baseline_flagmasks,
    };

    use super::write_uvfits;

    #[test]
    fn write_uvfits_mwa_ord() {
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

        // define our output flag file template
        let uvfits_out = tmp_dir.path().join("1297526432.uvfits");

        // Create an mwalib::CorrelatorContext for accessing visibilities.
        let context = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();

        // create a CxxAOFlagger object to perform AOFlagger operations
        let aoflagger = unsafe { cxx_aoflagger_new() };

        // Determine which timesteps and coarse channels we want to use
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();

        // Prepare our flagmasks with known bad antennae
        let mut baseline_flagmasks = init_baseline_flagmasks(
            &aoflagger,
            &context,
            img_coarse_chan_idxs,
            &img_timestep_idxs,
            Some(get_antenna_flags(&context)),
        );

        // generate imagesets for each baseline in the format required by aoflagger
        let baseline_imgsets = context_to_baseline_imgsets(
            &aoflagger,
            &context,
            img_coarse_chan_idxs,
            &img_timestep_idxs,
            Some(&mut baseline_flagmasks),
        );

        // use the default strategy file location for MWA
        let strategy_filename = &aoflagger.FindStrategyFileMWA();

        // run the strategy on the imagesets, and get the resulting flagmasks for each baseline
        flag_imgsets_existing(
            &aoflagger,
            strategy_filename,
            &baseline_imgsets,
            &mut baseline_flagmasks,
            true,
        );

        // write the visibilities to disk as .uvfits

        let baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();
        write_uvfits(
            uvfits_out.as_path(),
            &context,
            &baseline_idxs,
            &baseline_imgsets,
            &baseline_flagmasks,
            &img_timestep_idxs,
            img_coarse_chan_idxs,
            None,
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

        for (row_idx, exp_group_params, exp_vis) in expected.iter() {
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
