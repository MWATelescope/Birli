//! BYO Array Example
//! This example demonstrates how to use the `birli` library to preprocess and write out a UVFITS file
//! from a set of gpubox files and a metafits file.
//! It uses the `mwalib` and `marlu` libraries to read in the gpubox files and metafits file.
//! It then uses the `birli` library to preprocess the data and write out a UVFITS file.
//!
//! ```bash
//! cargo run --example preprocess \
//!     output.uvfits \
//!     tests/data/1196175296_mwa_ord/1196175296.metafits \
//!     tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits
//! ```

use std::path::PathBuf;

use birli::{
    flag_to_weight_array,
    flags::get_weight_factor,
    io::{read_mwalib, write_uvfits},
    marlu::{
        constants::{
            COTTER_MWA_HEIGHT_METRES, COTTER_MWA_LATITUDE_RADIANS, COTTER_MWA_LONGITUDE_RADIANS,
        },
        mwalib::CorrelatorContext,
        LatLngHeight, RADec,
    },
    FlagContext, PreprocessContext, VisSelection,
};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let uvfits_path = PathBuf::from(&args[1]);
    let metafits_path = PathBuf::from(&args[2]);
    let gpubox_filenames: Vec<PathBuf> = args[3..].iter().map(PathBuf::from).collect();

    let corr_ctx = CorrelatorContext::new(metafits_path, &gpubox_filenames)
        .expect("unable to get mwalib context");

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

    // generate weights
    let weight_factor = get_weight_factor(&corr_ctx);
    let mut weight_array = flag_to_weight_array(flag_array.view(), weight_factor);

    let prep_ctx = PreprocessContext {
        array_pos: LatLngHeight {
            longitude_rad: COTTER_MWA_LONGITUDE_RADIANS,
            latitude_rad: COTTER_MWA_LATITUDE_RADIANS,
            height_metres: COTTER_MWA_HEIGHT_METRES,
        },
        phase_centre: RADec::from_mwalib_phase_or_pointing(&corr_ctx.metafits_context),
        correct_van_vleck: true,
        correct_cable_lengths: false,
        correct_digital_gains: false,
        correct_geometry: false,
        draw_progress: false,
        passband_gains: None,
        calsols: None,
        #[cfg(feature = "aoflagger")]
        aoflagger_strategy: None,
    };

    prep_ctx
        .preprocess(
            &corr_ctx,
            jones_array.view_mut(),
            weight_array.view_mut(),
            flag_array.view_mut(),
            &vis_sel,
        )
        .unwrap();

    let (avg_time, avg_freq) = (1, 1);

    write_uvfits(
        &uvfits_path,
        &corr_ctx,
        jones_array.view(),
        weight_array.view(),
        &vis_sel.timestep_range,
        &vis_sel.coarse_chan_range,
        &vis_sel.baseline_idxs,
        Some(prep_ctx.array_pos),
        Some(prep_ctx.phase_centre),
        avg_time,
        avg_freq,
    )
    .unwrap();
}
