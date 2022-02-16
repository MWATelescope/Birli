//! Corrections that can be performed on visibility data

use crate::{
    ndarray::{parallel::prelude::*, Array2, Array3, Axis},
    BirliError, Jones,
};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use itertools::izip;
use log::trace;
use marlu::{
    constants::VEL_C, hifitime::Epoch, mwalib::CorrelatorContext, precession::precess_time,
    Complex, LatLngHeight, RADec, XyzGeodetic, UVW,
};
use std::{f64::consts::PI, ops::Range};

/// Perform cable length corrections, given an observation's
/// [`mwalib::CorrelatorContext`] and an ['ndarray::Array3`] of [`TestJones`]
/// visibilities
///
/// Cable lengths are determined by the difference between a baseline's rfInput
/// electrical lengths in the metafits. Complex visibilities are phase-shifted
/// by an angle determined by the electrical length, and the channel's
/// frequency.
///
/// # Examples
///
/// ```rust
/// use birli::{context_to_jones_array, correct_cable_lengths, mwalib::CorrelatorContext};
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
/// // Determine which timesteps and coarse channels we want to use
/// let sel_coarse_chan_idxs = &context.common_coarse_chan_indices;
/// let sel_timestep_idxs = &context.common_timestep_indices;
/// let baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();
///
/// let sel_timestep_range =
///     *sel_timestep_idxs.first().unwrap()..(*sel_timestep_idxs.last().unwrap() + 1);
/// let sel_coarse_chan_range =
///     *sel_coarse_chan_idxs.first().unwrap()..(*sel_coarse_chan_idxs.last().unwrap() + 1);
///
/// // read visibilities out of the gpubox files
/// let (mut jones_array, _) = context_to_jones_array(
///     &context,
///     &sel_timestep_range,
///     &sel_coarse_chan_range,
///     None,
///     false,
/// ).unwrap();
///
/// correct_cable_lengths(&context, &mut jones_array, &sel_coarse_chan_range, false);
/// ```
///
/// # Accuracy
///
/// This follows the Cotter implementation of cable length correction, however
/// there is a slower but more accurate version of the calculation which
/// uses f64 values for the sin_cos. According to benchmarks, the Cotter
/// implementation is about 32% faster (5.9 seconds vs 8.6) than the more
/// precise implementation, and they vary by about three parts in four million.
/// Therefore it was decided that the Cotter implementation was more favourable.
///
/// # Assumptions
///
/// - an Antenna's rfinput is the same for X and Y
pub fn correct_cable_lengths(
    context: &CorrelatorContext,
    jones_array: &mut Array3<Jones<f32>>,
    mwalib_coarse_chan_range: &Range<usize>,
    // TODO: allow subset of baselines
    // mwalib_baseline_idxs: &[usize],
    draw_progress: bool,
) {
    trace!("start correct_cable_lengths");

    let baselines = &context.metafits_context.baselines;
    let antennas = &context.metafits_context.antennas;

    let all_freqs_hz =
        context.get_fine_chan_freqs_hz_array(&mwalib_coarse_chan_range.clone().collect::<Vec<_>>());
    let jones_dims = jones_array.dim();

    let draw_target = if draw_progress {
        ProgressDrawTarget::stderr()
    } else {
        ProgressDrawTarget::hidden()
    };

    // Create a progress bar to show the status of the correction
    let correction_progress = ProgressBar::with_draw_target(jones_dims.2 as u64, draw_target);
    correction_progress.set_style(
        ProgressStyle::default_bar()
            .template(
                "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
            )
            .progress_chars("=> "),
    );
    correction_progress.set_message("cable corrections");

    jones_array
        .axis_iter_mut(Axis(2))
        .into_par_iter()
        .zip(baselines)
        .for_each(|(mut jones_baseline_view, baseline)| {
            if baseline.ant1_index == baseline.ant2_index {
                return;
            }

            let ant1 = &antennas[baseline.ant1_index];
            let ant2 = &antennas[baseline.ant2_index];

            let pol_lengths = [
                ant2.rfinput_x.electrical_length_m - ant1.rfinput_x.electrical_length_m,
                ant2.rfinput_y.electrical_length_m - ant1.rfinput_y.electrical_length_m,
                ant2.rfinput_y.electrical_length_m - ant1.rfinput_x.electrical_length_m,
                ant2.rfinput_x.electrical_length_m - ant1.rfinput_y.electrical_length_m,
            ];

            for (mut jones_chan_view, freq_hz) in jones_baseline_view
                .axis_iter_mut(Axis(1))
                .zip(all_freqs_hz.clone())
            {
                let pol_sin_cos = pol_lengths
                    .iter()
                    .map(|electrical_length_m| {
                        let angle: f64 = -2.0 * PI * electrical_length_m * (freq_hz as f64) / VEL_C;
                        let (sin_angle_f64, cos_angle_f64) = angle.sin_cos();
                        Complex::new(cos_angle_f64 as f32, sin_angle_f64 as f32)
                    })
                    .collect::<Vec<_>>();

                for jones in jones_chan_view.iter_mut() {
                    for (pol_complex, rotation) in jones.iter_mut().zip(pol_sin_cos.iter()) {
                        *pol_complex *= rotation;
                    }
                }
            }
            correction_progress.inc(1);
        });

    correction_progress.finish();

    trace!("end correct_cable_lengths");
}

/// Perform geometric corrections, given an observation's
/// [`mwalib::CorrelatorContext`] and an ['ndarray::Array3`] of [`TestJones`]
/// visibilities
///
/// Complex visibilities are phase-shifted by an angle determined by the length
/// of the w-coordinate for the baseline and the channel's frequency.
///
/// # Examples
///
/// ```rust
/// use birli::{context_to_jones_array, correct_geometry, mwalib::CorrelatorContext};
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
/// // Determine which timesteps and coarse channels we want to use
/// let sel_coarse_chan_idxs = &context.common_coarse_chan_indices;
/// let sel_timestep_idxs = &context.common_timestep_indices;
/// let baseline_idxs = (0..context.metafits_context.num_baselines).collect::<Vec<_>>();
///
/// let sel_timestep_range =
///     *sel_timestep_idxs.first().unwrap()..(*sel_timestep_idxs.last().unwrap() + 1);
/// let sel_coarse_chan_range =
///     *sel_coarse_chan_idxs.first().unwrap()..(*sel_coarse_chan_idxs.last().unwrap() + 1);
///
/// // read visibilities out of the gpubox files
/// let (mut jones_array, _) = context_to_jones_array(
///     &context,
///     &sel_timestep_range,
///     &sel_coarse_chan_range,
///     None,
///     false,
/// ).unwrap();
///
/// correct_geometry(
///     &context,
///     &mut jones_array,
///     &sel_timestep_range,
///     &sel_coarse_chan_range,
///     None,
///     None,
///     false,
/// );
/// ```
#[allow(clippy::too_many_arguments)]
pub fn correct_geometry(
    context: &CorrelatorContext,
    jones_array: &mut Array3<Jones<f32>>,
    mwalib_timestep_range: &Range<usize>,
    mwalib_coarse_chan_range: &Range<usize>,
    // TODO: allow subset of baselines
    // mwalib_baseline_idxs: &[usize],
    array_pos: Option<LatLngHeight>,
    phase_centre_ra: Option<RADec>,
    draw_progress: bool,
) {
    trace!("start correct_geometry");

    let array_pos = match array_pos {
        Some(pos) => pos,
        None => {
            // The results here are slightly different to those given by cotter.
            // This is at least partly due to different constants (the altitude is
            // definitely slightly different), but possibly also because ERFA is
            // more accurate than cotter's "homebrewed" Geodetic2XYZ.
            LatLngHeight::new_mwa()
        }
    };

    let timesteps = &context.timesteps[mwalib_timestep_range.clone()];

    let baselines = &context.metafits_context.baselines;

    let all_freqs_hz =
        context.get_fine_chan_freqs_hz_array(&mwalib_coarse_chan_range.clone().collect::<Vec<_>>());
    let jones_dims = jones_array.dim();

    let integration_time_s = context.metafits_context.corr_int_time_ms as f64 / 1000.0;

    let phase_centre_ra = match phase_centre_ra {
        Some(pc) => pc,
        None => RADec::from_mwalib_phase_or_pointing(&context.metafits_context),
    };
    let tiles_xyz_geod = XyzGeodetic::get_tiles(&context.metafits_context, array_pos.latitude_rad);

    // Create a progress bar to show the status of the correction
    let draw_target = if draw_progress {
        ProgressDrawTarget::stderr()
    } else {
        ProgressDrawTarget::hidden()
    };
    let correction_progress = ProgressBar::with_draw_target(jones_dims.0 as u64, draw_target);
    correction_progress.set_style(
        ProgressStyle::default_bar()
            .template(
                "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
            )
            .progress_chars("=> "),
    );
    correction_progress.set_message("geom corrections");

    jones_array
        .outer_iter_mut()
        .into_par_iter()
        .zip(timesteps)
        .for_each(|(mut jones_timestep_view, timestep)| {
            let epoch = Epoch::from_gpst_seconds(
                timestep.gps_time_ms as f64 / 1000.0 + integration_time_s / 2.0,
            );

            let prec_info = precess_time(
                phase_centre_ra,
                epoch,
                array_pos.longitude_rad,
                array_pos.latitude_rad,
            );

            let tiles_xyz_precessed = prec_info.precess_xyz_parallel(&tiles_xyz_geod);

            for (mut jones_baseline_view, baseline) in
                jones_timestep_view.axis_iter_mut(Axis(1)).zip(baselines)
            {
                let ant1_idx = baseline.ant1_index;
                let ant2_idx = baseline.ant2_index;

                let baseline_xyz_precessed =
                    tiles_xyz_precessed[ant1_idx] - tiles_xyz_precessed[ant2_idx];
                let uvw = UVW::from_xyz(baseline_xyz_precessed, prec_info.hadec_j2000);

                for (jones, freq_hz) in jones_baseline_view.iter_mut().zip(all_freqs_hz.clone()) {
                    let angle = -2.0 * PI * uvw.w * (freq_hz as f64) / VEL_C;
                    let (sin_angle_f64, cos_angle_f64) = angle.sin_cos();
                    *jones *= Complex::new(cos_angle_f64 as f32, sin_angle_f64 as f32);
                }
            }

            correction_progress.inc(1);
        });

    correction_progress.finish();

    trace!("end correct_geometry");
}

/// Apply corrections for digital gains for each coarse channel from values in metafits.
///
/// The channels provided in `jones_array` should correspond to the selected coarse channel range in
/// `sel_coarse_chan_range`, such that `jones_array.dim().1 = num_fine_chans_per_coarse * sel_coarse_chan_range.len()`
///
/// # Arguments
///
/// - `context` - The correlator context.
/// - `jones_array` - The array of Jones matrices to be corrected, [timestep][channel][baseleine].
/// - `sel_coarse_chan_range` - The range of mwalib coarse channels which are used in in the channel
///     dimension of the jones array.
///
/// # Assumptions
/// - the digital gains are provided in [`mwalib::Rfinput.digital_gains`] in the same order as the
///   coarse channel indices (increasing sky frequency)
///
/// # Errors
/// - Will throw BadArrayShape if:
///     - `jones_array.dim().1 != num_fine_chans_per_coarse * sel_coarse_chan_range.len()
///     - `jones_array.dim().2 != sel_baselines.len()`
pub fn correct_digital_gains(
    context: &CorrelatorContext,
    jones_array: &mut Array3<Jones<f32>>,
    sel_coarse_chan_range: &Range<usize>,
    sel_baselines: &[(usize, usize)],
) -> Result<(), BirliError> {
    let num_fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;

    let vis_dims = jones_array.dim();
    if vis_dims.1 != sel_coarse_chan_range.len() * num_fine_chans_per_coarse {
        return Err(BirliError::BadArrayShape {
            argument: "sel_coarse_chan_range".into(),
            function: "correct_digital_gains".into(),
            expected: format!(
                "(vis_dims.1={}) / (num_fine_chans_per_coarse={}) = {}",
                vis_dims.1,
                num_fine_chans_per_coarse,
                vis_dims.1 / num_fine_chans_per_coarse
            ),
            received: format!("{:?}", sel_coarse_chan_range.len()),
        });
    }
    if vis_dims.2 != sel_baselines.len() {
        return Err(BirliError::BadArrayShape {
            argument: "sel_baselines".into(),
            function: "correct_digital_gains".into(),
            expected: format!("vis_dims.2={}", vis_dims.2,),
            received: format!("{:?}", sel_baselines.len()),
        });
    }
    assert!(vis_dims.2 == sel_baselines.len());

    let gains = Array2::from_shape_fn(
        (
            context.metafits_context.num_ants,
            sel_coarse_chan_range.len(),
        ),
        |(ant_idx, coarse_chan_idx)| {
            let ant = context.metafits_context.antennas[ant_idx].clone();
            let mwalib_coarse_chan_idx = sel_coarse_chan_range.start + coarse_chan_idx;
            (
                ant.rfinput_x.digital_gains[mwalib_coarse_chan_idx],
                ant.rfinput_y.digital_gains[mwalib_coarse_chan_idx],
            )
        },
    );

    _correct_digital_gains(
        jones_array,
        &gains,
        sel_baselines,
        num_fine_chans_per_coarse,
    )
}

fn _correct_digital_gains(
    jones_array: &mut Array3<Jones<f32>>,
    gains: &Array2<(f64, f64)>,
    sel_baselines: &[(usize, usize)],
    num_fine_chans_per_coarse: usize,
) -> Result<(), BirliError> {
    // TODO: proper error handling
    let vis_dims = jones_array.dim();
    let gain_dims = gains.dim();
    if vis_dims.1 != gain_dims.1 * num_fine_chans_per_coarse {
        return Err(BirliError::BadArrayShape {
            argument: "gains".into(),
            function: "_correct_digital_gains".into(),
            expected: format!(
                "(_, n, _), where n = (vis_dims.1={}) / (num_fine_chans_per_coarse={}) = {}",
                vis_dims.1,
                num_fine_chans_per_coarse,
                vis_dims.1 / num_fine_chans_per_coarse
            ),
            received: format!("{:?}", gain_dims),
        });
    }
    assert!(vis_dims.1 == gain_dims.1 * num_fine_chans_per_coarse);

    // iterate through the selected baselines
    for (mut jones_array, (ant1_idx, ant2_idx)) in izip!(
        jones_array.axis_iter_mut(Axis(2)),
        sel_baselines.iter().cloned()
    ) {
        // iterate through the selected coarse channels
        for (mut jones_array, &(gain1x, gain1y), &(gain2x, gain2y)) in izip!(
            jones_array.axis_chunks_iter_mut(Axis(1), num_fine_chans_per_coarse),
            gains.index_axis(Axis(0), ant1_idx),
            gains.index_axis(Axis(0), ant2_idx),
        ) {
            // for all visibilities in the selected coarse channel, for all timesteps
            for jones in jones_array.iter_mut() {
                // promote
                let corrected = Jones::<f64>::from(*jones);
                // divide by gain and demote
                *jones = Jones::<f32>::from([
                    (corrected[0].re / gain1x / gain2x) as _,
                    (corrected[0].im / gain1x / gain2x) as _,
                    (corrected[1].re / gain1x / gain2y) as _,
                    (corrected[1].im / gain1x / gain2y) as _,
                    (corrected[2].re / gain1y / gain2x) as _,
                    (corrected[2].im / gain1y / gain2x) as _,
                    (corrected[3].re / gain1y / gain2y) as _,
                    (corrected[3].im / gain1y / gain2y) as _,
                ]);
            }
        }
    }

    Ok(())
}

/// Correct for pfb filter gains in each coarse channel by scaling coarse_band_gains to fit a coarse band.
///
/// # Arguments
///
/// - `jones_array` - The array of Jones matrices to be corrected, [timestep][channel][baseleine].
/// - `coarse_band_gains` - a slice of gains to be applied to each coarse channel.
/// - `num_fine_chans_per_coarse` - The number of fine channels in each coarse.
///
/// # Errors
///
/// Will throw `BirliError::BadArrayShape` if:
/// - The length of the channel axis in `jones_array` is not a multiple of num_fine_chans_per_coarse.
/// - The length of the coarse band gains is not a multiple of num_fine_chans_per_coarse, or vice versa.
pub fn correct_pfb_gains(
    jones_array: &mut Array3<Jones<f32>>,
    coarse_band_gains: &[f64],
    num_fine_chans_per_coarse: usize,
) -> Result<(), BirliError> {
    if jones_array.dim().1 % num_fine_chans_per_coarse != 0 {
        return Err(BirliError::BadArrayShape {
            argument: "jones_array".into(),
            function: "correct_pfb_gains".into(),
            expected: format!(
                "(_, n, _), where n is a multiple of num_fine_chans_per_coarse={}",
                num_fine_chans_per_coarse
            ),
            received: format!("{:?}", jones_array.dim()),
        });
    };

    let fscrunch = if num_fine_chans_per_coarse % coarse_band_gains.len() == 0 {
        num_fine_chans_per_coarse / coarse_band_gains.len()
    } else {
        return Err(BirliError::BadArrayShape {
            argument: "coarse_band_gains".into(),
            function: "correct_pfb_gains".into(),
            expected: format!(
                "n, where n is a multiple of num_fine_chans_per_coarse={}",
                num_fine_chans_per_coarse
            ),
            received: format!("{:?}", coarse_band_gains.len()),
        });
    };

    let scrunched_gains: Vec<f64> = if fscrunch == 1 {
        coarse_band_gains.to_vec()
    } else {
        coarse_band_gains
            .chunks(fscrunch)
            .map(|gains_chunk| gains_chunk.iter().sum::<f64>() / (gains_chunk.len() as f64))
            .collect()
    };

    for mut jones_array in jones_array.axis_chunks_iter_mut(Axis(1), num_fine_chans_per_coarse) {
        for (mut jones_array, &gain) in
            izip!(jones_array.axis_iter_mut(Axis(1)), scrunched_gains.iter())
        {
            for jones in jones_array.iter_mut() {
                let jones_f64 = Jones::<f64>::from(*jones);
                *jones = Jones::<f32>::from(jones_f64 / gain);
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    #![allow(clippy::float_cmp)]

    use super::{
        correct_cable_lengths, correct_digital_gains, correct_geometry, correct_pfb_gains, VEL_C,
    };
    use itertools::izip;
    use marlu::{
        hifitime::Epoch, mwalib::CorrelatorContext, precession::precess_time, Complex, Jones,
        LatLngHeight, RADec, XyzGeodetic, UVW,
    };
    use ndarray::{s, Array2, Array3, Axis};
    use std::f64::consts::PI;

    use crate::{
        approx::assert_abs_diff_eq,
        context_to_jones_array,
        corrections::_correct_digital_gains,
        flags::{get_coarse_chan_range, get_timestep_range},
        get_flaggable_timesteps, BirliError, TestJones,
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

    #[test]
    fn test_cable_length_corrections_mwax() {
        let context = get_mwax_context();

        let sel_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(sel_timestep_idxs.len(), 4);
        let sel_timestep_range =
            *sel_timestep_idxs.first().unwrap()..(*sel_timestep_idxs.last().unwrap() + 1);
        let sel_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let sel_coarse_chan_range =
            *sel_coarse_chan_idxs.first().unwrap()..(*sel_coarse_chan_idxs.last().unwrap() + 1);

        let all_freqs_hz = context.get_fine_chan_freqs_hz_array(sel_coarse_chan_idxs);

        // let sel_baseline_idxs: Vec<usize> = (0..context.metafits_context.num_baselines).collect();

        let (jones_array, _) = context_to_jones_array(
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            None,
            false,
        )
        .unwrap();

        let jones_array = jones_array.mapv(TestJones::from);

        // ts 0, chan 0 (cc 0, fc 0), baseline 0
        let viz_0_0_0 = *jones_array.get((0, 0, 0)).unwrap();
        assert_abs_diff_eq!(
            viz_0_0_0,
            TestJones::from([
                Complex::new(0x410000 as f32, 0x410001 as f32),
                Complex::new(0x410002 as f32, 0x410003 as f32),
                Complex::new(0x410004 as f32, 0x410005 as f32),
                Complex::new(0x410006 as f32, 0x410007 as f32),
            ])
        );

        // ts 0, chan 0 (cc 0, fc 0), baseline 1
        let viz_0_0_1 = *jones_array.get((0, 0, 1)).unwrap();
        assert_abs_diff_eq!(
            viz_0_0_1,
            TestJones::from([
                Complex::new(0x410010 as f32, 0x410011 as f32),
                Complex::new(0x410012 as f32, 0x410013 as f32),
                Complex::new(0x410014 as f32, 0x410015 as f32),
                Complex::new(0x410016 as f32, 0x410017 as f32),
            ])
        );

        // ts 3, chan 3 (cc 1, fc 1), baseline 1
        let viz_3_3_1 = *jones_array.get((3, 3, 1)).unwrap();
        assert_abs_diff_eq!(
            viz_3_3_1,
            TestJones::from([
                Complex::new(0x410718 as f32, 0x410719 as f32),
                Complex::new(0x41071a as f32, 0x41071b as f32),
                Complex::new(0x41071c as f32, 0x41071d as f32),
                Complex::new(0x41071e as f32, 0x41071f as f32),
            ])
        );

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

        let length_m_1_xx = length_1_2_x - length_1_1_x;
        let length_m_1_xy = length_1_2_y - length_1_1_x;
        let length_m_1_yx = length_1_2_x - length_1_1_y;
        let length_m_1_yy = length_1_2_y - length_1_1_y;

        // baseline 1, pol XX, chan 0 (cc 0, fc 0)
        let angle_1_xx_0: f64 = -2.0 * PI * length_m_1_xx * (all_freqs_hz[0] as f64) / VEL_C;
        let (sin_1_xx_0_f64, cos_1_xx_0_f64) = angle_1_xx_0.sin_cos();
        let (sin_1_xx_0, cos_1_xx_0) = (sin_1_xx_0_f64 as f32, cos_1_xx_0_f64 as f32);
        // baseline 1, pol XY, chan 0 (cc 0, fc 0)
        let angle_1_xy_0: f64 = -2.0 * PI * length_m_1_xy * (all_freqs_hz[0] as f64) / VEL_C;
        let (sin_1_xy_0_f64, cos_1_xy_0_f64) = angle_1_xy_0.sin_cos();
        let (sin_1_xy_0, cos_1_xy_0) = (sin_1_xy_0_f64 as f32, cos_1_xy_0_f64 as f32);
        // baseline 1, pol YX, chan 0 (cc 0, fc 0)
        let angle_1_yx_0: f64 = -2.0 * PI * length_m_1_yx * (all_freqs_hz[0] as f64) / VEL_C;
        let (sin_1_yx_0_f64, cos_1_yx_0_f64) = angle_1_yx_0.sin_cos();
        let (sin_1_yx_0, cos_1_yx_0) = (sin_1_yx_0_f64 as f32, cos_1_yx_0_f64 as f32);
        // baseline 1, pol YY, chan 0 (cc 0, fc 0)
        let angle_1_yy_0: f64 = -2.0 * PI * length_m_1_yy * (all_freqs_hz[0] as f64) / VEL_C;
        let (sin_1_yy_0_f64, cos_1_yy_0_f64) = angle_1_yy_0.sin_cos();
        let (sin_1_yy_0, cos_1_yy_0) = (sin_1_yy_0_f64 as f32, cos_1_yy_0_f64 as f32);
        // baseline 1, pol XX, chan 3 (cc 1, fc 1)
        let angle_1_xx_3: f64 = -2.0 * PI * length_m_1_xx * (all_freqs_hz[3] as f64) / VEL_C;
        let (sin_1_xx_3_f64, cos_1_xx_3_f64) = angle_1_xx_3.sin_cos();
        let (sin_1_xx_3, cos_1_xx_3) = (sin_1_xx_3_f64 as f32, cos_1_xx_3_f64 as f32);
        // baseline 1, pol XY, chan 3 (cc 1, fc 1)
        let angle_1_xy_3: f64 = -2.0 * PI * length_m_1_xy * (all_freqs_hz[3] as f64) / VEL_C;
        let (sin_1_xy_3_f64, cos_1_xy_3_f64) = angle_1_xy_3.sin_cos();
        let (sin_1_xy_3, cos_1_xy_3) = (sin_1_xy_3_f64 as f32, cos_1_xy_3_f64 as f32);
        // baseline 1, pol YX, chan 3 (cc 1, fc 1)
        let angle_1_yx_3: f64 = -2.0 * PI * length_m_1_yx * (all_freqs_hz[3] as f64) / VEL_C;
        let (sin_1_yx_3_f64, cos_1_yx_3_f64) = angle_1_yx_3.sin_cos();
        let (sin_1_yx_3, cos_1_yx_3) = (sin_1_yx_3_f64 as f32, cos_1_yx_3_f64 as f32);
        // baseline 1, pol YY, chan 3 (cc 1, fc 1)
        let angle_1_yy_3: f64 = -2.0 * PI * length_m_1_yy * (all_freqs_hz[3] as f64) / VEL_C;
        let (sin_1_yy_3_f64, cos_1_yy_3_f64) = angle_1_yy_3.sin_cos();
        let (sin_1_yy_3, cos_1_yy_3) = (sin_1_yy_3_f64 as f32, cos_1_yy_3_f64 as f32);

        // TODO: this is not great.
        let mut jones_array = jones_array.mapv(|e| Jones::from([e[0], e[1], e[2], e[3]]));
        correct_cable_lengths(&context, &mut jones_array, &sel_coarse_chan_range, false);

        let jones_array = jones_array.mapv(TestJones::from);

        // there should be no difference in baseline 0
        // ts 0, chan 0, baseline 0
        assert_abs_diff_eq!(*jones_array.get((0, 0, 0)).unwrap(), viz_0_0_0);

        ////
        // baseline 1 should be rotated
        ////
        // ts 0, chan 0 (cc 0, fc 0), baseline 1
        // -> pol XX
        let rot_1_xx_0_0_re = (cos_1_xx_0 * viz_0_0_1[0].re - sin_1_xx_0 * viz_0_0_1[0].im) as f32;
        let rot_1_xx_0_0_im = (sin_1_xx_0 * viz_0_0_1[0].re + cos_1_xx_0 * viz_0_0_1[0].im) as f32;
        // -> pol XY
        let rot_1_xy_0_0_re = (cos_1_xy_0 * viz_0_0_1[1].re - sin_1_xy_0 * viz_0_0_1[1].im) as f32;
        let rot_1_xy_0_0_im = (sin_1_xy_0 * viz_0_0_1[1].re + cos_1_xy_0 * viz_0_0_1[1].im) as f32;
        // -> pol YX
        let rot_1_yx_0_0_re = (cos_1_yx_0 * viz_0_0_1[2].re - sin_1_yx_0 * viz_0_0_1[2].im) as f32;
        let rot_1_yx_0_0_im = (sin_1_yx_0 * viz_0_0_1[2].re + cos_1_yx_0 * viz_0_0_1[2].im) as f32;
        // -> pol YY
        let rot_1_yy_0_0_re = (cos_1_yy_0 * viz_0_0_1[3].re - sin_1_yy_0 * viz_0_0_1[3].im) as f32;
        let rot_1_yy_0_0_im = (sin_1_yy_0 * viz_0_0_1[3].re + cos_1_yy_0 * viz_0_0_1[3].im) as f32;

        assert_abs_diff_eq!(
            *jones_array.get((0, 0, 1)).unwrap(),
            &TestJones::from([
                Complex::new(rot_1_xx_0_0_re, rot_1_xx_0_0_im),
                Complex::new(rot_1_xy_0_0_re, rot_1_xy_0_0_im),
                Complex::new(rot_1_yx_0_0_re, rot_1_yx_0_0_im),
                Complex::new(rot_1_yy_0_0_re, rot_1_yy_0_0_im),
            ])
        );

        // ts 3 (batch 1, scan 1), chan 3 (cc 1, fc 1), baseline 1
        // -> pol XX
        let rot_1_xx_3_3_re = (cos_1_xx_3 * viz_3_3_1[0].re - sin_1_xx_3 * viz_3_3_1[0].im) as f32;
        let rot_1_xx_3_3_im = (sin_1_xx_3 * viz_3_3_1[0].re + cos_1_xx_3 * viz_3_3_1[0].im) as f32;
        // -> pol XY
        let rot_1_xy_3_3_re = (cos_1_xy_3 * viz_3_3_1[1].re - sin_1_xy_3 * viz_3_3_1[1].im) as f32;
        let rot_1_xy_3_3_im = (sin_1_xy_3 * viz_3_3_1[1].re + cos_1_xy_3 * viz_3_3_1[1].im) as f32;
        // -> pol YX
        let rot_1_yx_3_3_re = (cos_1_yx_3 * viz_3_3_1[2].re - sin_1_yx_3 * viz_3_3_1[2].im) as f32;
        let rot_1_yx_3_3_im = (sin_1_yx_3 * viz_3_3_1[2].re + cos_1_yx_3 * viz_3_3_1[2].im) as f32;
        // -> pol YY
        let rot_1_yy_3_3_re = (cos_1_yy_3 * viz_3_3_1[3].re - sin_1_yy_3 * viz_3_3_1[3].im) as f32;
        let rot_1_yy_3_3_im = (sin_1_yy_3 * viz_3_3_1[3].re + cos_1_yy_3 * viz_3_3_1[3].im) as f32;

        assert_abs_diff_eq!(
            *jones_array.get((3, 3, 1)).unwrap(),
            &TestJones::from([
                Complex::new(rot_1_xx_3_3_re, rot_1_xx_3_3_im),
                Complex::new(rot_1_xy_3_3_re, rot_1_xy_3_3_im),
                Complex::new(rot_1_yx_3_3_re, rot_1_yx_3_3_im),
                Complex::new(rot_1_yy_3_3_re, rot_1_yy_3_3_im),
            ])
        );
    }

    #[test]
    fn test_cable_length_corrections_ord() {
        let context = get_mwa_ord_context();

        let sel_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(sel_timestep_idxs.len(), 4);
        let sel_timestep_range =
            *sel_timestep_idxs.first().unwrap()..(*sel_timestep_idxs.last().unwrap() + 1);
        let sel_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let sel_coarse_chan_range =
            *sel_coarse_chan_idxs.first().unwrap()..(*sel_coarse_chan_idxs.last().unwrap() + 1);

        let all_freqs_hz = context.get_fine_chan_freqs_hz_array(sel_coarse_chan_idxs);

        // let sel_baseline_idxs: Vec<usize> = (0..context.metafits_context.num_baselines).collect();

        let (jones_array, _) = context_to_jones_array(
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            None,
            false,
        )
        .unwrap();

        let jones_array = jones_array.mapv(TestJones::from);

        // ts 0, chan 0 (cc 0, fc 0), baseline 0
        let viz_0_0_0 = *jones_array.get((0, 0, 0)).unwrap();
        assert_abs_diff_eq!(
            viz_0_0_0,
            TestJones::from([
                Complex::new(0x10c5be as f32, -0x10c5bf as f32),
                Complex::new(0x10c5ae as f32, 0x10c5af as f32),
                Complex::new(0x10c5ae as f32, -0x10c5af as f32),
                Complex::new(0x10bec6 as f32, -0x10bec7 as f32),
            ])
        );

        // ts 0, chan 0 (cc 0, fc 0), baseline 5
        let viz_0_0_5 = *jones_array.get((0, 0, 5)).unwrap();
        assert_abs_diff_eq!(
            viz_0_0_5,
            TestJones::from([
                Complex::new(0x10f1ce as f32, -0x10f1cf as f32),
                Complex::new(0x10ea26 as f32, -0x10ea27 as f32),
                Complex::new(0x10f1be as f32, -0x10f1bf as f32),
                Complex::new(0x10ea16 as f32, -0x10ea17 as f32),
            ])
        );

        // ts 3, chan 3 (cc 1, fc 1), baseline 5
        let viz_3_3_5 = *jones_array.get((3, 3, 5)).unwrap();
        assert_abs_diff_eq!(
            viz_3_3_5,
            TestJones::from([
                Complex::new(0x0df3ce as f32, -0x0df3cf as f32),
                Complex::new(0x0dec26 as f32, -0x0dec27 as f32),
                Complex::new(0x0df3be as f32, -0x0df3bf as f32),
                Complex::new(0x0dec16 as f32, -0x0dec17 as f32),
            ])
        );

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

        let length_m_5_xx = length_5_2_x - length_5_1_x;
        let length_m_5_xy = length_5_2_y - length_5_1_x;
        let length_m_5_yx = length_5_2_x - length_5_1_y;
        let length_m_5_yy = length_5_2_y - length_5_1_y;

        // baseline 5, pol XX, chan 0 (cc 0, fc 0)
        let angle_5_xx_0: f64 = -2.0 * PI * length_m_5_xx * (all_freqs_hz[0] as f64) / VEL_C;
        let (sin_5_xx_0_f64, cos_5_xx_0_f64) = angle_5_xx_0.sin_cos();
        let (sin_5_xx_0, cos_5_xx_0) = (sin_5_xx_0_f64 as f32, cos_5_xx_0_f64 as f32);
        // baseline 5, pol XY, chan 0 (cc 0, fc 0)
        let angle_5_xy_0: f64 = -2.0 * PI * length_m_5_xy * (all_freqs_hz[0] as f64) / VEL_C;
        let (sin_5_xy_0_f64, cos_5_xy_0_f64) = angle_5_xy_0.sin_cos();
        let (sin_5_xy_0, cos_5_xy_0) = (sin_5_xy_0_f64 as f32, cos_5_xy_0_f64 as f32);
        // baseline 5, pol YX, chan 0 (cc 0, fc 0)
        let angle_5_yx_0: f64 = -2.0 * PI * length_m_5_yx * (all_freqs_hz[0] as f64) / VEL_C;
        let (sin_5_yx_0_f64, cos_5_yx_0_f64) = angle_5_yx_0.sin_cos();
        let (sin_5_yx_0, cos_5_yx_0) = (sin_5_yx_0_f64 as f32, cos_5_yx_0_f64 as f32);
        // baseline 5, pol YY, chan 0 (cc 0, fc 0)
        let angle_5_yy_0: f64 = -2.0 * PI * length_m_5_yy * (all_freqs_hz[0] as f64) / VEL_C;
        let (sin_5_yy_0_f64, cos_5_yy_0_f64) = angle_5_yy_0.sin_cos();
        let (sin_5_yy_0, cos_5_yy_0) = (sin_5_yy_0_f64 as f32, cos_5_yy_0_f64 as f32);
        // baseline 5, pol XX, chan 3 (cc 1, fc 1)
        let angle_5_xx_3: f64 = -2.0 * PI * length_m_5_xx * (all_freqs_hz[3] as f64) / VEL_C;
        let (sin_5_xx_3_f64, cos_5_xx_3_f64) = angle_5_xx_3.sin_cos();
        let (sin_5_xx_3, cos_5_xx_3) = (sin_5_xx_3_f64 as f32, cos_5_xx_3_f64 as f32);
        // baseline 5, pol XY, chan 3 (cc 1, fc 1)
        let angle_5_xy_3: f64 = -2.0 * PI * length_m_5_xy * (all_freqs_hz[3] as f64) / VEL_C;
        let (sin_5_xy_3_f64, cos_5_xy_3_f64) = angle_5_xy_3.sin_cos();
        let (sin_5_xy_3, cos_5_xy_3) = (sin_5_xy_3_f64 as f32, cos_5_xy_3_f64 as f32);
        // baseline 5, pol YX, chan 3 (cc 1, fc 1)
        let angle_5_yx_3: f64 = -2.0 * PI * length_m_5_yx * (all_freqs_hz[3] as f64) / VEL_C;
        let (sin_5_yx_3_f64, cos_5_yx_3_f64) = angle_5_yx_3.sin_cos();
        let (sin_5_yx_3, cos_5_yx_3) = (sin_5_yx_3_f64 as f32, cos_5_yx_3_f64 as f32);
        // baseline 5, pol YY, chan 3 (cc 1, fc 1)
        let angle_5_yy_3: f64 = -2.0 * PI * length_m_5_yy * (all_freqs_hz[3] as f64) / VEL_C;
        let (sin_5_yy_3_f64, cos_5_yy_3_f64) = angle_5_yy_3.sin_cos();
        let (sin_5_yy_3, cos_5_yy_3) = (sin_5_yy_3_f64 as f32, cos_5_yy_3_f64 as f32);

        // FIXME: aaaaaaa
        let mut jones_array = jones_array.mapv(|e| Jones::from([e[0], e[1], e[2], e[3]]));

        correct_cable_lengths(&context, &mut jones_array, &sel_coarse_chan_range, false);

        let jones_array = jones_array.mapv(TestJones::from);

        // there should be no difference in baseline 0
        // ts 0 (batch 0, scan 0), chan 0 (cc 0, fc 0), baseline 0
        assert_abs_diff_eq!(*jones_array.get((0, 0, 0)).unwrap(), viz_0_0_0);

        ////
        // baseline 5 should be rotated
        ////
        // ts 0 (batch 0, scan 0), chan 0 (cc 0, fc 0), baseline 5
        // -> pol XX
        let rot_5_xx_0_0_re = (cos_5_xx_0 * viz_0_0_5[0].re - sin_5_xx_0 * viz_0_0_5[0].im) as f32;
        let rot_5_xx_0_0_im = (sin_5_xx_0 * viz_0_0_5[0].re + cos_5_xx_0 * viz_0_0_5[0].im) as f32;
        // -> pol XY
        let rot_5_xy_0_0_re = (cos_5_xy_0 * viz_0_0_5[1].re - sin_5_xy_0 * viz_0_0_5[1].im) as f32;
        let rot_5_xy_0_0_im = (sin_5_xy_0 * viz_0_0_5[1].re + cos_5_xy_0 * viz_0_0_5[1].im) as f32;
        // -> pol YX
        let rot_5_yx_0_0_re = (cos_5_yx_0 * viz_0_0_5[2].re - sin_5_yx_0 * viz_0_0_5[2].im) as f32;
        let rot_5_yx_0_0_im = (sin_5_yx_0 * viz_0_0_5[2].re + cos_5_yx_0 * viz_0_0_5[2].im) as f32;
        // -> pol YY
        let rot_5_yy_0_0_re = (cos_5_yy_0 * viz_0_0_5[3].re - sin_5_yy_0 * viz_0_0_5[3].im) as f32;
        let rot_5_yy_0_0_im = (sin_5_yy_0 * viz_0_0_5[3].re + cos_5_yy_0 * viz_0_0_5[3].im) as f32;

        assert_abs_diff_eq!(
            *jones_array.get((0, 0, 5)).unwrap(),
            &TestJones::from([
                Complex::new(rot_5_xx_0_0_re, rot_5_xx_0_0_im),
                Complex::new(rot_5_xy_0_0_re, rot_5_xy_0_0_im),
                Complex::new(rot_5_yx_0_0_re, rot_5_yx_0_0_im),
                Complex::new(rot_5_yy_0_0_re, rot_5_yy_0_0_im),
            ])
        );

        // ts 3 (batch 1, scan 1), chan 3 (cc 1, fc 1), baseline 5
        // -> pol XX
        let rot_5_xx_3_3_re = (cos_5_xx_3 * viz_3_3_5[0].re - sin_5_xx_3 * viz_3_3_5[0].im) as f32;
        let rot_5_xx_3_3_im = (sin_5_xx_3 * viz_3_3_5[0].re + cos_5_xx_3 * viz_3_3_5[0].im) as f32;
        // -> pol XY
        let rot_5_xy_3_3_re = (cos_5_xy_3 * viz_3_3_5[1].re - sin_5_xy_3 * viz_3_3_5[1].im) as f32;
        let rot_5_xy_3_3_im = (sin_5_xy_3 * viz_3_3_5[1].re + cos_5_xy_3 * viz_3_3_5[1].im) as f32;
        // -> pol YX
        let rot_5_yx_3_3_re = (cos_5_yx_3 * viz_3_3_5[2].re - sin_5_yx_3 * viz_3_3_5[2].im) as f32;
        let rot_5_yx_3_3_im = (sin_5_yx_3 * viz_3_3_5[2].re + cos_5_yx_3 * viz_3_3_5[2].im) as f32;
        // -> pol YY
        let rot_5_yy_3_3_re = (cos_5_yy_3 * viz_3_3_5[3].re - sin_5_yy_3 * viz_3_3_5[3].im) as f32;
        let rot_5_yy_3_3_im = (sin_5_yy_3 * viz_3_3_5[3].re + cos_5_yy_3 * viz_3_3_5[3].im) as f32;

        assert_abs_diff_eq!(
            *jones_array.get((3, 3, 5)).unwrap(),
            &TestJones::from([
                Complex::new(rot_5_xx_3_3_re, rot_5_xx_3_3_im),
                Complex::new(rot_5_xy_3_3_re, rot_5_xy_3_3_im),
                Complex::new(rot_5_yx_3_3_re, rot_5_yx_3_3_im),
                Complex::new(rot_5_yy_3_3_re, rot_5_yy_3_3_im),
            ])
        );
    }

    #[test]
    fn test_geometric_corrections_ord() {
        let context = get_mwa_ord_context();

        let sel_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(sel_timestep_idxs.len(), 4);
        let sel_timestep_range =
            *sel_timestep_idxs.first().unwrap()..(*sel_timestep_idxs.last().unwrap() + 1);
        let sel_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let sel_coarse_chan_range =
            *sel_coarse_chan_idxs.first().unwrap()..(*sel_coarse_chan_idxs.last().unwrap() + 1);

        let all_freqs_hz = context.get_fine_chan_freqs_hz_array(sel_coarse_chan_idxs);

        let array_pos = LatLngHeight::new_mwa();

        let phase_centre_ra = RADec::from_mwalib_phase_or_pointing(&context.metafits_context);
        // let lst_rad = context.metafits_context.lst_rad;
        // let phase_centre_ha = phase_centre_ra.to_hadec(lst_rad);
        let tiles_xyz_geod = XyzGeodetic::get_tiles_mwa(&context.metafits_context);

        let (jones_array, _) = context_to_jones_array(
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            None,
            false,
        )
        .unwrap();

        let jones_array = jones_array.mapv(TestJones::from);

        // ts 0, chan 0 (cc 0, fc 0), baseline 0
        let viz_0_0_0 = *jones_array.get((0, 0, 0)).unwrap();
        assert_abs_diff_eq!(
            viz_0_0_0,
            TestJones::from([
                Complex::new(0x10c5be as f32, -0x10c5bf as f32),
                Complex::new(0x10c5ae as f32, 0x10c5af as f32),
                Complex::new(0x10c5ae as f32, -0x10c5af as f32),
                Complex::new(0x10bec6 as f32, -0x10bec7 as f32),
            ])
        );

        // ts 0, chan 0 (cc 0, fc 0), baseline 5
        let viz_0_0_5 = *jones_array.get((0, 0, 5)).unwrap();
        assert_abs_diff_eq!(
            viz_0_0_5,
            TestJones::from([
                Complex::new(0x10f1ce as f32, -0x10f1cf as f32),
                Complex::new(0x10ea26 as f32, -0x10ea27 as f32),
                Complex::new(0x10f1be as f32, -0x10f1bf as f32),
                Complex::new(0x10ea16 as f32, -0x10ea17 as f32),
            ])
        );

        // ts 3, chan 3 (cc 1, fc 1), baseline 5
        let viz_3_3_5 = *jones_array.get((3, 3, 5)).unwrap();
        assert_abs_diff_eq!(
            viz_3_3_5,
            TestJones::from([
                Complex::new(0x0df3ce as f32, -0x0df3cf as f32),
                Complex::new(0x0dec26 as f32, -0x0dec27 as f32),
                Complex::new(0x0df3be as f32, -0x0df3bf as f32),
                Complex::new(0x0dec16 as f32, -0x0dec17 as f32),
            ])
        );

        let integration_time_s = context.metafits_context.corr_int_time_ms as f64 / 1000.0;

        // timestep 0
        let timestep_0 = &context.timesteps[sel_timestep_idxs[0]];
        let epoch_0 = Epoch::from_gpst_seconds(
            timestep_0.gps_time_ms as f64 / 1000.0 + integration_time_s / 2.0,
        );
        let prec_info_0 = precess_time(
            phase_centre_ra,
            epoch_0,
            array_pos.longitude_rad,
            array_pos.latitude_rad,
        );
        let phase_centre_ha_j2000_0 = prec_info_0.hadec_j2000; // phase_centre_ra.to_hadec(prec_info_0.lmst_j2000);
        let tiles_xyz_precessed_0 = prec_info_0.precess_xyz_parallel(&tiles_xyz_geod);
        // timestep 3
        let timestep_3 = &context.timesteps[sel_timestep_idxs[3]];
        let epoch_3 = Epoch::from_gpst_seconds(
            timestep_3.gps_time_ms as f64 / 1000.0 + integration_time_s / 2.0,
        );
        let prec_info_3 = precess_time(
            phase_centre_ra,
            epoch_3,
            array_pos.longitude_rad,
            array_pos.latitude_rad,
        );
        let phase_centre_ha_j2000_3 = prec_info_3.hadec_j2000; // phase_centre_ra.to_hadec(prec_info_3.lmst_j2000);
        let tiles_xyz_precessed_3 = prec_info_3.precess_xyz_parallel(&tiles_xyz_geod);

        // baseline 5
        let bl_5 = &context.metafits_context.baselines[5];

        // ts 0, baseline 5
        let xyz_0_5 =
            tiles_xyz_precessed_0[bl_5.ant1_index] - tiles_xyz_precessed_0[bl_5.ant2_index];
        let w_0_5 = UVW::from_xyz(xyz_0_5, phase_centre_ha_j2000_0).w;
        // ts 3, baseline 5
        let xyz_3_5 =
            tiles_xyz_precessed_3[bl_5.ant1_index] - tiles_xyz_precessed_3[bl_5.ant2_index];
        let w_3_5 = UVW::from_xyz(xyz_3_5, phase_centre_ha_j2000_3).w;

        // ts 0, chan 0 (cc 0, fc 0), baseline 5
        let angle_0_0_5 = -2.0 * PI * w_0_5 * (all_freqs_hz[0] as f64) / VEL_C;
        let (sin_0_0_5_f64, cos_0_0_5_f64) = angle_0_0_5.sin_cos();
        let (sin_0_0_5, cos_0_0_5) = (sin_0_0_5_f64 as f32, cos_0_0_5_f64 as f32);

        // ts 3, chan 3 (cc 1, fc 1), baseline 5
        let angle_3_3_5 = -2.0 * PI * w_3_5 * (all_freqs_hz[3] as f64) / VEL_C;
        let (sin_3_3_5_f64, cos_3_3_5_f64) = angle_3_3_5.sin_cos();
        let (sin_3_3_5, cos_3_3_5) = (sin_3_3_5_f64 as f32, cos_3_3_5_f64 as f32);

        // let angle_5_3 = -2.0 * PI * w_5_3 * (all_freqs_hz[3] as f64) / VEL_C;

        let mut jones_array = jones_array.mapv(|e| Jones::from([e[0], e[1], e[2], e[3]]));

        correct_geometry(
            &context,
            &mut jones_array,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            None,
            None,
            false,
        );

        let jones_array = jones_array.mapv(TestJones::from);

        // there should be no difference in baseline 0
        // ts 0 (batch 0, scan 0), chan 0 (cc 0, fc 0), baseline 0
        assert_abs_diff_eq!(*jones_array.get((0, 0, 0)).unwrap(), viz_0_0_0);

        ////
        // baseline 5 should be rotated
        ////
        // ts 0 (batch 0, scan 0), chan 0 (cc 0, fc 0), baseline 5
        // -> pol XX
        let rot_5_xx_0_0_re = (cos_0_0_5 * viz_0_0_5[0].re - sin_0_0_5 * viz_0_0_5[0].im) as f32;
        let rot_5_xx_0_0_im = (sin_0_0_5 * viz_0_0_5[0].re + cos_0_0_5 * viz_0_0_5[0].im) as f32;
        // -> pol XY
        let rot_5_xy_0_0_re = (cos_0_0_5 * viz_0_0_5[1].re - sin_0_0_5 * viz_0_0_5[1].im) as f32;
        let rot_5_xy_0_0_im = (sin_0_0_5 * viz_0_0_5[1].re + cos_0_0_5 * viz_0_0_5[1].im) as f32;
        // -> pol YX
        let rot_5_yx_0_0_re = (cos_0_0_5 * viz_0_0_5[2].re - sin_0_0_5 * viz_0_0_5[2].im) as f32;
        let rot_5_yx_0_0_im = (sin_0_0_5 * viz_0_0_5[2].re + cos_0_0_5 * viz_0_0_5[2].im) as f32;
        // -> pol YY
        let rot_5_yy_0_0_re = (cos_0_0_5 * viz_0_0_5[3].re - sin_0_0_5 * viz_0_0_5[3].im) as f32;
        let rot_5_yy_0_0_im = (sin_0_0_5 * viz_0_0_5[3].re + cos_0_0_5 * viz_0_0_5[3].im) as f32;

        assert_abs_diff_eq!(
            *jones_array.get((0, 0, 5)).unwrap(),
            &TestJones::from([
                Complex::new(rot_5_xx_0_0_re, rot_5_xx_0_0_im),
                Complex::new(rot_5_xy_0_0_re, rot_5_xy_0_0_im),
                Complex::new(rot_5_yx_0_0_re, rot_5_yx_0_0_im),
                Complex::new(rot_5_yy_0_0_re, rot_5_yy_0_0_im),
            ])
        );

        // ts 3 (batch 1, scan 1), chan 3 (cc 1, fc 1), baseline 5
        // -> pol XX
        let rot_5_xx_3_3_re = (cos_3_3_5 * viz_3_3_5[0].re - sin_3_3_5 * viz_3_3_5[0].im) as f32;
        let rot_5_xx_3_3_im = (sin_3_3_5 * viz_3_3_5[0].re + cos_3_3_5 * viz_3_3_5[0].im) as f32;
        // -> pol XY
        let rot_5_xy_3_3_re = (cos_3_3_5 * viz_3_3_5[1].re - sin_3_3_5 * viz_3_3_5[1].im) as f32;
        let rot_5_xy_3_3_im = (sin_3_3_5 * viz_3_3_5[1].re + cos_3_3_5 * viz_3_3_5[1].im) as f32;
        // -> pol YX
        let rot_5_yx_3_3_re = (cos_3_3_5 * viz_3_3_5[2].re - sin_3_3_5 * viz_3_3_5[2].im) as f32;
        let rot_5_yx_3_3_im = (sin_3_3_5 * viz_3_3_5[2].re + cos_3_3_5 * viz_3_3_5[2].im) as f32;
        // -> pol YY
        let rot_5_yy_3_3_re = (cos_3_3_5 * viz_3_3_5[3].re - sin_3_3_5 * viz_3_3_5[3].im) as f32;
        let rot_5_yy_3_3_im = (sin_3_3_5 * viz_3_3_5[3].re + cos_3_3_5 * viz_3_3_5[3].im) as f32;

        assert_abs_diff_eq!(
            *jones_array.get((3, 3, 5)).unwrap(),
            &TestJones::from([
                Complex::new(rot_5_xx_3_3_re, rot_5_xx_3_3_im),
                Complex::new(rot_5_xy_3_3_re, rot_5_xy_3_3_im),
                Complex::new(rot_5_yx_3_3_re, rot_5_yx_3_3_im),
                Complex::new(rot_5_yy_3_3_re, rot_5_yy_3_3_im),
            ])
        );
    }

    #[test]
    fn test_geometric_corrections_mwax() {
        let context = get_mwax_context();

        let sel_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        let sel_timestep_range =
            *sel_timestep_idxs.first().unwrap()..(*sel_timestep_idxs.last().unwrap() + 1);
        let sel_coarse_chan_idxs = &context.common_coarse_chan_indices;
        let sel_coarse_chan_range =
            *sel_coarse_chan_idxs.first().unwrap()..(*sel_coarse_chan_idxs.last().unwrap() + 1);

        let all_freqs_hz = context.get_fine_chan_freqs_hz_array(sel_coarse_chan_idxs);

        let array_pos = LatLngHeight::new_mwa();

        let phase_centre_ra = RADec::from_mwalib_phase_or_pointing(&context.metafits_context);
        // let lst_rad = context.metafits_context.lst_rad;
        // let phase_centre_ha = phase_centre_ra.to_hadec(lst_rad);
        let tiles_xyz_geod = XyzGeodetic::get_tiles_mwa(&context.metafits_context);

        let (jones_array, _) = context_to_jones_array(
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            // sel_baseline_idxs.as_slice(),
            None,
            false,
        )
        .unwrap();

        let jones_array = jones_array.mapv(TestJones::from);

        // ts 0, chan 0 (cc 0, fc 0), baseline 0
        let viz_0_0_0 = *jones_array.get((0, 0, 0)).unwrap();
        assert_abs_diff_eq!(
            viz_0_0_0,
            TestJones::from([
                Complex::new(0x410000 as f32, 0x410001 as f32),
                Complex::new(0x410002 as f32, 0x410003 as f32),
                Complex::new(0x410004 as f32, 0x410005 as f32),
                Complex::new(0x410006 as f32, 0x410007 as f32),
            ])
        );

        // ts 0, chan 0 (cc 0, fc 0), baseline 5
        let viz_0_0_1 = *jones_array.get((0, 0, 1)).unwrap();
        assert_abs_diff_eq!(
            viz_0_0_1,
            TestJones::from([
                Complex::new(0x410010 as f32, 0x410011 as f32),
                Complex::new(0x410012 as f32, 0x410013 as f32),
                Complex::new(0x410014 as f32, 0x410015 as f32),
                Complex::new(0x410016 as f32, 0x410017 as f32),
            ])
        );

        // ts 3, chan 3 (cc 1, fc 1), baseline 5
        let viz_3_3_1 = *jones_array.get((3, 3, 1)).unwrap();
        assert_abs_diff_eq!(
            viz_3_3_1,
            TestJones::from([
                Complex::new(0x410718 as f32, 0x410719 as f32),
                Complex::new(0x41071a as f32, 0x41071b as f32),
                Complex::new(0x41071c as f32, 0x41071d as f32),
                Complex::new(0x41071e as f32, 0x41071f as f32),
            ])
        );

        let integration_time_s = context.metafits_context.corr_int_time_ms as f64 / 1000.0;

        // timestep 0
        let timestep_0 = &context.timesteps[sel_timestep_idxs[0]];
        let epoch_0 = Epoch::from_gpst_seconds(
            timestep_0.gps_time_ms as f64 / 1000.0 + integration_time_s / 2.0,
        );
        let prec_info_0 = precess_time(
            phase_centre_ra,
            epoch_0,
            array_pos.longitude_rad,
            array_pos.latitude_rad,
        );
        let phase_centre_ha_j2000_0 = prec_info_0.hadec_j2000; // phase_centre_ra.to_hadec(prec_info_0.lmst_j2000);
        let tiles_xyz_precessed_0 = prec_info_0.precess_xyz_parallel(&tiles_xyz_geod);
        // timestep 3
        let timestep_3 = &context.timesteps[sel_timestep_idxs[3]];
        let epoch_3 = Epoch::from_gpst_seconds(
            timestep_3.gps_time_ms as f64 / 1000.0 + integration_time_s / 2.0,
        );
        let prec_info_3 = precess_time(
            phase_centre_ra,
            epoch_3,
            array_pos.longitude_rad,
            array_pos.latitude_rad,
        );
        let phase_centre_ha_j2000_3 = prec_info_3.hadec_j2000; // phase_centre_ra.to_hadec(prec_info_3.lmst_j2000);
        let tiles_xyz_precessed_3 = prec_info_3.precess_xyz_parallel(&tiles_xyz_geod);

        // baseline 1
        let bl_1 = &context.metafits_context.baselines[1];

        // ts 0, baseline 1
        let xyz_0_1 =
            tiles_xyz_precessed_0[bl_1.ant1_index] - tiles_xyz_precessed_0[bl_1.ant2_index];
        let w_0_1 = UVW::from_xyz(xyz_0_1, phase_centre_ha_j2000_0).w;
        // ts 3, baseline 1
        let xyz_3_1 =
            tiles_xyz_precessed_3[bl_1.ant1_index] - tiles_xyz_precessed_3[bl_1.ant2_index];
        let w_3_1 = UVW::from_xyz(xyz_3_1, phase_centre_ha_j2000_3).w;

        // ts 0, chan 0 (cc 0, fc 0), baseline 1
        let angle_0_0_1 = -2.0 * PI * w_0_1 * (all_freqs_hz[0] as f64) / VEL_C;
        let (sin_0_0_1_f64, cos_0_0_1_f64) = angle_0_0_1.sin_cos();
        let (sin_0_0_1, cos_0_0_1) = (sin_0_0_1_f64 as f32, cos_0_0_1_f64 as f32);

        // ts 3, chan 3 (cc 1, fc 1), baseline 1
        let angle_3_3_1 = -2.0 * PI * w_3_1 * (all_freqs_hz[3] as f64) / VEL_C;
        let (sin_3_3_1_f64, cos_3_3_1_f64) = angle_3_3_1.sin_cos();
        let (sin_3_3_1, cos_3_3_1) = (sin_3_3_1_f64 as f32, cos_3_3_1_f64 as f32);

        // let angle_5_3 = -2.0 * PI * w_5_3 * (all_freqs_hz[3] as f64) / VEL_C;
        let mut jones_array = jones_array.mapv(|e| Jones::from([e[0], e[1], e[2], e[3]]));
        correct_geometry(
            &context,
            &mut jones_array,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            None,
            None,
            false,
        );
        let jones_array = jones_array.mapv(TestJones::from);
        // there should be no difference in baseline 0
        // ts 0 (batch 0, scan 0), chan 0 (cc 0, fc 0), baseline 0
        assert_abs_diff_eq!(*jones_array.get((0, 0, 0)).unwrap(), viz_0_0_0);

        ////
        // baseline 1 should be rotated
        ////
        // ts 0 (batch 0, scan 0), chan 0 (cc 0, fc 0), baseline 1
        // -> pol XX
        let rot_1_xx_0_0_re = (cos_0_0_1 * viz_0_0_1[0].re - sin_0_0_1 * viz_0_0_1[0].im) as f32;
        let rot_1_xx_0_0_im = (sin_0_0_1 * viz_0_0_1[0].re + cos_0_0_1 * viz_0_0_1[0].im) as f32;
        // -> pol XY
        let rot_1_xy_0_0_re = (cos_0_0_1 * viz_0_0_1[1].re - sin_0_0_1 * viz_0_0_1[1].im) as f32;
        let rot_1_xy_0_0_im = (sin_0_0_1 * viz_0_0_1[1].re + cos_0_0_1 * viz_0_0_1[1].im) as f32;
        // -> pol YX
        let rot_1_yx_0_0_re = (cos_0_0_1 * viz_0_0_1[2].re - sin_0_0_1 * viz_0_0_1[2].im) as f32;
        let rot_1_yx_0_0_im = (sin_0_0_1 * viz_0_0_1[2].re + cos_0_0_1 * viz_0_0_1[2].im) as f32;
        // -> pol YY
        let rot_1_yy_0_0_re = (cos_0_0_1 * viz_0_0_1[3].re - sin_0_0_1 * viz_0_0_1[3].im) as f32;
        let rot_1_yy_0_0_im = (sin_0_0_1 * viz_0_0_1[3].re + cos_0_0_1 * viz_0_0_1[3].im) as f32;

        assert_abs_diff_eq!(
            *jones_array.get((0, 0, 1)).unwrap(),
            &TestJones::from([
                Complex::new(rot_1_xx_0_0_re, rot_1_xx_0_0_im),
                Complex::new(rot_1_xy_0_0_re, rot_1_xy_0_0_im),
                Complex::new(rot_1_yx_0_0_re, rot_1_yx_0_0_im),
                Complex::new(rot_1_yy_0_0_re, rot_1_yy_0_0_im),
            ])
        );

        // ts 3 (batch 1, scan 1), chan 3 (cc 1, fc 1), baseline 1
        // -> pol XX
        let rot_1_xx_3_3_re = (cos_3_3_1 * viz_3_3_1[0].re - sin_3_3_1 * viz_3_3_1[0].im) as f32;
        let rot_1_xx_3_3_im = (sin_3_3_1 * viz_3_3_1[0].re + cos_3_3_1 * viz_3_3_1[0].im) as f32;
        // -> pol XY
        let rot_1_xy_3_3_re = (cos_3_3_1 * viz_3_3_1[1].re - sin_3_3_1 * viz_3_3_1[1].im) as f32;
        let rot_1_xy_3_3_im = (sin_3_3_1 * viz_3_3_1[1].re + cos_3_3_1 * viz_3_3_1[1].im) as f32;
        // -> pol YX
        let rot_1_yx_3_3_re = (cos_3_3_1 * viz_3_3_1[2].re - sin_3_3_1 * viz_3_3_1[2].im) as f32;
        let rot_1_yx_3_3_im = (sin_3_3_1 * viz_3_3_1[2].re + cos_3_3_1 * viz_3_3_1[2].im) as f32;
        // -> pol YY
        let rot_1_yy_3_3_re = (cos_3_3_1 * viz_3_3_1[3].re - sin_3_3_1 * viz_3_3_1[3].im) as f32;
        let rot_1_yy_3_3_im = (sin_3_3_1 * viz_3_3_1[3].re + cos_3_3_1 * viz_3_3_1[3].im) as f32;

        assert_abs_diff_eq!(
            *jones_array.get((3, 3, 1)).unwrap(),
            &TestJones::from([
                Complex::new(rot_1_xx_3_3_re, rot_1_xx_3_3_im),
                Complex::new(rot_1_xy_3_3_re, rot_1_xy_3_3_im),
                Complex::new(rot_1_yx_3_3_re, rot_1_yx_3_3_im),
                Complex::new(rot_1_yy_3_3_re, rot_1_yy_3_3_im),
            ])
        );
    }

    macro_rules! compare_jones {
        ($a:expr, $b:expr) => {
            assert_abs_diff_eq!(TestJones::<f32>::from($a), TestJones::<f32>::from($b));
        };
    }

    #[test]
    fn test_correct_digital_gains() {
        let context = get_mwa_ord_context();

        let sel_coarse_chan_range = get_coarse_chan_range(&context).unwrap();
        let sel_timestep_range = get_timestep_range(&context).unwrap();

        let (jones_array, _) = context_to_jones_array(
            &context,
            &sel_timestep_range,
            &sel_coarse_chan_range,
            // sel_baseline_idxs.as_slice(),
            None,
            false,
        )
        .unwrap();

        let sel_baseline_range = 0..2;
        let sel_baseline_idxs: Vec<usize> = sel_baseline_range.clone().collect();

        let sel_baselines = sel_baseline_idxs
            .iter()
            .map(|&idx| {
                let baseline = &context.metafits_context.baselines[idx];
                (baseline.ant1_index, baseline.ant2_index)
            })
            .collect::<Vec<_>>();

        let first_coarse_chan_idx = sel_coarse_chan_range.start;
        let last_coarse_chan_idx = sel_coarse_chan_range.end - 1;
        // first coarse chan, antenna 0
        let gain_first_0 =
            context.metafits_context.antennas[0].rfinput_x.digital_gains[first_coarse_chan_idx];
        // first coarse chan, antenna 1
        let gain_first_1 =
            context.metafits_context.antennas[1].rfinput_x.digital_gains[first_coarse_chan_idx];
        // last coarse chan, antenna 0
        let gain_last_0 =
            context.metafits_context.antennas[0].rfinput_x.digital_gains[last_coarse_chan_idx];
        // last coarse chan, antenna 1
        let gain_last_1 =
            context.metafits_context.antennas[1].rfinput_x.digital_gains[last_coarse_chan_idx];

        let max_chan = jones_array.dim().1 - 1;

        // ts 0, first chan, baseline 0 (0 vs 0)
        let jones_0_min_0 = jones_array[(0, 0, 0)];
        // ts 0, first chan, baseline 1 (0 vs 1)
        let jones_0_min_1 = jones_array[(0, 0, 1)];
        // ts 0, last chan, baseline 0 (0 vs 0)
        let jones_0_max_0 = jones_array[(0, max_chan, 0)];
        // ts 0, last chan, baseline 1 (0 vs 1)
        let jones_0_max_1 = jones_array[(0, max_chan, 1)];

        let mut out_jones_array = jones_array.slice(s![.., .., sel_baseline_range]).to_owned();

        correct_digital_gains(
            &context,
            &mut out_jones_array,
            &sel_coarse_chan_range,
            &sel_baselines,
        )
        .unwrap();

        compare_jones!(
            out_jones_array[(0, 0, 0)],
            jones_0_min_0 / ((gain_first_0 * gain_first_0) as f32)
        );
        compare_jones!(
            out_jones_array[(0, 0, 1)],
            jones_0_min_1 / ((gain_first_0 * gain_first_1) as f32)
        );
        compare_jones!(
            out_jones_array[(0, max_chan, 0)],
            jones_0_max_0 / ((gain_last_0 * gain_last_0) as f32)
        );
        compare_jones!(
            out_jones_array[(0, max_chan, 1)],
            jones_0_max_1 / ((gain_last_0 * gain_last_1) as f32)
        );
    }

    #[test]
    fn test_correct_digital_gains_bad_array_shape() {
        let context = get_mwa_ord_context();

        let sel_coarse_chan_range = get_coarse_chan_range(&context).unwrap();

        let mut jones_array = Array3::from_shape_fn((1, 1, 1), |_| Jones::nan());

        let sel_baseline_range = 0..2;
        let sel_baseline_idxs: Vec<usize> = sel_baseline_range.collect();

        let sel_baselines = sel_baseline_idxs
            .iter()
            .map(|&idx| {
                let baseline = &context.metafits_context.baselines[idx];
                (baseline.ant1_index, baseline.ant2_index)
            })
            .collect::<Vec<_>>();

        assert!(matches!(
            correct_digital_gains(
                &context,
                &mut jones_array,
                &sel_coarse_chan_range,
                &sel_baselines,
            ),
            Err(BirliError::BadArrayShape { .. })
        ));

        let channels =
            sel_coarse_chan_range.len() * context.metafits_context.num_corr_fine_chans_per_coarse;

        let mut jones_array = Array3::from_shape_fn((1, channels, 1), |_| Jones::nan());

        assert!(matches!(
            correct_digital_gains(
                &context,
                &mut jones_array,
                &sel_coarse_chan_range,
                &sel_baselines,
            ),
            Err(BirliError::BadArrayShape { .. })
        ));

        let gains = Array2::from_shape_fn((1, 1), |_| (1., 1.));

        assert!(matches!(
            _correct_digital_gains(&mut jones_array, &gains, &sel_baselines, 2),
            Err(BirliError::BadArrayShape { .. })
        ));
    }

    #[test]
    fn test_correct_pfb_gains() {
        let num_fine_chans_per_coarse = 2;
        // 2 coarse channels
        let mut jones_array: Array3<Jones<f32>> =
            Array3::from_shape_fn((2, 4, 2), |_| Jones::identity());
        let coarse_band_gains = vec![0.1, 0.2];

        correct_pfb_gains(
            &mut jones_array,
            &coarse_band_gains,
            num_fine_chans_per_coarse,
        )
        .unwrap();

        let expected: [Jones<f32>; 4] = [
            Jones::identity() / 0.1,
            Jones::identity() / 0.2,
            Jones::identity() / 0.1,
            Jones::identity() / 0.2,
        ];

        for (&expected, actuals) in izip!(expected.iter(), jones_array.axis_iter(Axis(1))) {
            for actual in actuals.iter() {
                compare_jones!(*expected, *actual);
            }
        }
    }
}
