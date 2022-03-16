//! Corrections that can be performed on visibility data
use crate::{
    ndarray::{parallel::prelude::*, Array2, Array3, Axis},
    BirliError, Jones,
};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use itertools::izip;
use log::trace;
use marlu::{
    constants::VEL_C,
    hifitime::Epoch,
    mwalib::{CorrelatorContext, MWAVersion},
    precession::precess_time,
    Complex, LatLngHeight, RADec, XyzGeodetic, UVW,
};
use std::{f64::consts::PI, ops::Range};

/// Perform cable length corrections, given an observation's
/// [`marlu::mwalib::CorrelatorContext`] and an [`ndarray::Array3`] of [`crate::TestJones`]
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
/// use birli::{correct_cable_lengths, mwalib::CorrelatorContext, VisSelection};
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
/// correct_cable_lengths(&corr_ctx, &mut jones_array, &vis_sel.coarse_chan_range, false);
/// ```
///
/// # Accuracy
///
/// This follows the Cotter implementation of cable length correction, however
/// there is a slower but more accurate version of the calculation which
/// uses f64 values for the `sin_cos`. According to benchmarks, the Cotter
/// implementation is about 32% faster (5.9 seconds vs 8.6) than the more
/// precise implementation, and they vary by about three parts in four million.
/// Therefore it was decided that the Cotter implementation was more favourable.
///
/// # Assumptions
///
/// - an Antenna's rfinput is the same for X and Y
pub fn correct_cable_lengths(
    corr_ctx: &CorrelatorContext,
    jones_array: &mut Array3<Jones<f32>>,
    coarse_chan_range: &Range<usize>,
    // TODO: allow subset of baselines
    // baseline_idxs: &[usize],
    draw_progress: bool,
) {
    trace!("start correct_cable_lengths");

    let baselines = &corr_ctx.metafits_context.baselines;
    let antennas = &corr_ctx.metafits_context.antennas;

    let all_freqs_hz =
        corr_ctx.get_fine_chan_freqs_hz_array(&coarse_chan_range.clone().collect::<Vec<_>>());
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

    // TODO: his whole thing is pretty inelegant, needs a re-do.
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
/// [`crate::mwalib::CorrelatorContext`] and an [`ndarray::Array3`] of [`crate::TestJones`]
/// visibilities
///
/// Complex visibilities are phase-shifted by an angle determined by the length
/// of the w-coordinate for the baseline and the channel's frequency.
///
/// # Examples
///
/// ```rust
/// use birli::{
///     FlagContext, correct_geometry, mwalib::CorrelatorContext, VisSelection,
///     correct_cable_lengths
/// };
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
/// // Determine which timesteps and coarse channels we want to use
/// let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
///
/// let sel_timestep_idxs = &corr_ctx.common_timestep_indices;
/// vis_sel.timestep_range =
///     *sel_timestep_idxs.first().unwrap()..(*sel_timestep_idxs.last().unwrap() + 1);
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
/// correct_cable_lengths(&corr_ctx, &mut jones_array, &vis_sel.coarse_chan_range, false);
///
/// correct_geometry(
///     &corr_ctx,
///     &mut jones_array,
///     &vis_sel.timestep_range,
///     &vis_sel.coarse_chan_range,
///     None,
///     None,
///     false,
/// );
/// ```
#[allow(clippy::too_many_arguments)]
pub fn correct_geometry(
    corr_ctx: &CorrelatorContext,
    jones_array: &mut Array3<Jones<f32>>,
    timestep_range: &Range<usize>,
    coarse_chan_range: &Range<usize>,
    // TODO: allow subset of baselines
    // baseline_idxs: &[usize],
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

    let timesteps = &corr_ctx.timesteps[timestep_range.clone()];

    let baselines = &corr_ctx.metafits_context.baselines;

    let all_freqs_hz =
        corr_ctx.get_fine_chan_freqs_hz_array(&coarse_chan_range.clone().collect::<Vec<_>>());
    let jones_dims = jones_array.dim();

    let integration_time_s = corr_ctx.metafits_context.corr_int_time_ms as f64 / 1000.0;

    let phase_centre_ra = match phase_centre_ra {
        Some(pc) => pc,
        None => RADec::from_mwalib_phase_or_pointing(&corr_ctx.metafits_context),
    };
    let tiles_xyz_geod = XyzGeodetic::get_tiles(&corr_ctx.metafits_context, array_pos.latitude_rad);

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
/// `coarse_chan_range`, such that `jones_array.dim().1 = num_fine_chans_per_coarse * coarse_chan_range.len()`
///
/// # Arguments
///
/// - `corr_ctx` - The correlator [`marlu::mwalib::CorrelatorContext`].
/// - `jones_array` - The array of Jones matrices to be corrected, [timestep][channel][baseleine].
/// - `coarse_chan_range` - The range of mwalib coarse channels which are used in in the channel
///     dimension of the jones array.
///
/// # Assumptions
/// - the digital gains are provided in [`marlu::mwalib::Rfinput.digital_gains`] in the same order as the
///   coarse channel indices (increasing sky frequency)
///
/// # Errors
/// - Will throw [`BirliError::BadArrayShape`] if:
///     - `jones_array.dim().1 != num_fine_chans_per_coarse * coarse_chan_range.len()`
///     - `jones_array.dim().2 != ant_pairs.len()`
pub fn correct_digital_gains(
    corr_ctx: &CorrelatorContext,
    jones_array: &mut Array3<Jones<f32>>,
    coarse_chan_range: &Range<usize>,
    ant_pairs: &[(usize, usize)],
    // TODO: take a VisSelection
) -> Result<(), BirliError> {
    let num_fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;

    let vis_dims = jones_array.dim();
    if vis_dims.1 != coarse_chan_range.len() * num_fine_chans_per_coarse {
        return Err(BirliError::BadArrayShape {
            argument: "coarse_chan_range".into(),
            function: "correct_digital_gains".into(),
            expected: format!(
                "(vis_dims.1={}) / (num_fine_chans_per_coarse={}) = {}",
                vis_dims.1,
                num_fine_chans_per_coarse,
                vis_dims.1 / num_fine_chans_per_coarse
            ),
            received: format!("{:?}", coarse_chan_range.len()),
        });
    }
    if vis_dims.2 != ant_pairs.len() {
        return Err(BirliError::BadArrayShape {
            argument: "ant_pairs".into(),
            function: "correct_digital_gains".into(),
            expected: format!("vis_dims.2={}", vis_dims.2,),
            received: format!("{:?}", ant_pairs.len()),
        });
    }
    assert!(vis_dims.2 == ant_pairs.len());

    let gains = Array2::from_shape_fn(
        (corr_ctx.metafits_context.num_ants, coarse_chan_range.len()),
        |(ant_idx, coarse_chan_idx)| {
            let ant = corr_ctx.metafits_context.antennas[ant_idx].clone();
            (
                ant.rfinput_x.digital_gains[coarse_chan_range.clone()][coarse_chan_idx],
                ant.rfinput_y.digital_gains[coarse_chan_range.clone()][coarse_chan_idx],
            )
        },
    );

    _correct_digital_gains(jones_array, &gains, ant_pairs, num_fine_chans_per_coarse)
}

fn _correct_digital_gains(
    jones_array: &mut Array3<Jones<f32>>,
    gains: &Array2<(f64, f64)>,
    ant_pairs: &[(usize, usize)],
    num_fine_chans_per_coarse: usize,
) -> Result<(), BirliError> {
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
        ant_pairs.iter().copied()
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

/// Correct for coarse pfb bandpass shape in each coarse channel by scaling `passband_gains` to
/// fit the coarse band.
///
/// # Arguments
///
/// - `jones_array` - The array of Jones matrices to be corrected, [timestep][channel][baseleine].
/// - `weight_array` - The array of weights to be corrected, same dimensions are `jones_array`.
/// - `passband_gains` - a slice of gains to be applied to each coarse channel.
/// - `num_fine_chans_per_coarse` - The number of fine channels in each coarse.
/// - `mwax` - Whether to emulate the MWAX correlator (true) or legacy (false)
///
/// # Errors
///
/// Will throw `BirliError::BadArrayShape` if:
/// - `num_fine_chans_per_coarse` is zero
/// - The length of the channel axis in `jones_array` is not a multiple of `num_fine_chans_per_coarse`.
/// - `jones_array` and `weight_array` have different shapes.
/// - The length of the coarse band gains is not a multiple of `num_fine_chans_per_coarse`, or vice versa.
///
pub fn correct_coarse_passband_gains(
    jones_array: &mut Array3<Jones<f32>>,
    weight_array: &mut Array3<f32>,
    passband_gains: &[f64],
    num_fine_chans_per_coarse: usize,
    scrunch_type: &ScrunchType,
) -> Result<(), BirliError> {
    if num_fine_chans_per_coarse == 0 {
        return Err(BirliError::BadArrayShape {
            argument: "num_fine_chans_per_coarse".into(),
            function: "correct_coarse_passband_gains".into(),
            expected: "a number greater than zero".into(),
            received: format!("{:?}", num_fine_chans_per_coarse),
        });
    }

    if jones_array.dim().1 % num_fine_chans_per_coarse != 0 {
        return Err(BirliError::BadArrayShape {
            argument: "jones_array".into(),
            function: "correct_coarse_passband_gains".into(),
            expected: format!(
                "(_, n, _), where n is a multiple of num_fine_chans_per_coarse={}",
                num_fine_chans_per_coarse
            ),
            received: format!("{:?}", jones_array.dim()),
        });
    };

    if weight_array.dim() != jones_array.dim() {
        return Err(BirliError::BadArrayShape {
            argument: "weight_array".into(),
            function: "correct_coarse_passband_gains".into(),
            expected: format!("same as jones_array.dim()={:?}", jones_array.dim()),
            received: format!("{:?}", weight_array.dim()),
        });
    };

    let fscrunch = if passband_gains.len() % num_fine_chans_per_coarse == 0 {
        passband_gains.len() / num_fine_chans_per_coarse
    } else {
        return Err(BirliError::BadArrayShape {
            argument: "passband_gains".into(),
            function: "correct_coarse_passband_gains".into(),
            expected: format!(
                "n, where n is a multiple of num_fine_chans_per_coarse={}",
                num_fine_chans_per_coarse
            ),
            received: format!("{:?}", passband_gains.len()),
        });
    };

    let scrunched_gains = scrunch_gains(passband_gains, fscrunch, scrunch_type);

    for (mut jones_array, mut weight_array) in izip!(
        jones_array.axis_chunks_iter_mut(Axis(1), num_fine_chans_per_coarse),
        weight_array.axis_chunks_iter_mut(Axis(1), num_fine_chans_per_coarse),
    ) {
        for (mut jones_array, mut weight_array, &gain) in izip!(
            jones_array.axis_iter_mut(Axis(1)),
            weight_array.axis_iter_mut(Axis(1)),
            scrunched_gains.iter()
        ) {
            for (jones, weight) in izip!(jones_array.iter_mut(), weight_array.iter_mut()) {
                let jones_f64 = Jones::<f64>::from(*jones);
                *jones = Jones::<f32>::from(jones_f64 / gain);
                *weight = (gain * (*weight as f64)) as f32;
            }
        }
    }

    Ok(())
}

#[derive(Clone)]
/// Possible types of scrunching that a correlator might do.
pub enum ScrunchType {
    /// used by legacy correlator
    Simple,
    /// used by MWAX
    CenterSymmetric,
}

impl ScrunchType {
    /// Return the corresponding scrunch version from the [`marlu::mwalib::MWAVersion`].
    ///
    /// # Errors
    ///
    /// Will throw [`BirliError::BadMWAVersion`] If you provide something other than
    /// [`marlu::mwalib::MWAVersion::CorrMWAXv2`], [`marlu::mwalib::MWAVersion::CorrLegacy`], or
    /// [`marlu::mwalib::MWAVersion::CorrOldLegacy`]
    pub fn from_mwa_version(ver: MWAVersion) -> Result<Self, BirliError> {
        match ver {
            MWAVersion::CorrMWAXv2 => Ok(Self::CenterSymmetric),
            MWAVersion::CorrLegacy | MWAVersion::CorrOldLegacy => Ok(Self::Simple),
            ver => Err(BirliError::BadMWAVersion {
                message:
                    "could not determine correlator scrunch type from provided mwalib::MWAVersion."
                        .into(),
                version: ver.to_string(),
            }),
        }
    }
}

/// Given a set of gains for each ultrafine channel, compute the gains for each fine channel,
/// is if these were averaged by given correlator type.
///
/// # Averaging
///
/// a scrunched channel contains fscrunch times the bandwidth of the ultrafine channel.
///
/// MWAX example with 12 ultrafine channels:
/// fscrunch|0| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |10 |11 |0|
/// --------|-|---|---|---|---|---|---|---|---|---|---|---|-|
/// 2 (e->e)| 0*|   1   |   2   |   3   |   4   |   5   | 0*|
/// 3 (o->e)|  0* |     1     |     2     |     3     |  0* |
/// 4 (e->o)|       0       |       1       |       2       |
///
/// MWAX example with 15 ultrafine channels:
/// fscrunch|0| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |10 |11 |12 |13 |14 |0|
/// --------|-|---|---|---|---|---|---|---|---|---|---|---|---|---|---| |
/// 3 (o->o)|     0*    |     1     |     2     |     3     |     4     |
/// 5 (o->o)|           0       |         1         |         2         |
///
/// For more details see: <https://wiki.mwatelescope.org/display/MP/MWA+Fine+Channel+Centre+Frequencies>
pub fn scrunch_gains(
    ultrafine_gains: &[f64],
    fscrunch: usize,
    scrunch_type: &ScrunchType,
) -> Vec<f64> {
    let scrunched_length = ultrafine_gains.len() / fscrunch;
    if fscrunch == 1 {
        ultrafine_gains.to_vec()
    } else {
        // in mwax, for an even number of scrunched channels, a given scrunched channel `n`:
        // - is centered on ultrafine channel `n * fscrunch`.
        // - contains power from the center channel, and some power from channels within
        //   `(fscrunch - 1) / 2` channels of the centre, wrapping around the coarse channel
        //
        // Scrunched channel zero will contain an equal amount power from the low and high ends of
        // the coarse channel.

        // in mwax, for an odd number of scrunched channels, a given scrunched channel `n`:
        // - contains half the power of ultrafine channel `n * fscrunch`.
        // - contains all of the power from ultrafine channels `n * fscrunch + (1..n-1)`
        // - contains half the power from ultrafine channel `n * fscrunch + n`
        // - centered on ultrafine channel `(n + 1/2) * fscrunch`.s
        //
        // Scrunched channel zero will contain an equal amount power from the low and high ends of
        // the coarse channel.
        // Here we use a cycle iterator, and calculate the first ultrafine channel which contributes
        // to scrunched channel zero, `first_zeroth_channel`.

        // get the relative index and weight of the window of ultrafine channels that contribute to a scrunched channel.
        #[rustfmt::skip]
        let window_offset_weights: Vec<(i32, f64)> = match (&scrunch_type, scrunched_length % 2, fscrunch % 2) {
            (ScrunchType::Simple, _, _) => (0..fscrunch).map(|w| (w as i32, 1./fscrunch as f64)).collect(),
            // even channels, even fscrunch: window length is fscrunch + 1, half-weighted edges
            (ScrunchType::CenterSymmetric, 0, 0) => (0..=fscrunch)
                .map(|w| (
                    (w as i32 - fscrunch as i32 / 2),
                    (if w == 0 || w == fscrunch { 0.5 } else { 1. }) / fscrunch as f64,
                ))
                .collect(),
            // even channels, odd fscrunch: window length is fscrunch, equal weights
            (ScrunchType::CenterSymmetric, 0, 1) => (0..fscrunch)
                .map(|w| (
                    (w as i32 - (fscrunch as i32 - 1) / 2),
                    1. / fscrunch as f64
                ))
                .collect(),
            // odd channels: window length is fscrunch + 1, half-weighted edges
            (ScrunchType::CenterSymmetric, 1, _) => (0..=fscrunch)
                .map(|w| (
                    (w as i32),
                    (if w == 0 || w == fscrunch { 0.5 } else { 1. }) / fscrunch as f64,
                ))
                .collect(),
            _ => unreachable!(),
        };
        // apply the weights to calculate the scrunched gains
        (0..scrunched_length)
            .map(|scrunched_chan| {
                window_offset_weights
                    .iter()
                    .fold(0., |acc, &(offset, weight)| {
                        // rem_euclid is basically mod but it correctly wraps around with negative numbers.
                        let ultrafine_chan = (((fscrunch * scrunched_chan) as i32 + offset)
                            .rem_euclid(ultrafine_gains.len() as i32))
                            as usize;
                        acc + ultrafine_gains[ultrafine_chan] * weight
                    })
            })
            .collect()
    }
}

#[cfg(test)]
#[allow(clippy::similar_names)]
mod tests {

    use super::{
        _correct_digital_gains, correct_cable_lengths, correct_coarse_passband_gains,
        correct_digital_gains, correct_geometry, scrunch_gains, VEL_C,
    };
    use float_cmp::assert_approx_eq;
    use itertools::izip;
    use marlu::{
        hifitime::Epoch, precession::precess_time, Complex, Jones, LatLngHeight, RADec,
        XyzGeodetic, UVW,
    };
    use ndarray::{s, Array2, Array3, Axis};
    use std::f64::consts::PI;

    use crate::{
        approx::assert_abs_diff_eq,
        compare_jones,
        corrections::ScrunchType,
        test_common::{get_mwa_ord_context, get_mwax_context},
        BirliError, TestJones, VisSelection,
    };

    #[test]
    fn test_cable_length_corrections_mwax() {
        let corr_ctx = get_mwax_context();
        let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

        // Create a blank array to store flags and visibilities
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        // read visibilities out of the gpubox files
        vis_sel
            .read_mwalib(
                &corr_ctx,
                jones_array.view_mut(),
                flag_array.view_mut(),
                false,
            )
            .unwrap();

        let coarse_chan_indices: Vec<_> = vis_sel.coarse_chan_range.clone().collect();
        let all_freqs_hz = corr_ctx.get_fine_chan_freqs_hz_array(&coarse_chan_indices);

        // ts 0, chan 0 (cc 0, fc 0), baseline 0
        let viz_0_0_0 = jones_array[(0, 0, 0)];
        compare_jones!(
            viz_0_0_0,
            Jones::from([
                Complex::new(0x410000 as f32, 0x410001 as f32),
                Complex::new(0x410002 as f32, 0x410003 as f32),
                Complex::new(0x410004 as f32, 0x410005 as f32),
                Complex::new(0x410006 as f32, 0x410007 as f32),
            ])
        );

        // ts 0, chan 0 (cc 0, fc 0), baseline 1
        let viz_0_0_1 = jones_array[(0, 0, 1)];
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
        let viz_3_3_1 = jones_array[(3, 3, 1)];
        compare_jones!(
            viz_3_3_1,
            Jones::from([
                Complex::new(0x410718 as f32, 0x410719 as f32),
                Complex::new(0x41071a as f32, 0x41071b as f32),
                Complex::new(0x41071c as f32, 0x41071d as f32),
                Complex::new(0x41071e as f32, 0x41071f as f32),
            ])
        );

        // baseline 1, input 1, pol x
        let length_1_1_x = &corr_ctx.metafits_context.antennas
            [corr_ctx.metafits_context.baselines[1].ant1_index]
            .rfinput_x
            .electrical_length_m;
        // baseline 1, input 1, pol y
        let length_1_1_y = &corr_ctx.metafits_context.antennas
            [corr_ctx.metafits_context.baselines[1].ant1_index]
            .rfinput_y
            .electrical_length_m;
        // baseline 1, input 2, pol x
        let length_1_2_x = &corr_ctx.metafits_context.antennas
            [corr_ctx.metafits_context.baselines[1].ant2_index]
            .rfinput_x
            .electrical_length_m;
        // baseline 1, input 2, pol y
        let length_1_2_y = &corr_ctx.metafits_context.antennas
            [corr_ctx.metafits_context.baselines[1].ant2_index]
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

        correct_cable_lengths(
            &corr_ctx,
            &mut jones_array,
            &vis_sel.coarse_chan_range,
            false,
        );

        // there should be no difference in baseline 0
        // ts 0, chan 0, baseline 0
        compare_jones!(jones_array[(0, 0, 0)], viz_0_0_0);

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

        compare_jones!(
            jones_array[(0, 0, 1)],
            Jones::from([
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

        compare_jones!(
            jones_array[(3, 3, 1)],
            Jones::from([
                Complex::new(rot_1_xx_3_3_re, rot_1_xx_3_3_im),
                Complex::new(rot_1_xy_3_3_re, rot_1_xy_3_3_im),
                Complex::new(rot_1_yx_3_3_re, rot_1_yx_3_3_im),
                Complex::new(rot_1_yy_3_3_re, rot_1_yy_3_3_im),
            ])
        );
    }

    #[test]
    fn test_cable_length_corrections_ord() {
        let corr_ctx = get_mwa_ord_context();
        let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

        // Create a blank array to store flags and visibilities
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        // read visibilities out of the gpubox files
        vis_sel
            .read_mwalib(
                &corr_ctx,
                jones_array.view_mut(),
                flag_array.view_mut(),
                false,
            )
            .unwrap();

        let coarse_chan_indices: Vec<_> = vis_sel.coarse_chan_range.clone().collect();
        let all_freqs_hz = corr_ctx.get_fine_chan_freqs_hz_array(&coarse_chan_indices);

        // ts 0, chan 0 (cc 0, fc 0), baseline 0
        let viz_0_0_0 = jones_array[(0, 0, 0)];
        compare_jones!(
            viz_0_0_0,
            Jones::from([
                Complex::new(0x10c5be as f32, -0x10c5bf as f32),
                Complex::new(0x10c5ae as f32, 0x10c5af as f32),
                Complex::new(0x10c5ae as f32, -0x10c5af as f32),
                Complex::new(0x10bec6 as f32, -0x10bec7 as f32),
            ])
        );

        // ts 0, chan 0 (cc 0, fc 0), baseline 5
        let viz_0_0_5 = jones_array[(0, 0, 5)];
        compare_jones!(
            viz_0_0_5,
            Jones::from([
                Complex::new(0x10f1ce as f32, -0x10f1cf as f32),
                Complex::new(0x10ea26 as f32, -0x10ea27 as f32),
                Complex::new(0x10f1be as f32, -0x10f1bf as f32),
                Complex::new(0x10ea16 as f32, -0x10ea17 as f32),
            ])
        );

        // ts 3, chan 3 (cc 1, fc 1), baseline 5
        let viz_3_3_5 = jones_array[(3, 3, 5)];
        compare_jones!(
            viz_3_3_5,
            Jones::from([
                Complex::new(0x0df3ce as f32, -0x0df3cf as f32),
                Complex::new(0x0dec26 as f32, -0x0dec27 as f32),
                Complex::new(0x0df3be as f32, -0x0df3bf as f32),
                Complex::new(0x0dec16 as f32, -0x0dec17 as f32),
            ])
        );

        // baseline 5, input 1, pol x
        let length_5_1_x = &corr_ctx.metafits_context.antennas
            [corr_ctx.metafits_context.baselines[5].ant1_index]
            .rfinput_x
            .electrical_length_m;
        // baseline 5, input 1, pol y
        let length_5_1_y = &corr_ctx.metafits_context.antennas
            [corr_ctx.metafits_context.baselines[5].ant1_index]
            .rfinput_y
            .electrical_length_m;
        // baseline 5, input 2, pol x
        let length_5_2_x = &corr_ctx.metafits_context.antennas
            [corr_ctx.metafits_context.baselines[5].ant2_index]
            .rfinput_x
            .electrical_length_m;
        // baseline 5, input 2, pol y
        let length_5_2_y = &corr_ctx.metafits_context.antennas
            [corr_ctx.metafits_context.baselines[5].ant2_index]
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

        correct_cable_lengths(
            &corr_ctx,
            &mut jones_array,
            &vis_sel.coarse_chan_range,
            false,
        );

        // there should be no difference in baseline 0
        // ts 0 (batch 0, scan 0), chan 0 (cc 0, fc 0), baseline 0
        compare_jones!(jones_array[(0, 0, 0)], viz_0_0_0);

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

        compare_jones!(
            jones_array[(0, 0, 5)],
            Jones::from([
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

        compare_jones!(
            jones_array[(3, 3, 5)],
            Jones::from([
                Complex::new(rot_5_xx_3_3_re, rot_5_xx_3_3_im),
                Complex::new(rot_5_xy_3_3_re, rot_5_xy_3_3_im),
                Complex::new(rot_5_yx_3_3_re, rot_5_yx_3_3_im),
                Complex::new(rot_5_yy_3_3_re, rot_5_yy_3_3_im),
            ])
        );
    }

    #[test]
    fn test_geometric_corrections_ord() {
        let corr_ctx = get_mwa_ord_context();
        let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        let coarse_chan_indices: Vec<_> = vis_sel.coarse_chan_range.clone().collect();
        let all_freqs_hz = corr_ctx.get_fine_chan_freqs_hz_array(&coarse_chan_indices);

        let array_pos = LatLngHeight::new_mwa();

        let phase_centre_ra = RADec::from_mwalib_phase_or_pointing(&corr_ctx.metafits_context);
        let tiles_xyz_geod = XyzGeodetic::get_tiles_mwa(&corr_ctx.metafits_context);

        // Create a blank array to store flags and visibilities
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        // read visibilities out of the gpubox files
        vis_sel
            .read_mwalib(
                &corr_ctx,
                jones_array.view_mut(),
                flag_array.view_mut(),
                false,
            )
            .unwrap();

        // ts 0, chan 0 (cc 0, fc 0), baseline 0
        let viz_0_0_0 = jones_array[(0, 0, 0)];
        compare_jones!(
            viz_0_0_0,
            Jones::from([
                Complex::new(0x10c5be as f32, -0x10c5bf as f32),
                Complex::new(0x10c5ae as f32, 0x10c5af as f32),
                Complex::new(0x10c5ae as f32, -0x10c5af as f32),
                Complex::new(0x10bec6 as f32, -0x10bec7 as f32),
            ])
        );

        // ts 0, chan 0 (cc 0, fc 0), baseline 5
        let viz_0_0_5 = jones_array[(0, 0, 5)];
        compare_jones!(
            viz_0_0_5,
            Jones::from([
                Complex::new(0x10f1ce as f32, -0x10f1cf as f32),
                Complex::new(0x10ea26 as f32, -0x10ea27 as f32),
                Complex::new(0x10f1be as f32, -0x10f1bf as f32),
                Complex::new(0x10ea16 as f32, -0x10ea17 as f32),
            ])
        );

        // ts 3, chan 3 (cc 1, fc 1), baseline 5
        let viz_3_3_5 = jones_array[(3, 3, 5)];
        compare_jones!(
            viz_3_3_5,
            Jones::from([
                Complex::new(0x0df3ce as f32, -0x0df3cf as f32),
                Complex::new(0x0dec26 as f32, -0x0dec27 as f32),
                Complex::new(0x0df3be as f32, -0x0df3bf as f32),
                Complex::new(0x0dec16 as f32, -0x0dec17 as f32),
            ])
        );

        let integration_time_s = corr_ctx.metafits_context.corr_int_time_ms as f64 / 1000.0;

        // timestep 0
        let timestep_0 = &corr_ctx.timesteps[vis_sel.timestep_range.clone()][0];
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
        let timestep_3 = &corr_ctx.timesteps[vis_sel.timestep_range.clone()][3];
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
        let bl_5 = &corr_ctx.metafits_context.baselines[5];

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

        correct_geometry(
            &corr_ctx,
            &mut jones_array,
            &vis_sel.timestep_range,
            &vis_sel.coarse_chan_range,
            None,
            None,
            false,
        );

        // there should be no difference in baseline 0
        // ts 0 (batch 0, scan 0), chan 0 (cc 0, fc 0), baseline 0
        compare_jones!(jones_array[(0, 0, 0)], viz_0_0_0);

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

        compare_jones!(
            jones_array[(0, 0, 5)],
            Jones::from([
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

        compare_jones!(
            jones_array[(3, 3, 5)],
            Jones::from([
                Complex::new(rot_5_xx_3_3_re, rot_5_xx_3_3_im),
                Complex::new(rot_5_xy_3_3_re, rot_5_xy_3_3_im),
                Complex::new(rot_5_yx_3_3_re, rot_5_yx_3_3_im),
                Complex::new(rot_5_yy_3_3_re, rot_5_yy_3_3_im),
            ])
        );
    }

    #[test]
    fn test_geometric_corrections_mwax() {
        let corr_ctx = get_mwax_context();
        let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

        // Create a blank array to store flags and visibilities
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        // read visibilities out of the gpubox files
        vis_sel
            .read_mwalib(
                &corr_ctx,
                jones_array.view_mut(),
                flag_array.view_mut(),
                false,
            )
            .unwrap();

        let coarse_chan_indices: Vec<_> = vis_sel.coarse_chan_range.clone().collect();
        let all_freqs_hz = corr_ctx.get_fine_chan_freqs_hz_array(&coarse_chan_indices);

        let array_pos = LatLngHeight::new_mwa();

        let phase_centre_ra = RADec::from_mwalib_phase_or_pointing(&corr_ctx.metafits_context);
        // let lst_rad = corr_ctx.metafits_context.lst_rad;
        // let phase_centre_ha = phase_centre_ra.to_hadec(lst_rad);
        let tiles_xyz_geod = XyzGeodetic::get_tiles_mwa(&corr_ctx.metafits_context);

        // ts 0, chan 0 (cc 0, fc 0), baseline 0
        let viz_0_0_0 = jones_array[(0, 0, 0)];
        compare_jones!(
            viz_0_0_0,
            Jones::from([
                Complex::new(0x410000 as f32, 0x410001 as f32),
                Complex::new(0x410002 as f32, 0x410003 as f32),
                Complex::new(0x410004 as f32, 0x410005 as f32),
                Complex::new(0x410006 as f32, 0x410007 as f32),
            ])
        );

        // ts 0, chan 0 (cc 0, fc 0), baseline 5
        let viz_0_0_1 = jones_array[(0, 0, 1)];
        compare_jones!(
            viz_0_0_1,
            Jones::from([
                Complex::new(0x410010 as f32, 0x410011 as f32),
                Complex::new(0x410012 as f32, 0x410013 as f32),
                Complex::new(0x410014 as f32, 0x410015 as f32),
                Complex::new(0x410016 as f32, 0x410017 as f32),
            ])
        );

        // ts 3, chan 3 (cc 1, fc 1), baseline 5
        let viz_3_3_1 = jones_array[(3, 3, 1)];
        compare_jones!(
            viz_3_3_1,
            Jones::from([
                Complex::new(0x410718 as f32, 0x410719 as f32),
                Complex::new(0x41071a as f32, 0x41071b as f32),
                Complex::new(0x41071c as f32, 0x41071d as f32),
                Complex::new(0x41071e as f32, 0x41071f as f32),
            ])
        );

        let integration_time_s = corr_ctx.metafits_context.corr_int_time_ms as f64 / 1000.0;

        // timestep 0
        let timestep_0 = &corr_ctx.timesteps[vis_sel.timestep_range.clone()][0];
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
        let timestep_3 = &corr_ctx.timesteps[vis_sel.timestep_range.clone()][3];
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
        let bl_1 = &corr_ctx.metafits_context.baselines[1];

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

        correct_geometry(
            &corr_ctx,
            &mut jones_array,
            &vis_sel.timestep_range,
            &vis_sel.coarse_chan_range,
            None,
            None,
            false,
        );
        // there should be no difference in baseline 0
        // ts 0 (batch 0, scan 0), chan 0 (cc 0, fc 0), baseline 0
        compare_jones!(jones_array[(0, 0, 0)], viz_0_0_0);

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

        compare_jones!(
            jones_array[(0, 0, 1)],
            Jones::from([
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

        compare_jones!(
            jones_array[(3, 3, 1)],
            Jones::from([
                Complex::new(rot_1_xx_3_3_re, rot_1_xx_3_3_im),
                Complex::new(rot_1_xy_3_3_re, rot_1_xy_3_3_im),
                Complex::new(rot_1_yx_3_3_re, rot_1_yx_3_3_im),
                Complex::new(rot_1_yy_3_3_re, rot_1_yy_3_3_im),
            ])
        );
    }

    #[test]
    fn test_correct_digital_gains() {
        let corr_ctx = get_mwa_ord_context();

        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
        let sel_baseline_range = 0..2;
        vis_sel.baseline_idxs = sel_baseline_range.clone().collect();

        // Create a blank array to store flags and visibilities
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
        let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
        // read visibilities out of the gpubox files
        vis_sel
            .read_mwalib(
                &corr_ctx,
                jones_array.view_mut(),
                flag_array.view_mut(),
                false,
            )
            .unwrap();

        let ant_pairs = vis_sel.get_ant_pairs(&corr_ctx.metafits_context);

        let first_coarse_chan_idx = vis_sel.coarse_chan_range.start;
        let last_coarse_chan_idx = vis_sel.coarse_chan_range.end - 1;
        // first coarse chan, antenna 0
        let gain_first_0 = corr_ctx.metafits_context.antennas[0]
            .rfinput_x
            .digital_gains[first_coarse_chan_idx];
        // first coarse chan, antenna 1
        let gain_first_1 = corr_ctx.metafits_context.antennas[1]
            .rfinput_x
            .digital_gains[first_coarse_chan_idx];
        // last coarse chan, antenna 0
        let gain_last_0 = corr_ctx.metafits_context.antennas[0]
            .rfinput_x
            .digital_gains[last_coarse_chan_idx];
        // last coarse chan, antenna 1
        let gain_last_1 = corr_ctx.metafits_context.antennas[1]
            .rfinput_x
            .digital_gains[last_coarse_chan_idx];

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
            &corr_ctx,
            &mut out_jones_array,
            &vis_sel.coarse_chan_range,
            &ant_pairs,
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
        let corr_ctx = get_mwa_ord_context();
        let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();

        let mut jones_array = Array3::from_shape_fn((1, 1, 1), |_| Jones::nan());

        let sel_baseline_range = 0..2;
        vis_sel.baseline_idxs = sel_baseline_range.collect();

        let ant_pairs = vis_sel.get_ant_pairs(&corr_ctx.metafits_context);

        assert!(matches!(
            correct_digital_gains(
                &corr_ctx,
                &mut jones_array,
                &vis_sel.coarse_chan_range,
                &ant_pairs,
            ),
            Err(BirliError::BadArrayShape { .. })
        ));

        let channels = vis_sel.coarse_chan_range.len()
            * corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;

        let mut jones_array = Array3::from_shape_fn((1, channels, 1), |_| Jones::nan());

        assert!(matches!(
            correct_digital_gains(
                &corr_ctx,
                &mut jones_array,
                &vis_sel.coarse_chan_range,
                &ant_pairs,
            ),
            Err(BirliError::BadArrayShape { .. })
        ));

        let gains = Array2::from_shape_fn((1, 1), |_| (1., 1.));

        assert!(matches!(
            _correct_digital_gains(&mut jones_array, &gains, &ant_pairs, 2),
            Err(BirliError::BadArrayShape { .. })
        ));
    }

    #[test]
    fn test_scrunch_gains_legacy() {
        let base: i32 = 2;
        let ultrafine_gains: Vec<f64> = (0..30).map(|x| (base.pow(x)) as _).collect();
        let expected_gains: Vec<f64> = (0..15)
            .map(|x| ((base.pow(2 * x)) + (base.pow(2 * x + 1))) as f64 / 2.)
            .collect();
        let scrunched_gains = scrunch_gains(&ultrafine_gains, 2, &ScrunchType::Simple);
        assert_eq!(scrunched_gains, expected_gains);
    }

    #[test]
    fn test_scrunch_gains_mwax_even_scrunch_even_channels() {
        let corr_ctx = get_mwax_context();
        let base: i32 = 2;
        let ultrafine_gains: Vec<f64> = (0..12).map(|x| (base.pow(x)) as _).collect();
        let expected_gains: Vec<f64> = (0..12 / 2)
            .map(|x| -> f64 {
                let left = ultrafine_gains
                    [(2 * x as i32 - 1).rem_euclid(ultrafine_gains.len() as i32) as usize]
                    as f64;
                let center = ultrafine_gains[2 * x as usize] as f64;
                let right = ultrafine_gains[(2 * x + 1) as usize] as f64;
                (left / 4. + center / 2. + right / 4.) as f64
            })
            .collect();
        let scrunched_gains = scrunch_gains(
            &ultrafine_gains,
            2,
            &ScrunchType::from_mwa_version(corr_ctx.metafits_context.mwa_version.unwrap()).unwrap(),
        );
        assert_eq!(scrunched_gains, expected_gains);
    }

    #[test]
    fn test_scrunch_gains_mwax_odd_scrunch_even_channels() {
        let base: i32 = 2;
        let ultrafine_gains: Vec<f64> = (0..12).map(|x| (base.pow(x)) as _).collect();
        let expected_gains: Vec<f64> = (0..12 / 3)
            .map(|x| -> f64 {
                let left = ultrafine_gains
                    [(3 * x as i32 - 1).rem_euclid(ultrafine_gains.len() as i32) as usize]
                    as f64;
                let center = ultrafine_gains[3 * x] as f64;
                let right = ultrafine_gains[3 * x + 1] as f64;
                (left / 3. + center / 3. + right / 3.) as f64
            })
            .collect();
        let scrunched_gains = scrunch_gains(&ultrafine_gains, 3, &ScrunchType::CenterSymmetric);
        assert_eq!(scrunched_gains, expected_gains);
    }

    #[test]
    fn test_scrunch_gains_mwax_even_scrunch_odd_channels() {
        let base: i32 = 2;
        let ultrafine_gains: Vec<f64> = (0..12).map(|x| (base.pow(x)) as _).collect();
        let expected_gains: Vec<f64> = (0..12 / 4)
            .map(|x| -> f64 {
                let left = ultrafine_gains[4 * x] as f64;
                let center1 = ultrafine_gains[4 * x + 1] as f64;
                let center2 = ultrafine_gains[4 * x + 2] as f64;
                let center3 = ultrafine_gains[4 * x + 3] as f64;
                let right = ultrafine_gains
                    [(4 * x as i32 + 4).rem_euclid(ultrafine_gains.len() as i32) as usize]
                    as f64;
                (left / 8. + center1 / 4. + center2 / 4. + center3 / 4. + right / 8.) as f64
            })
            .collect();
        let scrunched_gains = scrunch_gains(&ultrafine_gains, 4, &ScrunchType::CenterSymmetric);
        assert_eq!(scrunched_gains, expected_gains);
    }

    #[test]
    fn test_scrunch_gains_mwax_odd_scrunch_odd_channels() {
        let base: i32 = 2;
        let ultrafine_gains: Vec<f64> = (0..15).map(|x| (base.pow(x)) as _).collect();
        let expected_gains: Vec<f64> = (0..15 / 3)
            .map(|x| -> f64 {
                let left = ultrafine_gains[3 * x] as f64;
                let center1 = ultrafine_gains[3 * x + 1] as f64;
                let center2 = ultrafine_gains[3 * x + 2] as f64;
                let right = ultrafine_gains
                    [(3 * x as i32 + 3).rem_euclid(ultrafine_gains.len() as i32) as usize]
                    as f64;
                (left / 6. + center1 / 3. + center2 / 3. + right / 6.) as f64
            })
            .collect();
        let scrunched_gains = scrunch_gains(&ultrafine_gains, 3, &ScrunchType::CenterSymmetric);
        assert_eq!(scrunched_gains, expected_gains);
    }

    #[test]
    fn test_correct_coarse_passband_gains() {
        let num_fine_chans_per_coarse = 2;
        // 2 coarse channels
        let mut jones_array: Array3<Jones<f32>> =
            Array3::from_shape_fn((2, 4, 2), |_| Jones::identity());
        let mut weight_array: Array3<f32> = Array3::from_shape_fn((2, 4, 2), |_| 1.0);
        let passband_gains = vec![0.1, 0.2];

        correct_coarse_passband_gains(
            &mut jones_array,
            &mut weight_array,
            &passband_gains,
            num_fine_chans_per_coarse,
            &ScrunchType::Simple,
        )
        .unwrap();

        let exp_gains = [0.1, 0.2, 0.1, 0.2];

        for (&exp_gain, jones_array, weight_array) in izip!(
            exp_gains.iter(),
            jones_array.axis_iter(Axis(1)),
            weight_array.axis_iter(Axis(1))
        ) {
            for (&jones, &weight) in izip!(jones_array.iter(), weight_array.iter()) {
                compare_jones!(Jones::identity() / exp_gain, jones);
                assert_approx_eq!(f32, exp_gain, weight);
            }
        }
    }

    #[test]
    fn test_correct_coarse_passband_gains_good_fscrunch() {
        let num_fine_chans_per_coarse = 2;
        // 2 coarse channels
        let mut jones_array: Array3<Jones<f32>> =
            Array3::from_shape_fn((2, 4, 2), |_| Jones::identity());
        let mut weight_array: Array3<f32> = Array3::from_shape_fn((2, 4, 2), |_| 1.0);
        let passband_gains = vec![0.1, 0.2, 0.3, 0.4];

        correct_coarse_passband_gains(
            &mut jones_array,
            &mut weight_array,
            &passband_gains,
            num_fine_chans_per_coarse,
            &ScrunchType::Simple,
        )
        .unwrap();

        let exp_gains = [0.15, 0.35, 0.15, 0.35];

        for (&exp_gain, jones_array, weight_array) in izip!(
            exp_gains.iter(),
            jones_array.axis_iter(Axis(1)),
            weight_array.axis_iter(Axis(1))
        ) {
            for (&jones, &weight) in izip!(jones_array.iter(), weight_array.iter()) {
                compare_jones!(Jones::identity() / exp_gain, jones);
                assert_approx_eq!(f32, exp_gain, weight);
            }
        }
    }

    #[test]
    fn test_correct_coarse_passband_gains_high_fscrunch() {
        let num_fine_chans_per_coarse = 2;
        // 2 coarse channels
        let mut jones_array: Array3<Jones<f32>> =
            Array3::from_shape_fn((2, 4, 2), |_| Jones::identity());
        let mut weight_array: Array3<f32> = Array3::from_shape_fn((2, 4, 2), |_| 1.0);
        let passband_gains = vec![0.1, 0.2, 0.3, 100., 0.5, 0.6, 0.7, 100.];

        correct_coarse_passband_gains(
            &mut jones_array,
            &mut weight_array,
            &passband_gains,
            num_fine_chans_per_coarse,
            &ScrunchType::Simple,
        )
        .unwrap();

        let exp_gains = [25.15, 25.45, 25.15, 25.45];

        for (&exp_gain, jones_array, weight_array) in izip!(
            exp_gains.iter(),
            jones_array.axis_iter(Axis(1)),
            weight_array.axis_iter(Axis(1))
        ) {
            for (&jones, &weight) in izip!(jones_array.iter(), weight_array.iter()) {
                compare_jones!(Jones::identity() / exp_gain, jones);
                assert_approx_eq!(f32, exp_gain, weight);
            }
        }
    }

    #[test]
    fn test_correct_coarse_passband_gains_bad_array_shape() {
        let num_fine_chans_per_coarse = 2;
        let mut jones_array: Array3<Jones<f32>> =
            Array3::from_shape_fn((2, 4, 2), |_| Jones::identity());
        let mut weight_array: Array3<f32> = Array3::from_shape_fn((2, 4, 2), |_| 1.0);
        let passband_gains = vec![0.1, 0.2, 0.3, 100., 0.5, 0.6, 0.7, 100.];

        // test num_fine_chans_per_coarse=0
        assert!(matches!(
            correct_coarse_passband_gains(
                &mut jones_array,
                &mut weight_array,
                &passband_gains,
                0,
                &ScrunchType::Simple,
            ),
            Err(BirliError::BadArrayShape { .. })
        ));

        // test bad jones array shape
        assert!(matches!(
            correct_coarse_passband_gains(
                &mut jones_array,
                &mut weight_array,
                &passband_gains,
                3,
                &ScrunchType::Simple,
            ),
            Err(BirliError::BadArrayShape { .. })
        ));

        // test weight_array dimension mismatch with jones_array
        let mut bad_weight_array: Array3<f32> = Array3::from_shape_fn((2, 4, 3), |_| 1.0);
        assert!(matches!(
            correct_coarse_passband_gains(
                &mut jones_array,
                &mut bad_weight_array,
                &passband_gains,
                num_fine_chans_per_coarse,
                &ScrunchType::Simple,
            ),
            Err(BirliError::BadArrayShape { .. })
        ));

        // test bad gain shape
        let bad_passband_gains = vec![0.1, 0.2, 0.3];
        assert!(matches!(
            correct_coarse_passband_gains(
                &mut jones_array,
                &mut weight_array,
                &bad_passband_gains,
                num_fine_chans_per_coarse,
                &ScrunchType::Simple,
            ),
            Err(BirliError::BadArrayShape { .. })
        ));
    }
}
