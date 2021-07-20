//! Corrections that can be performed on visibility data

use crate::cxx_aoflagger::ffi::CxxImageSet;
use cxx::UniquePtr;
use itertools::izip;
use log::trace;
use mwalib::{CorrelatorContext, SPEED_OF_LIGHT_IN_VACUUM_M_PER_S};
use std::f64::consts::PI;

fn _correct_cable_length_buffers_cotter(
    freq_hz: &u32,
    electrical_length_m: &f64,
    buf_re: &mut [f32],
    buf_im: &mut [f32],
) {
    let angle: f64 =
        -2.0 * PI * electrical_length_m * (*freq_hz as f64) / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
    let (sin_angle_f64, cos_angle_f64) = angle.sin_cos();
    let (sin_angle, cos_angle) = (sin_angle_f64 as f32, cos_angle_f64 as f32);

    izip!(buf_re.iter_mut(), buf_im.iter_mut()).for_each(|(re, im)| {
        let vis_re = *re;
        *re = cos_angle * vis_re - sin_angle * *im;
        *im = sin_angle * vis_re + cos_angle * *im;
    })
}

/// A vector of all fine channel frequencies for the provided coarse channels
/// in order of provided `coarse_chan_idxs`, and ascending fine channel index
fn _get_all_freqs_hz(context: &CorrelatorContext, coarse_chan_idxs: &[usize]) -> Vec<u32> {
    let coarse_chans = context.coarse_chans.clone();
    let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
    let fine_chan_width_hz = context.metafits_context.corr_fine_chan_width_hz;
    coarse_chan_idxs
        .iter()
        .flat_map(|&coarse_chan_idx| {
            let coarse_chan = &coarse_chans[coarse_chan_idx];
            let chan_start_hz = coarse_chan.chan_start_hz;
            (0..fine_chans_per_coarse).map(move |fine_chan_idx| {
                chan_start_hz + (fine_chan_idx as u32 * fine_chan_width_hz)
            })
        })
        .collect()
}

/// Perform cable length corrections, given an observation's
/// [`mwalib::CorrelatorContext`] and a vector of [`CxxImageSet`]s for  each
/// baseline.
///
/// Cable lengths are determined by the difference between a baseline's rfInput
/// electrical lengths in the metafits. Complex visibilities are phase-shifted
/// by an angle determined by the electrical length, and the channel's
/// frequency.
///
/// # Examples
///
/// ```rust
/// use birli::{context_to_baseline_imgsets, cxx_aoflagger_new, correct_cable_lengths};
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
/// // Create an mwalib::CorrelatorContext for accessing visibilities.
/// let context = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
///
/// // Determine which timesteps and coarse channels we want to use
/// let img_coarse_chan_idxs = &context.common_coarse_chan_indices;
/// let img_timestep_idxs = &context.common_timestep_indices;
///
/// // create a CxxAOFlagger object to perform AOFlagger operations
/// let aoflagger = unsafe { cxx_aoflagger_new() };
///
/// // generate imagesets for each baseline in the format required by aoflagger
/// let mut baseline_imgsets = context_to_baseline_imgsets(
///     &aoflagger,
///     &context,
///     &img_coarse_chan_idxs,
///     &img_timestep_idxs,
///     None,
/// );
///
/// correct_cable_lengths(&context, &mut baseline_imgsets, &img_coarse_chan_idxs);
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
/// Here is the more accurate implementation of
/// [`_correct_cable_length_buffers_cotter`]
///
/// ```rust
/// use itertools::izip;
/// use std::f64::consts::PI;
///
/// const SPEED_OF_LIGHT_IN_VACUUM_M_PER_S: f64 = 299792458.0; // speed of light in m/s
///
/// fn _correct_cable_length_buffers_precise(
///     freq_hz: &u32,
///     electrical_length_m: &f64,
///     buf_re: &mut [f32],
///     buf_im: &mut [f32],
/// ) {
///     let angle: f64 = -2.0 * PI * electrical_length_m * (*freq_hz as f64) / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
///     let (sin_angle, cos_angle) = angle.sin_cos();
///     izip!(buf_re.iter_mut(), buf_im.iter_mut()).for_each(|(re, im)| {
///         let vis_re = *re as f64;
///         let vis_im = *im as f64;
///         *re = (cos_angle * vis_re - sin_angle * vis_im) as f32;
///         *im = (sin_angle * vis_re + cos_angle * vis_im) as f32;
///     })
/// }
/// ```
pub fn correct_cable_lengths(
    context: &CorrelatorContext,
    baseline_imgsets: &mut Vec<UniquePtr<CxxImageSet>>,
    img_coarse_chan_idxs: &[usize],
) {
    trace!("start correct_cable_lengths");

    let baselines = &context.metafits_context.baselines;
    let antennas = &context.metafits_context.antennas;
    // let coarse_chans = &context.coarse_chans;

    // let fine_chans_per_coarse = context.metafits_context.num_corr_fine_chans_per_coarse;
    // let fine_chan_width_hz = &context.metafits_context.corr_fine_chan_width_hz;

    // A vector of all fine channel frequencies for all coarse channels in ascending order.
    // let all_freqs_hz: Vec<u32> = coarse_chans
    //     .iter()
    //     .flat_map(|coarse_chan| {
    //         let chan_start_hz = coarse_chan.chan_start_hz;
    //         (0..fine_chans_per_coarse).map(move |fine_chan_idx| {
    //             chan_start_hz + (fine_chan_idx as u32 * fine_chan_width_hz)
    //         })
    //     })
    //     .collect();

    let all_freqs_hz = _get_all_freqs_hz(&context, img_coarse_chan_idxs);

    izip!(baseline_imgsets.iter_mut(), baselines).for_each(|(imgset, baseline)| {
        let imgset_stride: usize = imgset.HorizontalStride();
        let ant1 = &antennas[baseline.ant1_index];
        let ant2 = &antennas[baseline.ant2_index];

        // TODO: skip if baseline.ant1_index == baseline.ant2_index ?

        let pol_lengths = vec![
            ant2.rfinput_x.electrical_length_m - ant1.rfinput_x.electrical_length_m,
            ant2.rfinput_y.electrical_length_m - ant1.rfinput_x.electrical_length_m,
            ant2.rfinput_x.electrical_length_m - ant1.rfinput_y.electrical_length_m,
            ant2.rfinput_y.electrical_length_m - ant1.rfinput_y.electrical_length_m,
        ];
        pol_lengths
            .iter()
            .enumerate()
            .for_each(|(pol_idx, electrical_length_m)| {
                let imgset_buf_re: &mut [f32] = unsafe { imgset.ImageBufferMutUnsafe(2 * pol_idx) };
                let imgset_buf_im: &mut [f32] =
                    unsafe { imgset.ImageBufferMutUnsafe(2 * pol_idx + 1) };

                izip!(
                    all_freqs_hz.iter(),
                    imgset_buf_re.chunks_mut(imgset_stride),
                    imgset_buf_im.chunks_mut(imgset_stride)
                )
                .for_each(|(freq_hz, imgset_chunk_re, imgset_chunk_im)| {
                    // _correct_cable_length_buffers_precise(
                    _correct_cable_length_buffers_cotter(
                        freq_hz,
                        electrical_length_m,
                        imgset_chunk_re,
                        imgset_chunk_im,
                    )
                });
            });
    });
    trace!("end correct_cable_lengths");
}

#[cfg(test)]
mod tests {
    #![allow(clippy::float_cmp)]

    use super::{correct_cable_lengths, SPEED_OF_LIGHT_IN_VACUUM_M_PER_S};
    use float_cmp::{approx_eq, F32Margin};
    use mwalib::CorrelatorContext;

    use crate::{
        context_to_baseline_imgsets, corrections::_get_all_freqs_hz, cxx_aoflagger_new,
        get_flaggable_timesteps,
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

    macro_rules! test_imgset_val {
        ($imgset:expr, $imgset_idx:expr, $fr:expr, $ts:expr, $val:expr) => {
            let left = $imgset.ImageBuffer($imgset_idx)[$fr * $imgset.HorizontalStride() + $ts];
            let right = $val as f32;
            assert!(
                approx_eq!(f32, left, right, F32Margin::default()),
                "{} (0x{:08x}) != {} (0x{:08x})",
                left,
                left as i32,
                right,
                right as i32
            )
        };
    }

    #[test]
    fn test_cable_length_corrections_mwax() {
        let context = get_mwax_context();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(img_timestep_idxs.len(), 4);
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;

        let all_freqs_hz = _get_all_freqs_hz(&context, img_coarse_chan_idxs);

        let aoflagger = unsafe { cxx_aoflagger_new() };

        let width = img_timestep_idxs.len();
        let stride = (((width - 1) / 8) + 1) * 8;

        let mut baseline_imgsets = context_to_baseline_imgsets(
            &aoflagger,
            &context,
            &img_coarse_chan_idxs,
            &img_timestep_idxs,
            None,
        );

        let viz_0_xx_0_0_re = baseline_imgsets[0].ImageBuffer(0)[0];
        let viz_0_xx_0_0_im = baseline_imgsets[0].ImageBuffer(1)[0];
        let viz_1_xx_0_0_re = baseline_imgsets[1].ImageBuffer(0)[0];
        let viz_1_xx_0_0_im = baseline_imgsets[1].ImageBuffer(1)[0];
        let viz_1_xy_0_0_re = baseline_imgsets[1].ImageBuffer(2)[0];
        let viz_1_xy_0_0_im = baseline_imgsets[1].ImageBuffer(3)[0];
        let viz_1_yx_0_0_re = baseline_imgsets[1].ImageBuffer(4)[0];
        let viz_1_yx_0_0_im = baseline_imgsets[1].ImageBuffer(5)[0];
        let viz_1_yy_0_0_re = baseline_imgsets[1].ImageBuffer(6)[0];
        let viz_1_yy_0_0_im = baseline_imgsets[1].ImageBuffer(7)[0];
        let viz_1_yy_3_3_re = baseline_imgsets[1].ImageBuffer(6)[3 * stride + 3];
        let viz_1_yy_3_3_im = baseline_imgsets[1].ImageBuffer(7)[3 * stride + 3];

        assert_eq!(viz_0_xx_0_0_re as f32, 0x410000 as f32);
        assert_eq!(viz_0_xx_0_0_im as f32, 0x410001 as f32);
        assert_eq!(viz_1_xx_0_0_re as f32, 0x410010 as f32);
        assert_eq!(viz_1_xx_0_0_im as f32, 0x410011 as f32);
        assert_eq!(viz_1_xy_0_0_re as f32, 0x410012 as f32);
        assert_eq!(viz_1_xy_0_0_im as f32, 0x410013 as f32);
        assert_eq!(viz_1_yx_0_0_re as f32, 0x410014 as f32);
        assert_eq!(viz_1_yx_0_0_im as f32, 0x410015 as f32);
        assert_eq!(viz_1_yy_0_0_re as f32, 0x410016 as f32);
        assert_eq!(viz_1_yy_0_0_im as f32, 0x410017 as f32);
        assert_eq!(viz_1_yy_3_3_re as f32, 0x41071e as f32);
        assert_eq!(viz_1_yy_3_3_im as f32, 0x41071f as f32);

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

        // baseline 1, pol XX, cc 0, fc 0
        let length_m_1_xx = length_1_2_x - length_1_1_x;
        let angle_1_xx_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_1_xx * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
        let (sin_1_xx_0_f64, cos_1_xx_0_f64) = angle_1_xx_0.sin_cos();
        let (sin_1_xx_0, cos_1_xx_0) = (sin_1_xx_0_f64 as f32, cos_1_xx_0_f64 as f32);
        // baseline 1, pol XY, cc 0, fc 0
        let length_m_1_xy = length_1_2_y - length_1_1_x;
        let angle_1_xy_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_1_xy * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
        let (sin_1_xy_0_f64, cos_1_xy_0_f64) = angle_1_xy_0.sin_cos();
        let (sin_1_xy_0, cos_1_xy_0) = (sin_1_xy_0_f64 as f32, cos_1_xy_0_f64 as f32);
        // baseline 1, pol YX, cc 0, fc 0
        let length_m_1_yx = length_1_2_x - length_1_1_y;
        let angle_1_yx_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_1_yx * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
        let (sin_1_yx_0_f64, cos_1_yx_0_f64) = angle_1_yx_0.sin_cos();
        let (sin_1_yx_0, cos_1_yx_0) = (sin_1_yx_0_f64 as f32, cos_1_yx_0_f64 as f32);
        // baseline 1, pol YY, cc 0, fc 0
        let length_m_1_yy = length_1_2_y - length_1_1_y;
        let angle_1_yy_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_1_yy * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
        let (sin_1_yy_0_f64, cos_1_yy_0_f64) = angle_1_yy_0.sin_cos();
        let (sin_1_yy_0, cos_1_yy_0) = (sin_1_yy_0_f64 as f32, cos_1_yy_0_f64 as f32);
        // baseline 1, pol YY, cc 1, fc 1
        let angle_1_yy_3: f64 =
            -2.0 * std::f64::consts::PI * length_m_1_yy * (all_freqs_hz[3] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
        let (sin_1_yy_3_f64, cos_1_yy_3_f64) = angle_1_yy_3.sin_cos();
        let (sin_1_yy_3, cos_1_yy_3) = (sin_1_yy_3_f64 as f32, cos_1_yy_3_f64 as f32);

        correct_cable_lengths(&context, &mut baseline_imgsets, img_coarse_chan_idxs);

        // there should be no difference in baseline 0
        test_imgset_val!(baseline_imgsets[0], 0, 0, 0, viz_0_xx_0_0_re);
        test_imgset_val!(baseline_imgsets[0], 1, 0, 0, viz_0_xx_0_0_im);

        ////
        // baseline 1 should be rotated
        ////
        // baseline 1, pol xx, cc 0, fc 0, ts 0
        let rot_1_xx_0_0_re = (cos_1_xx_0 * viz_1_xx_0_0_re - sin_1_xx_0 * viz_1_xx_0_0_im) as f32;
        let rot_1_xx_0_0_im = (sin_1_xx_0 * viz_1_xx_0_0_re + cos_1_xx_0 * viz_1_xx_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[1], 0, 0, 0, rot_1_xx_0_0_re);
        test_imgset_val!(baseline_imgsets[1], 1, 0, 0, rot_1_xx_0_0_im);
        // baseline 1, pol xy, cc 0, fc 0, ts 0
        let rot_1_xy_0_0_re = (cos_1_xy_0 * viz_1_xy_0_0_re - sin_1_xy_0 * viz_1_xy_0_0_im) as f32;
        let rot_1_xy_0_0_im = (sin_1_xy_0 * viz_1_xy_0_0_re + cos_1_xy_0 * viz_1_xy_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[1], 2, 0, 0, rot_1_xy_0_0_re);
        test_imgset_val!(baseline_imgsets[1], 3, 0, 0, rot_1_xy_0_0_im);
        // baseline 1, pol yx, cc 0, fc 0, ts 0
        let rot_1_yx_0_0_re = (cos_1_yx_0 * viz_1_yx_0_0_re - sin_1_yx_0 * viz_1_yx_0_0_im) as f32;
        let rot_1_yx_0_0_im = (sin_1_yx_0 * viz_1_yx_0_0_re + cos_1_yx_0 * viz_1_yx_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[1], 4, 0, 0, rot_1_yx_0_0_re);
        test_imgset_val!(baseline_imgsets[1], 5, 0, 0, rot_1_yx_0_0_im);
        // baseline 1, pol yy, cc 0, fc 0, ts 0
        let rot_1_yy_0_0_re = (cos_1_yy_0 * viz_1_yy_0_0_re - sin_1_yy_0 * viz_1_yy_0_0_im) as f32;
        let rot_1_yy_0_0_im = (sin_1_yy_0 * viz_1_yy_0_0_re + cos_1_yy_0 * viz_1_yy_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[1], 6, 0, 0, rot_1_yy_0_0_re);
        test_imgset_val!(baseline_imgsets[1], 7, 0, 0, rot_1_yy_0_0_im);
        // baseline 1, pol yy, cc 1, fc 1, ts 1
        let rot_1_yy_3_3_re = (cos_1_yy_3 * viz_1_yy_3_3_re - sin_1_yy_3 * viz_1_yy_3_3_im) as f32;
        let rot_1_yy_3_3_im = (sin_1_yy_3 * viz_1_yy_3_3_re + cos_1_yy_3 * viz_1_yy_3_3_im) as f32;
        test_imgset_val!(baseline_imgsets[1], 6, 3, 3, rot_1_yy_3_3_re);
        test_imgset_val!(baseline_imgsets[1], 7, 3, 3, rot_1_yy_3_3_im);
    }

    #[test]
    fn test_cable_length_corrections_ord() {
        let context = get_mwa_ord_context();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        assert_eq!(img_timestep_idxs.len(), 4);
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;

        let all_freqs_hz = _get_all_freqs_hz(&context, img_coarse_chan_idxs);

        let aoflagger = unsafe { cxx_aoflagger_new() };

        let width = img_timestep_idxs.len();
        let stride = (((width - 1) / 8) + 1) * 8;

        let mut baseline_imgsets = context_to_baseline_imgsets(
            &aoflagger,
            &context,
            &img_coarse_chan_idxs,
            &img_timestep_idxs,
            None,
        );

        let viz_0_xx_0_0_re = baseline_imgsets[0].ImageBuffer(0)[0];
        let viz_0_xx_0_0_im = baseline_imgsets[0].ImageBuffer(1)[0];
        let viz_5_xx_0_0_re = baseline_imgsets[5].ImageBuffer(0)[0];
        let viz_5_xx_0_0_im = baseline_imgsets[5].ImageBuffer(1)[0];
        let viz_5_xy_0_0_re = baseline_imgsets[5].ImageBuffer(2)[0];
        let viz_5_xy_0_0_im = baseline_imgsets[5].ImageBuffer(3)[0];
        let viz_5_yx_0_0_re = baseline_imgsets[5].ImageBuffer(4)[0];
        let viz_5_yx_0_0_im = baseline_imgsets[5].ImageBuffer(5)[0];
        let viz_5_yy_0_0_re = baseline_imgsets[5].ImageBuffer(6)[0];
        let viz_5_yy_0_0_im = baseline_imgsets[5].ImageBuffer(7)[0];
        let viz_5_yy_3_3_re = baseline_imgsets[5].ImageBuffer(6)[3 * stride + 3];
        let viz_5_yy_3_3_im = baseline_imgsets[5].ImageBuffer(7)[3 * stride + 3];

        assert_eq!(viz_0_xx_0_0_re as f32, 0x10c5be as f32);
        assert_eq!(viz_0_xx_0_0_im as f32, 0x10c5bf as f32);
        assert_eq!(viz_5_xx_0_0_re as f32, 0x10f1ce as f32);
        assert_eq!(viz_5_xx_0_0_im as f32, -0x10f1cf as f32);
        assert_eq!(viz_5_xy_0_0_re as f32, 0x10ea26 as f32);
        assert_eq!(viz_5_xy_0_0_im as f32, -0x10ea27 as f32);
        assert_eq!(viz_5_yx_0_0_re as f32, 0x10f1be as f32);
        assert_eq!(viz_5_yx_0_0_im as f32, -0x10f1bf as f32);
        assert_eq!(viz_5_yy_0_0_re as f32, 0x10ea16 as f32);
        assert_eq!(viz_5_yy_0_0_im as f32, -0x10ea17 as f32);
        assert_eq!(viz_5_yy_3_3_re as f32, 0x0dec16 as f32);
        assert_eq!(viz_5_yy_3_3_im as f32, -0x0dec17 as f32);

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

        // baseline 1, pol XX, cc 0, fc 0
        let length_m_5_xx = length_5_2_x - length_5_1_x;
        let angle_5_xx_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_5_xx * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
        let (sin_5_xx_0_f64, cos_5_xx_0_f64) = angle_5_xx_0.sin_cos();
        let (sin_5_xx_0, cos_5_xx_0) = (sin_5_xx_0_f64 as f32, cos_5_xx_0_f64 as f32);
        // baseline 1, pol XY, cc 0, fc 0
        let length_m_5_xy = length_5_2_y - length_5_1_x;
        let angle_5_xy_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_5_xy * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
        let (sin_5_xy_0_f64, cos_5_xy_0_f64) = angle_5_xy_0.sin_cos();
        let (sin_5_xy_0, cos_5_xy_0) = (sin_5_xy_0_f64 as f32, cos_5_xy_0_f64 as f32);
        // baseline 1, pol YX, cc 0, fc 0
        let length_m_5_yx = length_5_2_x - length_5_1_y;
        let angle_5_yx_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_5_yx * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
        let (sin_5_yx_0_f64, cos_5_yx_0_f64) = angle_5_yx_0.sin_cos();
        let (sin_5_yx_0, cos_5_yx_0) = (sin_5_yx_0_f64 as f32, cos_5_yx_0_f64 as f32);
        // baseline 1, pol YY, cc 0, fc 0
        let length_m_5_yy = length_5_2_y - length_5_1_y;
        let angle_5_yy_0: f64 =
            -2.0 * std::f64::consts::PI * length_m_5_yy * (all_freqs_hz[0] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
        let (sin_5_yy_0_f64, cos_5_yy_0_f64) = angle_5_yy_0.sin_cos();
        let (sin_5_yy_0, cos_5_yy_0) = (sin_5_yy_0_f64 as f32, cos_5_yy_0_f64 as f32);
        // baseline 1, pol YY, cc 1, fc 1
        let angle_5_yy_3: f64 =
            -2.0 * std::f64::consts::PI * length_m_5_yy * (all_freqs_hz[3] as f64)
                / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
        let (sin_5_yy_3_f64, cos_5_yy_3_f64) = angle_5_yy_3.sin_cos();
        let (sin_5_yy_3, cos_5_yy_3) = (sin_5_yy_3_f64 as f32, cos_5_yy_3_f64 as f32);

        correct_cable_lengths(&context, &mut baseline_imgsets, img_coarse_chan_idxs);

        // there should be no difference in baseline 0
        test_imgset_val!(baseline_imgsets[0], 0, 0, 0, viz_0_xx_0_0_re);
        test_imgset_val!(baseline_imgsets[0], 1, 0, 0, viz_0_xx_0_0_im);

        ////
        // baseline 1 should be rotated
        ////
        // baseline 1, pol xx, cc 0, fc 0, ts 0
        let rot_5_xx_0_0_re = (cos_5_xx_0 * viz_5_xx_0_0_re - sin_5_xx_0 * viz_5_xx_0_0_im) as f32;
        let rot_5_xx_0_0_im = (sin_5_xx_0 * viz_5_xx_0_0_re + cos_5_xx_0 * viz_5_xx_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[5], 0, 0, 0, rot_5_xx_0_0_re);
        test_imgset_val!(baseline_imgsets[5], 1, 0, 0, rot_5_xx_0_0_im);
        // baseline 1, pol xy, cc 0, fc 0, ts 0
        let rot_5_xy_0_0_re = (cos_5_xy_0 * viz_5_xy_0_0_re - sin_5_xy_0 * viz_5_xy_0_0_im) as f32;
        let rot_5_xy_0_0_im = (sin_5_xy_0 * viz_5_xy_0_0_re + cos_5_xy_0 * viz_5_xy_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[5], 2, 0, 0, rot_5_xy_0_0_re);
        test_imgset_val!(baseline_imgsets[5], 3, 0, 0, rot_5_xy_0_0_im);
        // baseline 1, pol yx, cc 0, fc 0, ts 0
        let rot_5_yx_0_0_re = (cos_5_yx_0 * viz_5_yx_0_0_re - sin_5_yx_0 * viz_5_yx_0_0_im) as f32;
        let rot_5_yx_0_0_im = (sin_5_yx_0 * viz_5_yx_0_0_re + cos_5_yx_0 * viz_5_yx_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[5], 4, 0, 0, rot_5_yx_0_0_re);
        test_imgset_val!(baseline_imgsets[5], 5, 0, 0, rot_5_yx_0_0_im);
        // baseline 1, pol yy, cc 0, fc 0, ts 0
        let rot_5_yy_0_0_re = (cos_5_yy_0 * viz_5_yy_0_0_re - sin_5_yy_0 * viz_5_yy_0_0_im) as f32;
        let rot_5_yy_0_0_im = (sin_5_yy_0 * viz_5_yy_0_0_re + cos_5_yy_0 * viz_5_yy_0_0_im) as f32;
        test_imgset_val!(baseline_imgsets[5], 6, 0, 0, rot_5_yy_0_0_re);
        test_imgset_val!(baseline_imgsets[5], 7, 0, 0, rot_5_yy_0_0_im);
        // baseline 1, pol yy, cc 1, fc 1, ts 1
        let rot_5_yy_3_3_re = (cos_5_yy_3 * viz_5_yy_3_3_re - sin_5_yy_3 * viz_5_yy_3_3_im) as f32;
        let rot_5_yy_3_3_im = (sin_5_yy_3 * viz_5_yy_3_3_re + cos_5_yy_3 * viz_5_yy_3_3_im) as f32;
        test_imgset_val!(baseline_imgsets[5], 6, 3, 3, rot_5_yy_3_3_re);
        test_imgset_val!(baseline_imgsets[5], 7, 3, 3, rot_5_yy_3_3_im);
    }
}
