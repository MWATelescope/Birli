//! Calibrating visibilities.

use crate::ndarray::{ArrayView1, ArrayView2, ArrayViewMut3, Axis};
use itertools::izip;
use marlu::Jones;
use thiserror::Error;

#[derive(Error, Debug)]
/// Errors that can occur when calibrating visibilities.
pub enum CalibrationError {
    #[error("The provided calibration solution has {calsol_chans} channels, but the data has {data_chans} channels")]
    /// When the number of channels in the calibration solution does not match the number of channels in the data.
    ChannelSizeMismatch {
        /// The argument name within the funciton
        calsol_chans: usize,
        /// The number of channels in the data
        data_chans: usize,
    },

    #[error("bad array shape supplied to argument {argument} of function {function}. expected {expected}, received {received}")]
    /// Error for bad array shape in provided argument
    BadArrayShape {
        /// The argument name within the funciton
        argument: String,
        /// The function name
        function: String,
        /// The expected shape
        expected: String,
        /// The shape that was received instead
        received: String,
    },
}

/// Returns whether every channel in a tile's calibration solution contains at least one NaN.
fn tile_calsol_is_nan_flagged(tile_row: ArrayView1<Jones<f64>>) -> bool {
    !tile_row.is_empty() && tile_row.iter().all(|j| j.any_nan())
}

/// Get a per-tile mask for tiles whose calibration solutions contain NaN in every channel.
pub fn get_calsol_nan_flagged_tiles(calsols: ArrayView2<Jones<f64>>) -> Vec<bool> {
    calsols
        .axis_iter(Axis(0))
        .map(tile_calsol_is_nan_flagged)
        .collect()
}

/// Flag antennas whose calibration solutions contain NaN in every channel.
pub fn flag_antennas_with_nan_calsols(antenna_flags: &mut [bool], calsols: ArrayView2<Jones<f64>>) {
    for (tile_idx, tile_row) in calsols.axis_iter(Axis(0)).enumerate() {
        if tile_idx >= antenna_flags.len() {
            break;
        }
        if tile_calsol_is_nan_flagged(tile_row) {
            antenna_flags[tile_idx] = true;
        }
    }
}

/// apply a direction independent calibration solution for a single timeblock to the given
/// visibility data
///
/// Baselines involving a tile marked in `flagged_tiles`, or a calibration solution containing
/// NaN, are flagged and calibration is not applied.
///
/// # Errors
///
/// calsols should have the same number of channels as `vis_array`, `flag_array`, `weight_array` etc.
///
#[allow(clippy::too_many_arguments)]
pub fn apply_di_calsol(
    // a two dimensional array of jones matrix calibration solutions with
    // dimensions `[tile][channel]`
    calsols: ArrayView2<Jones<f64>>,
    // dimensions `[timestep][channel][baselines]`
    mut vis_array: ArrayViewMut3<Jones<f32>>,
    // dimensions `[timestep][channel][baselines]`
    mut weight_array: ArrayViewMut3<f32>,
    // dimensions `[timestep][channel][baselines]`
    // todo: setting both flags and weights is redundant, but it's not clear how to rip this out
    mut flag_array: ArrayViewMut3<bool>,
    // The tile index pairs for each selected baseline
    sel_baselines: &[(usize, usize)],
    // Per-tile flag mask indexed by antenna/tile number (e.g. from [`FlagContext::antenna_flags`])
    flagged_tiles: &[bool],
    // If we are emulating cotter we don't flag NaNs
    emulate_cotter: bool,
) -> Result<(), CalibrationError> {
    let di_dims = calsols.dim();
    let vis_dims = vis_array.dim();
    let weight_dims = weight_array.dim();
    let flag_dims = flag_array.dim();
    if weight_dims != vis_dims {
        return Err(CalibrationError::BadArrayShape {
            argument: "weight_array".into(),
            function: "apply_di_calsol".into(),
            expected: format!("{vis_dims:?}"),
            received: format!("{weight_dims:?}"),
        });
    }
    if flag_dims != vis_dims {
        return Err(CalibrationError::BadArrayShape {
            argument: "flag_array".into(),
            function: "apply_di_calsol".into(),
            expected: format!("{vis_dims:?}"),
            received: format!("{flag_dims:?}"),
        });
    }

    if (vis_dims.1 as f64 / di_dims.1 as f64).fract().abs() > 0.01 {
        return Err(CalibrationError::ChannelSizeMismatch {
            calsol_chans: di_dims.1,
            data_chans: vis_dims.1,
        });
    }
    let channel_ratio = (vis_dims.1 as f64 / di_dims.1 as f64).round() as usize;

    // time axis
    for (mut vis_array, mut weight_array, mut flag_array) in izip!(
        vis_array.axis_iter_mut(Axis(0)),
        weight_array.axis_iter_mut(Axis(0)),
        flag_array.axis_iter_mut(Axis(0)),
    ) {
        // baseline axis
        for (&(ant1_idx, ant2_idx), mut vis_array, mut weight_array, mut flag_array) in izip!(
            sel_baselines.iter(),
            vis_array.axis_iter_mut(Axis(1)),
            weight_array.axis_iter_mut(Axis(1)),
            flag_array.axis_iter_mut(Axis(1)),
        ) {
            let baseline_tile_flagged = flagged_tiles.get(ant1_idx).copied().unwrap_or(false)
                || flagged_tiles.get(ant2_idx).copied().unwrap_or(false);

            // channel axis (chunked by channel_ratio)
            for (&sol1, &sol2, mut vis_chunk, mut weight_chunk, mut flag_chunk) in izip!(
                calsols.index_axis(Axis(0), ant1_idx),
                calsols.index_axis(Axis(0), ant2_idx),
                vis_array.axis_chunks_iter_mut(Axis(0), channel_ratio),
                weight_array.axis_chunks_iter_mut(Axis(0), channel_ratio),
                flag_array.axis_chunks_iter_mut(Axis(0), channel_ratio),
            ) {
                // apply the calibration solution to all visibilities in the chunk
                for (vis, weight, flag) in izip!(
                    vis_chunk.iter_mut(),
                    weight_chunk.iter_mut(),
                    flag_chunk.iter_mut()
                ) {
                    // ignore nans if we are emulating cotter
                    if baseline_tile_flagged
                        || (!emulate_cotter && (sol1.any_nan() || sol2.any_nan()))
                    {
                        *flag = true;
                        if *weight > 0. {
                            *weight = -*weight;
                        }
                        continue;
                    }

                    // promote
                    let vis_f64 = Jones::<f64>::from(*vis);

                    // demote J1 * D * J2^H
                    *vis = Jones::<f32>::from(sol1 * vis_f64 * sol2.h());

                    // if the data now contains a NaN, flag it
                    // Cotter doesn't do this, so skip when emulating it.
                    if !emulate_cotter && vis.any_nan() {
                        *flag = true;
                        if *weight > 0. {
                            *weight = -*weight;
                        }
                    }
                }
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use crate::{compare_jones, Complex};

    use ndarray::{array, Array2, Array3};

    use super::*;

    /// Test the calsols are correctly applied in the antenna axis.
    #[test]
    fn test_apply_calsols_antenna() {
        let sel_baselines = vec![(0, 0), (0, 1), (1, 1)];
        let num_times = 1;

        let calsols = Array2::from_shape_fn((2, 1), |(i, _)| Jones::identity() * (i + 1) as f64);
        let shape = (num_times, calsols.dim().1, sel_baselines.len());
        let mut vis_array = Array3::from_shape_fn(shape, |(_, _, bl)| {
            Jones::<f32>::identity() * (bl + 1) as f32
        });
        let mut flag_array = Array3::from_shape_fn(shape, |_| false);
        let mut weight_array = Array3::from_shape_fn(shape, |_| 1_f32);
        let emulate_cotter = false;
        apply_di_calsol(
            calsols.view(),
            vis_array.view_mut(),
            weight_array.view_mut(),
            flag_array.view_mut(),
            &sel_baselines,
            &[],
            emulate_cotter,
        )
        .unwrap();

        compare_jones!(
            vis_array[(0, 0, 0)],
            calsols[(0, 0)] * (Jones::<f64>::identity() * 1.) * calsols[(0, 0)].h()
        );
        compare_jones!(
            vis_array[(0, 0, 1)],
            calsols[(0, 0)] * (Jones::<f64>::identity() * 2.) * calsols[(1, 0)].h()
        );
        compare_jones!(
            vis_array[(0, 0, 2)],
            calsols[(1, 0)] * (Jones::<f64>::identity() * 3.) * calsols[(1, 0)].h()
        );
    }

    /// Test the calsols are correctly applied in the channel axis.
    #[test]
    fn test_apply_calsols_chan() {
        let sel_baselines = vec![(0, 0)];
        let num_times = 1;

        let calsols =
            Array2::from_shape_fn((1, 2), |(_, c)| Jones::identity() * (c * 2 + 1) as f64);
        let shape = (num_times, calsols.dim().1, sel_baselines.len());
        let mut vis_array = Array3::from_shape_fn(shape, |(_, c, _)| {
            Jones::<f32>::identity() * (c * 2 + 2) as f32
        });
        let mut flag_array = Array3::from_shape_fn(shape, |_| false);
        let mut weight_array = Array3::from_shape_fn(shape, |_| 1_f32);
        let emulate_cotter = false;
        apply_di_calsol(
            calsols.view(),
            vis_array.view_mut(),
            weight_array.view_mut(),
            flag_array.view_mut(),
            &sel_baselines,
            &[],
            emulate_cotter,
        )
        .unwrap();

        compare_jones!(
            vis_array[(0, 0, 0)],
            calsols[(0, 0)] * (Jones::<f64>::identity() * 2.) * calsols[(0, 0)].h()
        );
        compare_jones!(
            vis_array[(0, 1, 0)],
            calsols[(0, 1)] * (Jones::<f64>::identity() * 4.) * calsols[(0, 1)].h()
        );
    }

    /// Test the calsols are correctly applied in the channel axis, when vis chans = 2 * cal chans.
    #[test]
    fn test_apply_calsols_chan_uneven() {
        let sel_baselines = vec![(0, 0)];
        let num_times = 1;

        let calsols =
            Array2::from_shape_fn((1, 2), |(_, c)| Jones::identity() * (c * 2 + 1) as f64);
        let shape = (num_times, calsols.dim().1 * 2, sel_baselines.len());
        let mut vis_array = Array3::from_shape_fn(shape, |(_, c, _)| {
            Jones::<f32>::identity() * (c * 2 + 2) as f32
        });
        let mut flag_array = Array3::from_shape_fn(shape, |_| false);
        let mut weight_array = Array3::from_shape_fn(shape, |_| 1_f32);
        let emulate_cotter = false;
        apply_di_calsol(
            calsols.view(),
            vis_array.view_mut(),
            weight_array.view_mut(),
            flag_array.view_mut(),
            &sel_baselines,
            &[],
            emulate_cotter,
        )
        .unwrap();

        compare_jones!(
            vis_array[(0, 0, 0)],
            calsols[(0, 0)] * (Jones::<f64>::identity() * 2.) * calsols[(0, 0)].h()
        );
        compare_jones!(
            vis_array[(0, 1, 0)],
            calsols[(0, 0)] * (Jones::<f64>::identity() * 4.) * calsols[(0, 0)].h()
        );
        compare_jones!(
            vis_array[(0, 2, 0)],
            calsols[(0, 1)] * (Jones::<f64>::identity() * 6.) * calsols[(0, 1)].h()
        );
        compare_jones!(
            vis_array[(0, 3, 0)],
            calsols[(0, 1)] * (Jones::<f64>::identity() * 8.) * calsols[(0, 1)].h()
        );
    }

    /// Test the calsols are correctly applied in to all timesteps.
    #[test]
    fn test_apply_calsols_time() {
        let sel_baselines = vec![(0, 0)];
        let num_times = 2;

        let calsols = Array2::from_shape_fn((1, 1), |_| Jones::identity() * 2.);
        let shape = (num_times, calsols.dim().1, sel_baselines.len());
        let mut vis_array = Array3::from_shape_fn(shape, |(t, _, _)| {
            Jones::<f32>::identity() * (t * 2 + 2) as f32
        });
        let mut flag_array = Array3::from_shape_fn(shape, |_| false);
        let mut weight_array = Array3::from_shape_fn(shape, |_| 1_f32);
        let emulate_cotter = false;
        apply_di_calsol(
            calsols.view(),
            vis_array.view_mut(),
            weight_array.view_mut(),
            flag_array.view_mut(),
            &sel_baselines,
            &[],
            emulate_cotter,
        )
        .unwrap();

        compare_jones!(
            vis_array[(0, 0, 0)],
            calsols[(0, 0)] * (Jones::<f64>::identity() * 2.) * calsols[(0, 0)].h()
        );
        compare_jones!(
            vis_array[(1, 0, 0)],
            calsols[(0, 0)] * (Jones::<f64>::identity() * 4.) * calsols[(0, 0)].h()
        );
    }

    /// Test the calsols are correctly applied based on real values from cotter debugger.
    #[test]
    fn test_apply_calsols_real() {
        let sel_baselines = vec![(0, 1)];

        let calsols: Array2<Jones<f64>> = array![
            // -exec p solA[solChannel]
            [
                Jones::from([
                    Complex::new(-0.05711880819681107, 0.8909723224701427),
                    Complex::new(0., 0.),
                    Complex::new(0., 0.),
                    Complex::new(-0.3190681285208096, 0.8975262420831493)
                ]),
                Jones::from([
                    Complex::new(-0.05790403500446751, 0.8906022388084277),
                    Complex::new(0., 0.),
                    Complex::new(0., 0.),
                    Complex::new(-0.31938558050469074, 0.8973555420886708)
                ]),
            ],
            // -exec p solB[solChannel]
            [
                Jones::from([
                    Complex::new(0.7738792841865286, 0.4448506027871696),
                    Complex::new(0., 0.),
                    Complex::new(0., 0.),
                    Complex::new(0.218178442910526, 0.8469966867353856)
                ]),
                Jones::from([
                    Complex::new(0.7727769657690016, 0.4451541611407178),
                    Complex::new(0., 0.),
                    Complex::new(0., 0.),
                    Complex::new(0.21786624664314946, 0.8466270165385981)
                ]),
            ],
        ];
        let shape = (1, calsols.dim().1, sel_baselines.len());

        let mut vis_array = array![[
            // -exec p dataAsDouble
            [Jones::<f32>::from([
                Complex::new(24.25, 1.),
                Complex::new(85.5, 81.75),
                Complex::new(35.25, -2.),
                Complex::new(154.5, 9.625)
            ])],
            [Jones::<f32>::from([
                Complex::new(58.25, -67.),
                Complex::new(3.875, -12.375),
                Complex::new(-36., 75.75),
                Complex::new(17.375, 75.625)
            ])],
        ]];
        let exp_vis_array = array![[
            // -exec p dataAsDouble
            [Jones::<f32>::from([
                Complex::new(7.8246384, 17.68882),
                Complex::new(43.610638, 81.43078),
                Complex::new(7.043186, 29.182451),
                Complex::new(102.209915, 78.65481)
            ])],
            [Jones::<f32>::from([
                Complex::new(68.32589, 18.026802),
                Complex::new(5.8807054, -8.232894),
                Complex::new(-68.7944, -18.519669),
                Complex::new(-23.242767, 60.28708)
            ])],
        ]];
        let mut flag_array = Array3::from_shape_fn(shape, |_| false);
        let mut weight_array = Array3::from_shape_fn(shape, |_| 1_f32);
        let emulate_cotter = false;
        apply_di_calsol(
            calsols.view(),
            vis_array.view_mut(),
            weight_array.view_mut(),
            flag_array.view_mut(),
            &sel_baselines,
            &[],
            emulate_cotter,
        )
        .unwrap();

        compare_jones!(vis_array[(0, 0, 0)], exp_vis_array[(0, 0, 0)]);
        compare_jones!(vis_array[(0, 1, 0)], exp_vis_array[(0, 1, 0)]);
    }

    /// Tiles with all-NaN calibration solutions flag all baselines involving that tile.
    #[test]
    fn test_apply_calsols_nan_flagged_tile() {
        let sel_baselines = vec![(0, 0), (0, 1), (1, 1)];
        let num_times = 1;
        let num_chans = 2;

        let mut calsols = Array2::from_shape_fn((2, num_chans), |(_, _)| Jones::identity());
        for chan in 0..num_chans {
            calsols[(1, chan)] = Jones::nan();
        }

        let shape = (num_times, num_chans, sel_baselines.len());
        let mut vis_array = Array3::from_shape_fn(shape, |(_, _, bl)| {
            Jones::<f32>::identity() * (bl + 1) as f32
        });
        let orig_vis_array = vis_array.clone();
        let mut flag_array = Array3::from_shape_fn(shape, |_| false);
        let mut weight_array = Array3::from_shape_fn(shape, |_| 1_f32);
        let emulate_cotter = false;

        apply_di_calsol(
            calsols.view(),
            vis_array.view_mut(),
            weight_array.view_mut(),
            flag_array.view_mut(),
            &sel_baselines,
            &[false, true],
            emulate_cotter,
        )
        .unwrap();

        // auto on unflagged tile 0 is calibrated
        compare_jones!(
            vis_array[(0, 0, 0)],
            calsols[(0, 0)] * (Jones::<f64>::identity() * 1.) * calsols[(0, 0)].h()
        );
        assert!(!flag_array[(0, 0, 0)]);

        // baselines involving all-NaN tile 1 are flagged without modifying visibilities
        for chan in 0..num_chans {
            assert!(flag_array[(0, chan, 1)]);
            assert!(flag_array[(0, chan, 2)]);
            assert_eq!(weight_array[(0, chan, 1)], -1.);
            assert_eq!(weight_array[(0, chan, 2)], -1.);
            compare_jones!(vis_array[(0, chan, 1)], orig_vis_array[(0, chan, 1)]);
            compare_jones!(vis_array[(0, chan, 2)], orig_vis_array[(0, chan, 2)]);
        }
    }

    /// Partially NaN calibration solutions flag only affected channels.
    #[test]
    fn test_apply_calsols_partial_nan_channel() {
        let sel_baselines = vec![(0, 1)];
        let num_times = 1;
        let num_chans = 2;

        let mut calsols = Array2::from_shape_fn((2, num_chans), |(_, _)| Jones::identity() * 2.);
        calsols[(1, 1)] = Jones::nan();

        let shape = (num_times, num_chans, sel_baselines.len());
        let mut vis_array = Array3::from_shape_fn(shape, |_| Jones::<f32>::identity());
        let orig_vis_array = vis_array.clone();
        let mut flag_array = Array3::from_shape_fn(shape, |_| false);
        let mut weight_array = Array3::from_shape_fn(shape, |_| 1_f32);
        let emulate_cotter = false;

        apply_di_calsol(
            calsols.view(),
            vis_array.view_mut(),
            weight_array.view_mut(),
            flag_array.view_mut(),
            &sel_baselines,
            &[],
            emulate_cotter,
        )
        .unwrap();

        assert!(!flag_array[(0, 0, 0)]);
        assert!(flag_array[(0, 1, 0)]);
        compare_jones!(vis_array[(0, 0, 0)], orig_vis_array[(0, 0, 0)] * 4.);
        compare_jones!(vis_array[(0, 1, 0)], orig_vis_array[(0, 1, 0)]);
    }

    /// User-flagged tiles skip calibration even when solutions are valid.
    #[test]
    fn test_apply_calsols_user_flagged_tile() {
        let sel_baselines = vec![(0, 0), (0, 1)];
        let num_times = 1;

        let calsols = Array2::from_shape_fn((2, 1), |(_, _)| Jones::identity() * 2.);
        let shape = (num_times, calsols.dim().1, sel_baselines.len());
        let mut vis_array = Array3::from_shape_fn(shape, |(_, _, bl)| {
            Jones::<f32>::identity() * (bl + 1) as f32
        });
        let orig_vis_array = vis_array.clone();
        let mut flag_array = Array3::from_shape_fn(shape, |_| false);
        let mut weight_array = Array3::from_shape_fn(shape, |_| 1_f32);
        let emulate_cotter = false;

        apply_di_calsol(
            calsols.view(),
            vis_array.view_mut(),
            weight_array.view_mut(),
            flag_array.view_mut(),
            &sel_baselines,
            &[false, true],
            emulate_cotter,
        )
        .unwrap();

        compare_jones!(
            vis_array[(0, 0, 0)],
            calsols[(0, 0)] * (Jones::<f64>::identity() * 1.) * calsols[(0, 0)].h()
        );
        assert!(!flag_array[(0, 0, 0)]);

        assert!(flag_array[(0, 0, 1)]);
        assert_eq!(weight_array[(0, 0, 1)], -1.);
        compare_jones!(vis_array[(0, 0, 1)], orig_vis_array[(0, 0, 1)]);
    }

    #[test]
    fn test_get_calsol_nan_flagged_tiles() {
        let mut calsols = Array2::from_shape_fn((3, 2), |(_, _)| Jones::identity());
        calsols[(1, 0)] = Jones::nan();
        calsols[(1, 1)] = Jones::nan();
        calsols[(2, 0)] = Jones::from([
            Complex::new(f64::NAN, 0.),
            Complex::new(0., 0.),
            Complex::new(0., 0.),
            Complex::new(0., 0.),
        ]);

        let flagged = get_calsol_nan_flagged_tiles(calsols.view());
        assert_eq!(flagged, vec![false, true, false]);
    }
}
