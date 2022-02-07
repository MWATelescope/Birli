//! Calibrating visibilities.

use crate::ndarray::{ArrayView2, ArrayViewMut3, Axis};
use itertools::izip;
use marlu::{Jones, MarluVisContext};
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

/// apply a direction independent calibration solution for a single timeblock to the given
/// visibility data
///
/// # Errors
///
/// di_jones should have the same number of channels as vis_array, flag_array, weight_array etc.
pub fn apply_di_calsol(
    // a two dimensional array of jones matrix calibration solutions with
    // dimensions `[tile][channel]`
    di_jones: ArrayView2<Jones<f64>>,
    // dimensions `[timestep][channel][baselines]`
    mut vis_array: ArrayViewMut3<Jones<f32>>,
    // dimensions `[timestep][channel][baselines]`
    mut weight_array: ArrayViewMut3<f32>,
    // dimensions `[timestep][channel][baselines]`
    // todo: setting both flags and weights is redundant, but it's not clear how to rip this out
    mut flag_array: ArrayViewMut3<bool>,
    context: MarluVisContext,
) -> Result<(), CalibrationError> {
    let di_dims = di_jones.dim();
    let vis_dims = vis_array.dim();
    let weight_dims = weight_array.dim();
    let flag_dims = flag_array.dim();
    if weight_dims != vis_dims {
        return Err(CalibrationError::BadArrayShape {
            argument: "weight_array".into(),
            function: "apply_di_calsol".into(),
            expected: format!("{:?}", vis_dims),
            received: format!("{:?}", weight_dims),
        });
    }
    if flag_dims != vis_dims {
        return Err(CalibrationError::BadArrayShape {
            argument: "flag_array".into(),
            function: "apply_di_calsol".into(),
            expected: format!("{:?}", vis_dims),
            received: format!("{:?}", flag_dims),
        });
    }

    if di_dims.1 != vis_dims.1 {
        return Err(CalibrationError::ChannelSizeMismatch {
            calsol_chans: di_dims.1,
            data_chans: vis_dims.1,
        });
    }

    // TODO: cotter allows for n * di_dims.1 = vis_dims.1.

    // time axis
    for (mut vis_array, mut weight_array, mut flag_array) in izip!(
        vis_array.axis_iter_mut(Axis(0)),
        weight_array.axis_iter_mut(Axis(0)),
        flag_array.axis_iter_mut(Axis(0)),
    ) {
        // baseline axis
        for ((ant1_idx, ant2_idx), mut vis_array, mut weight_array, mut flag_array) in izip!(
            context.sel_baselines.clone().into_iter(),
            vis_array.axis_iter_mut(Axis(1)),
            weight_array.axis_iter_mut(Axis(1)),
            flag_array.axis_iter_mut(Axis(1)),
        ) {
            // channel axis
            for (&sol1, &sol2, vis, weight, flag) in izip!(
                di_jones.index_axis(Axis(0), ant1_idx),
                di_jones.index_axis(Axis(0), ant2_idx),
                vis_array.iter_mut(),
                weight_array.iter_mut(),
                flag_array.iter_mut(),
            ) {
                // apply the calibration solution

                // promote and divide by promoted weight
                let vis_f64 = Jones::<f64>::from(*vis) / (*weight as f64);

                // demote J1 * D * J2^H
                *vis = Jones::<f32>::from(sol1 * vis_f64 * sol2.h());

                // if the data now contains a NaN, flag it
                // todo: not sure about this
                if vis.any_nan() {
                    *flag = true;
                    if *weight > 0. {
                        *weight = -*weight;
                    }
                }
            }
        }
    }

    Ok(())
}
