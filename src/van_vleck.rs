//! Van Vleck corrections for MWA data.
//!
//! Derived from pyuvdata [mwa_corr_fits.py](https://github.com/RadioAstronomySoftwareGroup/pyuvdata/blob/main/src/pyuvdata/uvdata/mwa_corr_fits.py).
//! Original license: BSD 2-Clause Simplified
//!
//! ```txt
//! Copyright (c) 2018, Radio Astronomy Software Group
//! All rights reserved.
//!
//! Redistribution and use in source and binary forms, with or without
//! modification, are permitted provided that the following conditions are met:
//!
//! * Redistributions of source code must retain the above copyright notice, this
//!   list of conditions and the following disclaimer.
//!
//! * Redistributions in binary form must reproduce the above copyright notice,
//!   this list of conditions and the following disclaimer in the documentation
//!   and/or other materials provided with the distribution.
//!
//! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
//! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//! ```
//!
#![allow(mixed_script_confusables)]

use crate::{
    ndarray::{parallel::prelude::*, prelude::*},
    Jones,
};
use errorfunctions::RealErrorFunctions;
use itertools::{zip_eq, Itertools};
use marlu::{io::error::BadArrayShape, jones, mwalib::CorrelatorContext, rayon};
use rayon::prelude::*;
use std::{collections::HashMap, f64::consts::PI};
use thiserror::Error;

#[derive(Error, Debug)]
/// Error for Passband Corrections
pub enum VanVleckCorrection {
    #[error(transparent)]
    /// Error for bad array shape in provided argument
    BadArrayShape(#[from] BadArrayShape),

    // bad number of samples
    #[error("invalid number of correlator samples {nsamples}, check metadata")]
    BadNSamples { nsamples: u32 },
}

/// Perform Van Vleck corrections given an observation's
/// [`mwalib::CorrelatorContext`](crate::mwalib::CorrelatorContext) and an
/// [`ndarray::Array3`](crate::ndarray::Array3) of [`marlu::Jones`] visibilities (time, frequency, baseline)
///
/// Credit: <https://github.com/RadioAstronomySoftwareGroup/pyuvdata/blob/3b412bc5397dd2a87a3bb6239e9c584b52306603/src/pyuvdata/uvdata/mwa_corr_fits.py#L775>
///
/// Original license: BSD 2-Clause Simplified
///
/// # Examples
///
/// ```rust
/// use birli::{correct_cable_lengths, mwalib::CorrelatorContext, VisSelection, io::read_mwalib};
///
/// // define our input files
/// let metafits_path = "tests/data/1297526432_mwax/1297526432.metafits";
/// let gpufits_paths = vec![
///     "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits",
/// ];
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
/// let ant_pairs = vis_sel.get_ant_pairs(&corr_ctx.metafits_context);
///
/// van_vleck_correction(&corr_ctx, jones_array.view_mut(), &ant_pairs);
/// ```
///
/// # Errors
/// - [`VanVleckCorrection::BadArrayShape`] - If the shape of the provided `ant_pairs` argument is incorrect.
pub fn correct_van_vleck(
    corr_ctx: &CorrelatorContext,
    mut jones_array: ArrayViewMut3<Jones<f32>>,
    // baseline_idxs: &[usize],
    ant_pairs: &[(usize, usize)],
    // TODO: flagged_ant_idxs: &[usize],
) -> Result<(), VanVleckCorrection> {
    let vis_dims = jones_array.dim();
    if vis_dims.2 != ant_pairs.len() {
        return Err(VanVleckCorrection::BadArrayShape(BadArrayShape {
            argument: "ant_pairs",
            function: "correct_van_vleck",
            expected: format!("vis_dims.2={}", vis_dims.2,),
            received: format!("{:?}", ant_pairs.len()),
        }));
    }

    // scale the data by 1/2N, where N is the number of samples per fine channel
    //  = channel width (Hz)
    //  * integration time (s)
    // circular symmetry gives a factor of two
    // e.g. 40kHz * 2s = 160_000
    let nsamples = corr_ctx.metafits_context.corr_fine_chan_width_hz
        * corr_ctx.metafits_context.corr_int_time_ms as u32
        / 2_000;

    dbg!(&nsamples);
    if nsamples < 1 {
        return Err(VanVleckCorrection::BadNSamples { nsamples });
    }

    // let mut jones_double = jones_array.mapv(|j| Jones::<f64>::from(j) / (2.0 * nsamples as f64));

    // TODO: figure out antenna flags
    // TODO: ensure not mwax

    // map between antenna index (from ant_pairs) and unflagged antenna index
    // let auto_map: HashMap<usize, usize> = ant_pairs
    //     .iter()
    //     .filter_map(|&(ant1, ant2)| {
    //         // TODO: filter flagged antennas
    //         if ant1 == ant2 {
    //             Some(ant1)
    //         } else {
    //             None
    //         }
    //     })
    //     .enumerate()
    //     .map(|(i, a)| (a, i))
    //     .collect();
    // let good_autos = ant_pairs
    //     .iter()
    //     .enumerate()
    //     .filter_map(|(i, (ant1, ant2))| {
    //         // TODO: filter flagged antennas
    //         if ant1 == ant2 {
    //             Some(ant1)
    //         } else {
    //             None
    //         }
    //     });

    // let mut auto_mask = auto_map.values().copied().collect::<Vec<usize>>();
    // auto_mask.sort_unstable();

    // ant_pair indices which are unflagged autocorrelations
    let (unflagged_auto_mask, unflagged_autos): (Vec<_>, Vec<_>) = ant_pairs
        .iter()
        .enumerate()
        .filter_map(|(i, (ant1, ant2))| {
            // TODO: filter flagged antennas
            if ant1 == ant2 {
                Some((i, ant1))
            } else {
                None
            }
        })
        .unzip();
    let n_autos = unflagged_autos.len();
    dbg!(&n_autos);

    // // mask for ant_pairs which are cross-correlations
    // let cross_mask = ant_pairs
    //     .iter()
    //     .enumerate()
    //     .filter_map(|(i, &(ant1, ant2))| {
    //         if ant1 == ant2 { None } else { Some(i) }
    //     })
    //     .collect::<Vec<usize>>();
    // let n_crosses = cross_mask.len();
    // dbg!(&n_crosses);
    // // lookup table between index in cross mask and ant pair

    // new mutable copy of jones_array, only autos
    // TODO: pyuvdata does `[timestep][baseleine][channel]` -> `[baseleine][timestep * channel]`
    // but we do `[timestep][channel][baseleine]`.
    // let mut jones_autos = jones_double.select(Axis(2), &auto_mask);
    // let autos_dim = jones_autos.dim();
    // dbg!(&jones_autos.dim());

    // partitioning of auto jones matrix:
    // - `sighat_xxr_yyr` <- xx real and yy real, for correction with with van_vleck_autos
    // - `sighat_xyr_xyi` <- xy real and xy imaginary are corrected with van_vleck_crosses using the corrected autos
    // - xx imag and yy imag are assumed to be zero
    let (sighat_xxr_yyr, _sighat_xyr_xyi): (Vec<_>, Vec<_>) = jones_array
        .select(Axis(2), &unflagged_auto_mask)
        .iter()
        .flat_map(|j| {
            let j_f64 = Jones::<f64>::from(j) / (nsamples as f64);
            [
                (j_f64[0].re.sqrt(), j_f64[1].re.sqrt()), // left: xx, right: xy
                (j_f64[3].re.sqrt(), j_f64[2].re.sqrt()), // left: yy, right: yx
            ]
        })
        .unzip();

    let sighat_xxr_yyr = jones_array
        .select(Axis(2), &unflagged_auto_mask)
        .iter()
        .flat_map(|j| {
            let j_f64 = Jones::<f64>::from(j) / (nsamples as f64);
            [
                j_f64[0].re.sqrt(), // xx
                j_f64[3].re.sqrt(), // yy
            ]
        })
        .collect_vec();
    dbg!(&sighat_xxr_yyr);

    // partitioning of auto jones matrix into vectors of [time * freq * n_autos]:
    // - `sighat_xxr`, `sighat_yyr` <- sqrt of xx real and yy real, for correction with with van_vleck_autos
    // - `auto_xyr`, `auto_xyi` <- xy real and xy imaginary are corrected with van_vleck_crosses using the corrected autos
    // - xx imag and yy imag are assumed to be zero
    // - i think yxr + i * yxi is assumed to be the the complex conjugate of xyr + i*yxi
    // let (sighat_xxr, sighat_yyr, auto_xyr, auto_xyi): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) =
    //     jones_array
    //         .select(Axis(2), &unflagged_auto_mask)
    //         .iter()
    //         .map(|j| {
    //             let j_f64 = Jones::<f64>::from(j);
    //             (
    //                 j_f64[0].re.sqrt(), // sqrt xx real
    //                 j_f64[3].re.sqrt(), // sqrt yy real
    //                 j_f64[1].re,        // xy real
    //                 j_f64[1].im,        // xy imag
    //             )
    //         })
    //         .unzip();

    // let max_sighat = sighat.iter().fold(0.0_f64, |acc, &j| acc.max(j));
    // dbg!(&max_sighat);

    // correct autos
    let sigma_xxr_yyr = Array::from(van_vleck_autos(&sighat_xxr_yyr))
        .into_shape((vis_dims.0, vis_dims.1, n_autos, 2))
        .unwrap();

    // reassign the corrected values back to jones_double
    // jones_double
    //     .axis_iter_mut(Axis(0))
    //     .zip_eq(sigma.axis_iter(Axis(0)))
    //     .for_each(|(mut j_fb, s_fb)| {
    //         j_fb.axis_iter_mut(Axis(0))
    //             .zip_eq(s_fb.axis_iter(Axis(0)))
    //             .for_each(|(mut j_b, s_b)| {
    //                 unflagged_auto_mask
    //                     .iter()
    //                     .zip_eq(s_b.axis_iter(Axis(0)))
    //                     .for_each(|(&i, s)| {
    //                         j_b[i][0].re = s[0];
    //                         j_b[i][0].im = 0.0;
    //                         j_b[i][3].re = s[1];
    //                         j_b[i][3].im = 0.0;
    //                     });
    //             });
    //     });

    jones_array
        .axis_iter_mut(Axis(2))
        .zip_eq(ant_pairs.iter())
        .for_each(|(mut j_tf, &(ant1, ant2))| {
            match (
                unflagged_autos.binary_search(&ant1),
                unflagged_autos.binary_search(&ant2),
            ) {
                // unflagged autocorrelation.
                // - update jones xx real, yy real from sigma_xxr_yyr
                // - set xx imag, yy imag to zero
                // - correct yx with yy real, xx real
                (Ok(s_idx1), _) if ant1 == ant2 => {
                    let sigma_tf = sigma_xxr_yyr.index_axis(Axis(2), s_idx1);
                    sigma_tf
                        .indexed_iter()
                        .for_each(|((t_idx, f_idx, p_idx), &s)| {
                            // map sigma_xxr_yyr pol index 0,1 to jones index 0,3 (xx, yy) by multiplying by 3
                            let jp_idx = p_idx * 3;
                            let re = (nsamples as f64) * s.powi(2);
                            println!(
                                "j@({:?},{:?}).tf[({:?},{:?})][{:?}]={:?}<-{:?}",
                                ant1,
                                ant2,
                                t_idx,
                                f_idx,
                                jp_idx,
                                j_tf[(t_idx, f_idx)][jp_idx].re,
                                re
                            );
                            j_tf[(t_idx, f_idx)][jp_idx].re = re as f32;
                            // j_tf[(t_idx, f_idx)][p_idx * 3].im = 0.0; ?
                        });
                }
                // unflagged cross correlation
                // (Ok(s_idx1), Ok(s_idx2)) => {
                //     let sigma_xxr_yyr_tf1 = sigma_xxr_yyr.index_axis(Axis(2), s_idx1).into_owned();
                //     let sigma_xxr_yyr_tf2 = sigma_xxr_yyr.index_axis(Axis(2), s_idx2).into_owned();

                //     let (khat_xxr_yyr_tf, j_xyr_xyi_tf) = j_tf
                //         .iter()
                //         .flat_map(|j| {
                //             let j_f64 = Jones::<f64>::from(j);
                //             [
                //                 (j_f64[0].re.sqrt(), j_f64[0].re.sqrt()), // left: xx, right: xy
                //                 (j_f64[3].re.sqrt(), j_f64[1].re.sqrt()), // left: yy, right: yx
                //             ]
                //         })
                //         .unzip();

                //     // TODO: correct cross xy, yx
                //     // let s_xxr_yyr_tf1 = sighat_xyr_xyi.index_axis(Axis(2), s_idx1).into_owned();
                //     // let s_xxr_yyr_tf2 = sighat_xyr_xyi.index_axis(Axis(2), s_idx2).into_owned();
                // }
                // either antenna is flagged
                _ => {}
            };
        });

    // jones_double.select(Axis(2), &cross_mask).axis_iter(Axis(2)).enumerate()

    Ok(())
}

/// Use Newton's method to solve the inverse of `sighat_vector`.
fn van_vleck_auto(s: f64) -> Option<f64> {
    let tol = 1e-10;
    let niter = 10;
    let mut guess = s;
    let mut delta = sighat(guess) - s;
    let mut count = 0;
    while delta.abs() > tol {
        guess -= delta / sighat_prime(guess);
        delta = sighat(guess) - s;
        count += 1;
        if count > niter {
            log::warn!("Van Vleck correction did not converge for sigma={s}");
            return None;
        }
    }
    Some(guess)
}

/// Apply `van_vleck_auto` over `sighats` in parallel.
pub fn van_vleck_autos(hat: &[f64]) -> Vec<f64> {
    hat.par_iter()
        .map(|&sighat| {
            van_vleck_auto(sighat).map_or(sighat, |sigma| {
                println!("sigma={sigma} <- sighat={sighat}");
                sigma
            })
        })
        .collect()
}

#[allow(clippy::doc_markdown)]
/// Generate quantized $\hat\sigma$ using Van Vleck relation.
/// $$\hat\sigma(\sigma_i) = \left[ 7^2 - \sum_{k=0}^6(2k+1){\rm erf}\left({ k+0.5 \over \sigma_i\sqrt{2} }\right) \right]^{1/2}$$
/// substitute $K=(k+0.5)$
fn sighat(sigma: f64) -> f64 {
    let mut sum: f64 = 0.0;
    for k in 0..=6 {
        let k_ = k as f64 + 0.5;
        sum += 2.0 * k_ * (k_ / (sigma * std::f64::consts::SQRT_2)).erf();
    }
    (49.0_f64 - sum).sqrt()
}

const SQRT_TAU: f64 = 2.5066282746310002;

#[allow(clippy::doc_markdown)]
/// The derivative of the $\hat σ$ function .
/// $\hat σ^\prime(σ_i) = \sum_{k=0}^6{(2k+1)(k+0.5)e^{-(k+0.5)^2/2σ_i^2}\over\sqrt{2π}σ_i^2}$
/// substitute $K=(k+0.5)^2$, $S=σ_i^2$:
/// $$ = \sum_{k=0}^6{(2K^2e^{-K^2/2S}\over\sqrt{2π}S}$$
fn sighat_prime(sigma: f64) -> f64 {
    let s_ = sigma.powi(2);
    let mut sum = 0.0;
    for k in 0..=6 {
        let k_ = (k as f64 + 0.5).powi(2);
        sum += 2.0 * k_ * (-k_ / (2.0 * s_)).exp() / (SQRT_TAU * s_);
    }
    sum / sighat(sigma)
}

#[cfg(test)]
mod vv_auto_tests {
    use float_cmp::assert_approx_eq;

    use super::*;

    const SIGHATS: [f64; 20] = [
        1.3732557118031588, // sqrr=1.17185994 , sq=1.88583125
        1.4567407971221236,
        1.58477324876463,
        1.7205649508228396,
        1.826940748902383,
        1.8929606440705524,
        1.925808271869243,
        1.932247719626032,
        1.94109505176846,
        1.9421363881046048,
        1.9405717585289137,
        1.945186366392691,
        1.9506393182749087,
        1.9506457264198438,
        1.945944500750214,
        1.9444102576359754,
        1.9511054558890455,
        1.9488121382011145,
        1.939882406229821,
        1.9340307650086646,
    ];
    const SIGMAS: [f64; 20] = [
        1.3425715134733938, // 1.15869388
        1.427852482209185,
        1.5582670082555274,
        1.6962213882104307,
        1.80413614011039,
        1.87109216839722,
        1.9044119839802796,
        1.9109450441433622,
        1.9199216944258406,
        1.9209783088033163,
        1.9193907283568603,
        1.9240731081035445,
        1.9296064755666014,
        1.9296129784329366,
        1.9248424008595775,
        1.9232855835622369,
        1.930079504724327,
        1.927752308498216,
        1.9186912731345944,
        1.9127540839654953,
    ];
    const SIGHATS_PRIME: [f64; 20] = [
        0.9776527939739493,
        0.9801533937619006,
        0.9831603111567421,
        0.9852926971585774,
        0.9860246973943583,
        0.9859337319315375,
        0.9857012789829934,
        0.9856397675064799,
        0.9855463948958575,
        0.9855347215188939,
        0.985552206321147,
        0.9854996947568739,
        0.9854339402996796,
        0.9854338606488524,
        0.9854907938706629,
        0.9855087264060677,
        0.9854281317512257,
        0.9854564230798217,
        0.9855598061732679,
        0.9856217802550746,
    ];

    #[test]
    fn test_sighat() {
        let result = SIGMAS.iter().map(|&s| sighat(s)).collect_vec();
        for (&r, &s) in result.iter().zip(SIGHATS.iter()) {
            assert_approx_eq!(f64, r, s, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_sighat_prime() {
        let result = SIGMAS.iter().map(|&s| sighat_prime(s)).collect_vec();
        for (&r, &s) in result.iter().zip(SIGHATS_PRIME.iter()) {
            assert_approx_eq!(f64, r, s, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_van_vleck_autos() {
        let result = van_vleck_autos(&SIGHATS);
        for (&r, &s) in result.iter().zip(SIGMAS.iter()) {
            assert_approx_eq!(f64, r, s, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_correct_van_vleck_autos_bad_nsamples() {
        let vis_dims = (1, 1, 3);
        let mut corr_ctx = CorrelatorContext::new(
            "tests/data/1297526432_mwax/1297526432.metafits",
            &["tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits"],
        )
        .unwrap();
        corr_ctx.metafits_context.corr_fine_chan_width_hz = 1;
        corr_ctx.metafits_context.corr_int_time_ms = 1;

        let mut jones_array = Array3::<Jones<f32>>::zeros(vis_dims);

        let ant_pairs = vec![(0, 0), (0, 1), (1, 1)];

        // assert error is BadNSamples
        assert!(correct_van_vleck(&corr_ctx, jones_array.view_mut(), &ant_pairs).is_err());
    }
    #[test]
    fn test_correct_van_vleck_autos_good() {
        let vis_dims = (1, 1, 3);
        let mut corr_ctx = CorrelatorContext::new(
            "tests/data/1297526432_mwax/1297526432.metafits",
            &["tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits"],
        )
        .unwrap();
        corr_ctx.metafits_context.corr_fine_chan_width_hz = 1;
        corr_ctx.metafits_context.corr_int_time_ms = 2_000;

        let mut jones_array = Array3::<Jones<f32>>::zeros(vis_dims);

        jones_array[(0, 0, 0)][0].re = SIGHATS[0].powi(2) as f32;
        jones_array[(0, 0, 0)][1].re = SIGHATS[1].powi(2) as f32;
        jones_array[(0, 0, 0)][2].re = SIGHATS[2].powi(2) as f32;
        jones_array[(0, 0, 0)][3].re = SIGHATS[3].powi(2) as f32;

        jones_array[(0, 0, 1)][0].re = SIGHATS[4].powi(2) as f32;
        jones_array[(0, 0, 1)][1].re = SIGHATS[5].powi(2) as f32;
        jones_array[(0, 0, 1)][2].re = SIGHATS[6].powi(2) as f32;
        jones_array[(0, 0, 1)][3].re = SIGHATS[7].powi(2) as f32;

        jones_array[(0, 0, 2)][0].re = SIGHATS[8].powi(2) as f32;
        jones_array[(0, 0, 2)][1].re = SIGHATS[9].powi(2) as f32;
        jones_array[(0, 0, 2)][2].re = SIGHATS[10].powi(2) as f32;
        jones_array[(0, 0, 2)][3].re = SIGHATS[11].powi(2) as f32;

        let ant_pairs = vec![(0, 0), (0, 1), (1, 1)];

        correct_van_vleck(&corr_ctx, jones_array.view_mut(), &ant_pairs).unwrap();

        assert_approx_eq!(f32, jones_array[(0, 0, 0)][0].re, SIGMAS[0].powi(2) as f32);
        // assert_approx_eq!(f32, jones_array[(0, 0, 0)][1].re, SIGMAS[1].powi(2) as f32);
        // assert_approx_eq!(f32, jones_array[(0, 0, 0)][2].re, SIGMAS[2].powi(2) as f32);
        assert_approx_eq!(f32, jones_array[(0, 0, 0)][3].re, SIGMAS[3].powi(2) as f32);

        // assert_approx_eq!(f32, jones_array[(0, 0, 1)][0].re, SIGMAS[4].powi(2) as f32);
        // assert_approx_eq!(f32, jones_array[(0, 0, 1)][1].re, SIGMAS[5].powi(2) as f32);
        // assert_approx_eq!(f32, jones_array[(0, 0, 1)][2].re, SIGMAS[6].powi(2) as f32);
        // assert_approx_eq!(f32, jones_array[(0, 0, 1)][3].re, SIGMAS[7].powi(2) as f32);

        assert_approx_eq!(f32, jones_array[(0, 0, 2)][0].re, SIGMAS[8].powi(2) as f32);
        // assert_approx_eq!(f32, jones_array[(0, 0, 2)][1].re, SIGMAS[9].powi(2) as f32);
        // assert_approx_eq!(f32, jones_array[(0, 0, 2)][2].re, SIGMAS[10].powi(2) as f32);
        assert_approx_eq!(f32, jones_array[(0, 0, 2)][3].re, SIGMAS[11].powi(2) as f32);
    }
}

// /////// //
// CROSSES //
// /////// //

// Evaluate the bivariate normal probability density function
// over a grid of $x_i$ and $y_k$ where $ρ$ is the correlation coefficient.
// $σ$ are the standard deviations, which are used to normalize $x,y$
// a lot of MWA specific stuff involving symmetries of the x grid.
//
// takes a slice of ρ because we're integrating over ρ.
//
// assumes σ_x,σ_y > 0 and -1 ≤ ρ ≤ 1
//
// write back result to &mut r_
//
// Definition of bivariate normal pdf
// $$
// f(x,y,σ_x,σ_y,ρ) = {1 \over 2π σ_x σ_y\sqrt{1-ρ^2}} \exp\left[
//     -{1 \over 2(1-ρ^2)}\left(
//         {x^2 \over σ_x^2}
//         + {y^2 \over σ_y^2}
//         - 2ρ {xy \over σ_x σ_y}
//     \right)
// \right]
// $$
// From Price's theorem
// $$
// {d\hat κ \over dρ} = σ_x σ_y \sum_{i=1}^{n-1} \sum_{k=1}^{n-1} Δ h_i Δ h_k f(a_i,a_k,σ,ρ)
// $$
//
// for legacy MWA RRI Rx, 4-bit symmetric quantisation pattern, change summation $\sum_{i=-7}^6$
//
// $$\begin{align}
// \text{} a_i&=(i+{1\over 2}), Δ h_i=1 \\
// \text{define }
// 	x_i &= {i+{1\over 2} \over σ_x},
// 	y_k = {k+{1\over 2} \over σ_y},
// 	r = 1-ρ^2 \\
// 	{d\hat κ \over dρ} &= {1 \over 2π \sqrt{r}} \sum_{i=-7}^{6} \sum_{k=-7}^{6} \exp\left[
// 	    -{x_i^2 + y_k^2 - 2ρ x_i y_k \over 2r}
// 	\right]
// \end{align}$$
// also make use of symmetry
// $$
// x_i = {i+{1\over 2} \over σ_x} = -{(-i-1)+{1\over 2} \over σ_x} = -x_{-i-1}
// $$
// and
// $$\begin{align}
// \text{if } f(a,b) &= a^2 + b^2 - 2 a b \\
// 	f(-a,-b) &= a^2 + b^2 - 2 a b = f(a,b)\\
// 	f(-a,b) &= a^2 + b^2 + 2 a b = f(a,-b)\\
// \end{align}$$
// so
// $$\begin{align}
// \sum_{i=-7}^{6} \sum_{k=-7}^{6} f(x_i,y_k)
// 	&= \sum_{i=-7}^{6} \sum_{k=0}^{6} \left[
// 		f(x_i,y_k) + f(x_i,-y_k)
// 	\right] \\
// 	&= \sum_{i=0}^{6} \sum_{k=0}^{6} \left[
// 		f(x_i,y_k) + f(x_i,-y_k) + f(-x_i,y_k) + f(-x_i,-y_k)
// 	\right] \\
// 	&= \sum_{i=0}^{6} \sum_{k=0}^{6} \left[
// 		2f(x_i,y_k) + 2f(x_i,-y_k)
// 	\right] \\
// \end{align}$$
// finally
// $$\begin{align}
// \text{define }
// 	α_{ik}(ρ) &= {x_i^2 + y_i^2 \over 2(1-ρ^2)},
// 	β_{ik}(ρ) = {x_i y_i \over 1-ρ^2}\\
// 	{d\hat κ \over dρ} &= {1 \over π \sqrt{r}} \sum_{i=0}^{6} \sum_{k=0}^{6} \left[
// 		\exp(-α_{ik} + β_{ik}ρ)
// 		+ \exp(-α_{ik} - β_{ik}ρ)
// 	\right] \\
// 	&= {1 \over π \sqrt{r}} \sum_{i=0}^{6} \sum_{k=0}^{6} \left[
// 		e^{-α_{ik}}\cdot e^{β_{ik}ρ}
// 		+ e^{-α_{ik}}\cdot e^{-β_{ik}ρ}
// 	\right] \\
// 	&= {1 \over π \sqrt{r}} \sum_{i=0}^{6} \sum_{k=0}^{6} e^{-α_{ik}} \left(
// 		e^{β_{ik}ρ} + e^{-β_{ik}ρ}
// 	\right) \\
// 	&= {1 \over π \sqrt{r}} \sum_{i=0}^{6} \sum_{k=0}^{6} {2 \cosh (β_{ik}ρ) \over e^{α_{ik}}} \\
// 	&= {1 \over π \sqrt{r}} \sum_{i=0}^{6} \sum_{k=0}^{6} {2 \cosh ({x_i y_i ρ \over r}) \over
// 		e^{(x_i^2 + y_i^2) / 2r}
// 	} \\
// \end{align}$$
//
// https://www.geogebra.org/m/pO4JcWPz
fn pdf(x_: &[f64], y_: &[f64], ρ_: &[f64], r_: &mut [f64]) {
    // iterate over ρ_, and store result in r_
    for (ρ, r) in zip_eq(ρ_.iter(), r_.iter_mut()) {
        *r = 0.0;
        // (1-ρ^2) should be > 0
        let d: f64 = (1.0 - ρ.powi(2)).max(1e-20);
        // if d.abs() <= 1e-20 {
        //     println!("d is small ρ={}, d={}", ρ, d);
        // }
        for (x, y) in x_.iter().cartesian_product(y_.iter()) {
            let numer = 2.0 * (*ρ * x * y / d).cosh();
            // if f64::is_infinite(numer) {
            //     println!("numer is infinite: x={}, y={}, ρ={}, d={}", x, y, ρ, d);
            // }
            let denom = (((x).powi(2) + (y).powi(2)) / (2.0 * d)).exp();
            // if denom.abs() < 1e-20 {
            //     println!("denom is zero: x={}, y={}, ρ={}, d={}", x, y, ρ, d);
            // }
            *r += numer / denom;
        }
        *r /= PI * d.sqrt();
    }
}

// Use Simpson's rule to evaluate $\int_a^b \sum_i \sum_k f(x_i,y_k,ρ) dρ$ with `n` intervals
// in this case, `f` is the bivariate normal pdf, evaluated at $x_i,y_k$ already normalized by σ
#[allow(clippy::too_many_arguments)]
fn simpsons_rule<F>(f: F, a: f64, b: f64, n: usize, x_: &[f64], y_: &[f64]) -> f64
where
    F: Fn(&[f64], &[f64], &[f64], &mut [f64]),
{
    debug_assert!(n % 2 == 0, "n must be even: {n:?}");
    debug_assert!(b < 1.0, "b (which is ρ) must be < 1: {b:?}");
    debug_assert!(b > -1.0, "b (which is ρ) must be > -1: {b:?}");
    let h = (b - a) / n as f64;
    let mut r_ = vec![0.0; n + 1]; // result
    let ρ_ = (0..=n).map(|i| a + i as f64 * h).collect_vec(); // integration points

    // compute the pdf at ρ and store result in r
    f(x_, y_, &ρ_, &mut r_);

    let mut sum = r_[0] + r_[n];

    for (i, r) in r_[1..n].iter().enumerate() {
        sum += if i % 2 == 1 { 2.0 * *r } else { 4.0 * *r };
    }

    sum * h / 3.0
}

// generate $\hat κ$ from rho using the Van Vleck relation
//
// sigma is encoded in the x and y grid
fn corrcorrect_simp(rho: f64, x_: &[f64], y_: &[f64]) -> f64 {
    // number of intervals in Simpson integral
    let n = 10;
    simpsons_rule(pdf, 0.0, rho, n, x_, y_)
}

// calculate $d \hat κ / dρ$
//
// sigma is encoded in the x and y grid
fn corrcorrect_prime(rho: f64, x_: &[f64], y_: &[f64]) -> f64 {
    let mut r_ = vec![0.0; 1];
    pdf(x_, y_, &[rho], &mut r_);
    r_[0]
}

// solve a single van vleck cross correlation using the Van Vleck relation for legacy MWA
fn van_vleck_cross_int(khat: f64, sigma_x: f64, sigma_y: f64) -> Option<f64> {
    debug_assert!(sigma_x > 0.0, "σ_x must be > 0: {sigma_x:?}");
    debug_assert!(sigma_y > 0.0, "σ_y must be > 0: {sigma_y:?}");
    let sign = khat.signum();
    let khat = khat.abs();

    // from legacy MWA RRI Rx symmetric quantisation pattern
    let x_ = (0..7).map(|i| (i as f64 + 0.5) / sigma_x).collect_vec();
    let y_ = (0..7).map(|i| (i as f64 + 0.5) / sigma_y).collect_vec();

    let tol = 1e-10; // it's 1e-10 for the autos
    let niter = 10;
    let mut guess = khat / (sigma_x * sigma_y);
    debug_assert!(guess >= 0.0, "|ρ| must be >= 0: {guess:?}");
    debug_assert!(guess < 1.0, "|ρ| must be < 1: {guess:?}");
    let mut delta = corrcorrect_simp(guess, &x_, &y_) - khat;
    let mut count = 0;
    while delta.abs() > tol {
        guess -= delta / corrcorrect_prime(guess, &x_, &y_);
        delta = corrcorrect_simp(guess, &x_, &y_) - khat;
        count += 1;
        if count > niter {
            log::warn!("Van Vleck correction did not converge for khat={khat}");
            return None;
        }
    }
    Some(sign * guess * sigma_x * sigma_y)
}

/// Apply `van_vleck_cross_int` over $\hat κ, σ_1, σ_2$ in parallel.
pub fn van_vleck_crosses_int(k_: &[f64], σ1_: &[f64], σ2_: &[f64]) -> Vec<f64> {
    k_.par_iter()
        .zip_eq(σ1_.par_iter())
        .zip_eq(σ2_.par_iter())
        .map(|((k, σ1), σ2)| van_vleck_cross_int(*k, *σ1, *σ2).unwrap_or(*k))
        .collect()
}
// pub fn van_vleck_crosses_int(*, k_arr, sig1_arr, sig2_arr, cheby_approx):

#[cfg(test)]
mod vv_cross_tests {
    use super::*;
    use float_cmp::assert_approx_eq;

    #[test]
    fn test_simpsons_identity() {
        fn test_function(_: &[f64], _: &[f64], ρ_: &[f64], r_: &mut [f64]) {
            for (ρ, r) in zip_eq(ρ_.iter(), r_.iter_mut()) {
                *r = *ρ;
            }
        }

        let a = 0.0;
        let b = 0.9;
        let n = 4; // number of intervals
        let x_ = &[0.0; 1]; // grid x
        let y_ = &[0.0; 1]; // grid y

        let result = simpsons_rule(test_function, a, b, n, x_, y_);
        let expected = b.powi(2) / 2.0; // Analytical result of $∫_0^b x dx$

        assert_approx_eq!(f64, result, expected);
    }

    #[test]
    // $$\begin{align}
    // x_i = y_k = 0 &\implies
    //     α_{ik}(ρ) = 0 =
    //     β_{ik}(ρ)\\
    //  &\implies {1 \over π \sqrt{r}} \sum_{i=0}^{6} \sum_{k=0}^{6} {2 \cosh (β_{ik}ρ) \over e^{α_{ik}}
    //     } = {2 \over π \sqrt{r}} \\
    // &\implies ∫_0^ρ {1\over π\sqrt{1−ζ}^2} dζ = {2\over π} ​\arcsin(ρ)
    // \end{align}$$
    fn test_simpsons_pdf() {
        let a = 0.0;
        let b = 0.5;
        let n = 10; // number of intervals
        let x_ = &[0.0; 1]; // grid x
        let y_ = &[0.0; 1]; // grid y

        let result = simpsons_rule(pdf, a, b, n, x_, y_);
        let expected = (b).asin() * 2.0 / PI; // Analytical result of $∫_0^b 2/π√(1−ζ^2) dζ$

        assert_approx_eq!(f64, result, expected, epsilon = 1e-6);
    }

    #[test]
    fn test_corrcorrect_simp() {
        let sigma_x = 1.03637188;
        let sigma_y = 0.98278517;
        let rho = 0.01021078;
        let expexted_khat = 0.01040000;
        let x_ = (0..7).map(|i| (i as f64 + 0.5) / sigma_x).collect_vec();
        let y_ = (0..7).map(|i| (i as f64 + 0.5) / sigma_y).collect_vec();
        // println!("{:?}", x_);
        let khat = corrcorrect_simp(rho, &x_, &y_);
        assert_approx_eq!(f64, khat, expexted_khat, epsilon = 1e-6);
    }

    #[test]
    fn test_corrcorrect_prime() {
        let sigma_x = 1.03637188;
        let sigma_y = 0.98278517;
        let rho = 0.01021078;
        let expexted_result = 1.0185308989;
        let x_ = (0..7).map(|i| (i as f64 + 0.5) / sigma_x).collect_vec();
        let y_ = (0..7).map(|i| (i as f64 + 0.5) / sigma_y).collect_vec();
        let khat = corrcorrect_prime(rho, &x_, &y_);
        assert_approx_eq!(f64, khat, expexted_result, epsilon = 1e-6);
    }

    // test data:
    // ```python
    // test_dir = environ['BIRLI_TEST_DIR']
    // obsid = 1090701368
    // vis_dir = path.join(test_dir, f'{obsid}/raw')
    // metafits = path.join(vis_dir, f'{obsid}.metafits')
    // raw_path = path.join(vis_dir, f'{obsid}_20140729203555_gpubox01_00.fits')
    // out_dir = '.'

    // uvd_raw = UVData()
    // uvd_vvcheby = UVData()
    // uvd_vvnocheby = UVData()

    // nocorrect = {
    //     'use_aoflagger_flags': False,
    //     'remove_dig_gains': False,
    //     'remove_coarse_band': False,
    //     'correct_cable_len': False,
    //     'flag_small_auto_ants': False,
    //     'phase_to_pointing_center': False,
    //     'propagate_coarse_flags': False,
    //     'flag_init': False,
    //     'remove_flagged_ants': False,
    //     'fix_autos': False,
    // }
    // nocheck = {
    //     'run_check': False,
    //     'check_extra': False,
    //     'run_check_acceptability': False,
    //     'strict_uvw_antpos_check': False,
    //     'check_autos': False,
    // }
    // nofix = {
    //     'run_check': False,
    //     'check_extra': False,
    //     'run_check_acceptability': False,
    //     'strict_uvw_antpos_check': False,
    //     'check_autos': False,
    //     'fix_autos': False,
    // }

    // uvd_raw.read_mwa_corr_fits( [metafits, raw_path], correct_van_vleck=False, **nocorrect )
    // uvd_vvcheby.read_mwa_corr_fits([metafits, raw_path], correct_van_vleck=True, cheby_approx=True, **nocorrect)
    // uvd_vvnocheby.read_mwa_corr_fits([metafits, raw_path], correct_van_vleck=True, cheby_approx=False, **nocorrect)

    // ant1 = uvd_raw.get_ants()[0]
    // ant2 = uvd_raw.get_ants()[1]
    // nsamples = (uvd_raw.channel_width * uvd_raw.integration_time[0]) * 2
    // print(f"{nsamples=} {uvd_raw.channel_width=} {uvd_raw.integration_time=}")
    // data11_vvnocheby = np.sqrt(uvd_vvnocheby.get_data(*(ant1, ant1), 'xx')[0].real / nsamples)
    // data22_vvnocheby = np.sqrt(uvd_vvnocheby.get_data(*(ant2, ant2), 'xx')[0].real / nsamples)
    // data12_raw = uvd_raw.get_data(*(ant1, ant2), 'xx')[0].real / nsamples
    // data12_vvnocheby = uvd_vvnocheby.get_data(*(ant1, ant2), 'xx')[0].real / nsamples
    // print(f"{data11_vvnocheby=}")
    // print(f"{data22_vvnocheby=}")
    // print(f"{data12_raw=}")
    // print(f"{data12_vvnocheby=}")

    // ```

    #[rustfmt::skip]
    const SIGMAS1: [f64; 480] = [
        1.0363718780, 1.0314391258, 1.0454504621, 1.0310997382, 1.0505078145,
        1.0533122368, 1.0522198764, 1.0393106702, 1.0414253065, 1.0616810575,
        1.0527424511, 1.0355030998, 1.0464543319, 1.0460720189, 1.0587571336,
        1.0558251125, 1.0379627495, 1.0465498883, 1.0520297849, 1.0491266218,
        1.0444934987, 1.0470275396, 1.0372399286, 1.0523624225, 1.0590876583,
        1.0381554168, 1.0473617659, 1.0542137678, 1.0514117501, 1.0558251125,
        1.0473617659, 1.0525049493, 1.0417613302, 1.0495554623, 1.0547353546,
        1.0536919228, 1.0610686443, 1.0365166036, 1.0608801383, 1.0554461939,
        1.0662863910, 1.0755773647, 1.0459286155, 1.0561565547, 1.0581902797,
        1.0475049731, 1.0431043423, 1.0364201221, 1.0546405398, 1.0463109809,
        1.0516970421, 1.0497936314, 1.0659580985, 1.0428166996, 1.0617752436,
        1.0466932066, 1.0463587667, 1.0386850674, 1.0598427562, 1.0679263400,
        1.0557303956, 1.0479344772, 1.0545931292, 1.0583320216, 1.0454026348,
        1.0700778792, 1.0546405398, 1.0624813729, 1.0597012163, 1.0483638054,
        1.0462631928, 1.0415693299, 1.0455461103, 1.0547353546, 1.0428166996,
        1.0537393739, 1.0591820750, 1.0340051595, 1.0496507364, 1.0608330066,
        1.0692832495, 1.0465498883, 1.0543086209, 1.0512215125, 1.0533597051,
        1.0546879482, 1.0508409339, 1.0476481607, 1.0407048905, 1.0403204647,
        1.0544983016, 1.0562985696, 1.0497936314, 1.0504602174, 1.0584265057,
        1.0551619156, 1.0542137678, 1.0581902797, 1.0799382701, 1.0566298633,
        1.0621048290, 1.0602201035, 1.0500317464, 1.0506030023, 1.0542137678,
        1.0561092122, 1.0587571336, 1.0469797843, 1.0506981814, 1.0760421309,
        1.0525524539, 1.0555883043, 1.0517445832, 1.0541189061, 1.0622931177,
        1.0549249585, 1.0402243360, 1.0568191274, 1.0499365069, 1.0363236317,
        1.0363718780, 1.0314391258, 1.0454504621, 1.0310997382, 1.0505078145,
        1.0533122368, 1.0522198764, 1.0393106702, 1.0414253065, 1.0616810575,
        1.0527424511, 1.0355030998, 1.0464543319, 1.0460720189, 1.0587571336,
        1.0558251125, 1.0379627495, 1.0465498883, 1.0520297849, 1.0491266218,
        1.0444934987, 1.0470275396, 1.0372399286, 1.0523624225, 1.0590876583,
        1.0381554168, 1.0473617659, 1.0542137678, 1.0514117501, 1.0558251125,
        1.0473617659, 1.0525049493, 1.0417613302, 1.0495554623, 1.0547353546,
        1.0536919228, 1.0610686443, 1.0365166036, 1.0608801383, 1.0554461939,
        1.0662863910, 1.0755773647, 1.0459286155, 1.0561565547, 1.0581902797,
        1.0475049731, 1.0431043423, 1.0364201221, 1.0546405398, 1.0463109809,
        1.0516970421, 1.0497936314, 1.0659580985, 1.0428166996, 1.0617752436,
        1.0466932066, 1.0463587667, 1.0386850674, 1.0598427562, 1.0679263400,
        1.0557303956, 1.0479344772, 1.0545931292, 1.0583320216, 1.0454026348,
        1.0700778792, 1.0546405398, 1.0624813729, 1.0597012163, 1.0483638054,
        1.0462631928, 1.0415693299, 1.0455461103, 1.0547353546, 1.0428166996,
        1.0537393739, 1.0591820750, 1.0340051595, 1.0496507364, 1.0608330066,
        1.0692832495, 1.0465498883, 1.0543086209, 1.0512215125, 1.0533597051,
        1.0546879482, 1.0508409339, 1.0476481607, 1.0407048905, 1.0403204647,
        1.0544983016, 1.0562985696, 1.0497936314, 1.0504602174, 1.0584265057,
        1.0551619156, 1.0542137678, 1.0581902797, 1.0799382701, 1.0566298633,
        1.0621048290, 1.0602201035, 1.0500317464, 1.0506030023, 1.0542137678,
        1.0561092122, 1.0587571336, 1.0469797843, 1.0506981814, 1.0760421309,
        1.0525524539, 1.0555883043, 1.0517445832, 1.0541189061, 1.0622931177,
        1.0549249585, 1.0402243360, 1.0568191274, 1.0499365069, 1.0363236317,
        1.0772031690, 1.0568664381, 1.0652542737, 1.0480298987, 1.0611628847,
        1.0531223425, 1.0501269773, 1.0480776062, 1.0528849264, 1.0518871937,
        1.0639862160, 1.0562512334, 1.0554461939, 1.0512690751, 1.0310997382,
        1.0531698193, 1.0382998938, 1.0582847764, 1.0455939311, 1.0503174132,
        1.0339084437, 1.0485068758, 1.0418093247, 1.0438710020, 1.0495554623,
        1.0494125349, 1.0455939311, 1.0302265141, 1.0502221995, 1.0459286155,
        1.0387332041, 1.0357445002, 1.0418573170, 1.0365648410, 1.0437752004,
        1.0393106702, 1.0437272962, 1.0398878157, 1.0447806798, 1.0375291174,
        1.0480776062, 1.0301294435, 1.0402724014, 1.0532647664, 1.0521248349,
        1.0426248937, 1.0429125893, 1.0261903678, 1.0375291174, 1.0461676102,
        1.0491266218, 1.0319722235, 1.0338600823, 1.0487929579, 1.0332795701,
        1.0408970502, 1.0419053071, 1.0415213244, 1.0403204647, 1.0436314814,
        1.0495554623, 1.0589932332, 1.0537393739, 1.0390219773, 1.0356479468,
        1.0816961992, 1.0432002056, 1.0479821891, 1.0481730146, 1.0282833610,
        1.0499365069, 1.0352616431, 1.0276023894, 1.0457851924, 1.0390700984,
        1.0453548052, 1.0462154026, 1.0515544058, 1.0446370991, 1.0567244996,
        1.0422891485, 1.0302750460, 1.0344402688, 1.0384443506, 1.0423850867,
        1.0539765976, 1.0485068758, 1.0369024397, 1.0328439718, 1.0404646410,
        1.0543086209, 1.0484591878, 1.0412332443, 1.0413292799, 1.0493172392,
        1.0422411760, 1.0470275396, 1.0554461939, 1.1023913403, 1.0307602390,
        1.0458330023, 1.0385406441, 1.0456417497, 1.0454026348, 1.0381554168,
        1.0485068758, 1.0387813386, 1.0398878157, 1.0387813386, 1.0603615741,
        1.0537868229, 1.0507457677, 1.0408970502, 1.0293525492, 1.0506981814,
        1.0428646456, 1.0451156246, 1.0443498786, 1.0372881323, 1.0395030877,
        1.0772031690, 1.0568664381, 1.0652542737, 1.0480298987, 1.0611628847,
        1.0531223425, 1.0501269773, 1.0480776062, 1.0528849264, 1.0518871937,
        1.0639862160, 1.0562512334, 1.0554461939, 1.0512690751, 1.0310997382,
        1.0531698193, 1.0382998938, 1.0582847764, 1.0455939311, 1.0503174132,
        1.0339084437, 1.0485068758, 1.0418093247, 1.0438710020, 1.0495554623,
        1.0494125349, 1.0455939311, 1.0302265141, 1.0502221995, 1.0459286155,
        1.0387332041, 1.0357445002, 1.0418573170, 1.0365648410, 1.0437752004,
        1.0393106702, 1.0437272962, 1.0398878157, 1.0447806798, 1.0375291174,
        1.0480776062, 1.0301294435, 1.0402724014, 1.0532647664, 1.0521248349,
        1.0426248937, 1.0429125893, 1.0261903678, 1.0375291174, 1.0461676102,
        1.0491266218, 1.0319722235, 1.0338600823, 1.0487929579, 1.0332795701,
        1.0408970502, 1.0419053071, 1.0415213244, 1.0403204647, 1.0436314814,
        1.0495554623, 1.0589932332, 1.0537393739, 1.0390219773, 1.0356479468,
        1.0816961992, 1.0432002056, 1.0479821891, 1.0481730146, 1.0282833610,
        1.0499365069, 1.0352616431, 1.0276023894, 1.0457851924, 1.0390700984,
        1.0453548052, 1.0462154026, 1.0515544058, 1.0446370991, 1.0567244996,
        1.0422891485, 1.0302750460, 1.0344402688, 1.0384443506, 1.0423850867,
        1.0539765976, 1.0485068758, 1.0369024397, 1.0328439718, 1.0404646410,
        1.0543086209, 1.0484591878, 1.0412332443, 1.0413292799, 1.0493172392,
        1.0422411760, 1.0470275396, 1.0554461939, 1.1023913403, 1.0307602390,
        1.0458330023, 1.0385406441, 1.0456417497, 1.0454026348, 1.0381554168,
        1.0485068758, 1.0387813386, 1.0398878157, 1.0387813386, 1.0603615741,
        1.0537868229, 1.0507457677, 1.0408970502, 1.0293525492, 1.0506981814,
        1.0428646456, 1.0451156246, 1.0443498786, 1.0372881323, 1.0395030877,
    ];

    #[rustfmt::skip]
    const SIGMAS2: [f64; 480] = [
        0.9827851685, 0.9815124495, 0.9634140900, 0.9782467435, 0.9789620476,
        0.9697250642, 0.9655395930, 0.9717853147, 0.9822253753, 0.9771728057,
        0.9671952768, 0.9598784894, 0.9736358127, 0.9783489619, 0.9716309474,
        0.9551788959, 0.9659537800, 0.9612318723, 0.9692093173, 0.9676604274,
        0.9587839780, 0.9772751364, 0.9798809570, 0.9733790094, 0.9722996952,
        0.9703435987, 0.9809519300, 0.9742518643, 0.9789620476, 0.9923037243,
        0.9571137443, 0.9574793568, 0.9714765556, 0.9656949340, 0.9772239724,
        0.9681769994, 0.9791152590, 0.9734303755, 0.9528204071, 0.9711676984,
        0.9727624051, 0.9686416786, 0.9707557360, 0.9510345589, 0.9562252460,
        0.9669884713, 0.9749188143, 0.9683835510, 0.9541313984, 0.9694672251,
        0.9962262189, 0.9772751364, 0.9675570799, 0.9898821560, 0.9629988107,
        0.9627391702, 0.9852749290, 0.9809009578, 0.9701374643, 0.9721454096,
        0.9717338616, 0.9727624051, 0.9774286122, 0.9782978540, 0.9782467435,
        0.9896295685, 0.9605033646, 0.9663677895, 0.9714250862, 0.9755340557,
        0.9637773126, 0.9786044609, 0.9601389035, 0.9825816443, 0.9608156498,
        0.9594616799, 0.9691061350, 0.9681769994, 0.9690545398, 0.9721968408,
        0.9800850420, 0.9677120970, 0.9796768294, 0.9712191815, 0.9741492161,
        0.9746110479, 0.9694672251, 0.9718367650, 0.9636735487, 0.9812577076,
        0.9787066419, 0.9595137910, 0.9820726491, 0.9688997374, 0.9716824059,
        0.9745084375, 0.9660055409, 0.9625314074, 1.0073066440, 0.9700343807,
        0.9795237057, 0.9731735181, 0.9668850520, 0.9767121852, 0.9640885373,
        0.9724025386, 0.9694672251, 0.9744058164, 0.9658502499, 0.9664712641,
        0.9649697956, 0.9681769994, 0.9725053712, 0.9769681123, 0.9659020163,
        0.9819708184, 0.9677637639, 0.9671952768, 0.9765073952, 0.9701890020,
        0.9499824918, 0.9614399159, 0.9500877509, 0.9550741979, 0.9392373321,
        0.9410986982, 0.9488238707, 0.9679704037, 0.9589925527, 0.9519804253,
        0.9620637769, 0.9592010820, 0.9463967159, 0.9572182193, 0.9349688669,
        0.9509819832, 0.9520329458, 0.9494033580, 0.9419483834, 0.9434864876,
        0.9569570104, 0.9424259947, 0.9422668179, 0.9486657673, 0.9388113589,
        0.9433274897, 0.9375856078, 0.9378522117, 0.9457625208, 0.9365718209,
        0.9408330142, 0.9344874445, 0.9405140944, 0.9520854635, 0.9469248875,
        0.9570092579, 0.9432744844, 0.9379055234, 0.9500351228, 0.9395566853,
        0.9241032598, 0.9306808264, 0.9423198798, 0.9426912298, 0.9405672552,
        0.9391840960, 0.9468720836, 0.9399823209, 0.9484022030, 0.9424259947,
        0.9397163212, 0.9382786204, 0.9468192767, 0.9542361999, 0.9380654403,
        0.9253468688, 0.9377988970, 0.9546552908, 0.9392373321, 0.9454981472,
        0.9455510279, 0.9332024375, 0.9388646161, 0.9476110704, 0.9394502463,
        0.9722996952, 0.9361446349, 0.9360378079, 0.9455510279, 0.9329345049,
        0.9364650427, 0.9444928530, 0.9407267196, 0.9401950664, 0.9476110704,
        0.9447575078, 0.9545505354, 0.9380654403, 0.9537120773, 0.9506664679,
        0.9414705299, 0.9424790477, 0.9620637769, 0.9471888628, 0.9400886996,
        0.9471360736, 0.9426381888, 0.9459211094, 0.9432214762, 0.9461325196,
        0.9516127001, 0.9621677144, 0.9433804920, 0.9422668179, 0.9473999901,
        0.9373722701, 0.9518228461, 0.9443869703, 0.9950209444, 0.9475583047,
        0.9390776148, 0.9535547843, 0.9387580986, 0.9500351228, 0.9481912988,
        0.9323447823, 0.9505086709, 0.9483494814, 0.9516652409, 0.9415767405,
        0.9430624336, 0.9286370714, 0.9414174200, 0.9531876664, 0.9281523651,
        0.9509294046, 0.9462910463, 0.9442281240, 0.9312179112, 0.9468192767,
        0.9827851685, 0.9815124495, 0.9634140900, 0.9782467435, 0.9789620476,
        0.9697250642, 0.9655395930, 0.9717853147, 0.9822253753, 0.9771728057,
        0.9671952768, 0.9598784894, 0.9736358127, 0.9783489619, 0.9716309474,
        0.9551788959, 0.9659537800, 0.9612318723, 0.9692093173, 0.9676604274,
        0.9587839780, 0.9772751364, 0.9798809570, 0.9733790094, 0.9722996952,
        0.9703435987, 0.9809519300, 0.9742518643, 0.9789620476, 0.9923037243,
        0.9571137443, 0.9574793568, 0.9714765556, 0.9656949340, 0.9772239724,
        0.9681769994, 0.9791152590, 0.9734303755, 0.9528204071, 0.9711676984,
        0.9727624051, 0.9686416786, 0.9707557360, 0.9510345589, 0.9562252460,
        0.9669884713, 0.9749188143, 0.9683835510, 0.9541313984, 0.9694672251,
        0.9962262189, 0.9772751364, 0.9675570799, 0.9898821560, 0.9629988107,
        0.9627391702, 0.9852749290, 0.9809009578, 0.9701374643, 0.9721454096,
        0.9717338616, 0.9727624051, 0.9774286122, 0.9782978540, 0.9782467435,
        0.9896295685, 0.9605033646, 0.9663677895, 0.9714250862, 0.9755340557,
        0.9637773126, 0.9786044609, 0.9601389035, 0.9825816443, 0.9608156498,
        0.9594616799, 0.9691061350, 0.9681769994, 0.9690545398, 0.9721968408,
        0.9800850420, 0.9677120970, 0.9796768294, 0.9712191815, 0.9741492161,
        0.9746110479, 0.9694672251, 0.9718367650, 0.9636735487, 0.9812577076,
        0.9787066419, 0.9595137910, 0.9820726491, 0.9688997374, 0.9716824059,
        0.9745084375, 0.9660055409, 0.9625314074, 1.0073066440, 0.9700343807,
        0.9795237057, 0.9731735181, 0.9668850520, 0.9767121852, 0.9640885373,
        0.9724025386, 0.9694672251, 0.9744058164, 0.9658502499, 0.9664712641,
        0.9649697956, 0.9681769994, 0.9725053712, 0.9769681123, 0.9659020163,
        0.9819708184, 0.9677637639, 0.9671952768, 0.9765073952, 0.9701890020,
        0.9499824918, 0.9614399159, 0.9500877509, 0.9550741979, 0.9392373321,
        0.9410986982, 0.9488238707, 0.9679704037, 0.9589925527, 0.9519804253,
        0.9620637769, 0.9592010820, 0.9463967159, 0.9572182193, 0.9349688669,
        0.9509819832, 0.9520329458, 0.9494033580, 0.9419483834, 0.9434864876,
        0.9569570104, 0.9424259947, 0.9422668179, 0.9486657673, 0.9388113589,
        0.9433274897, 0.9375856078, 0.9378522117, 0.9457625208, 0.9365718209,
        0.9408330142, 0.9344874445, 0.9405140944, 0.9520854635, 0.9469248875,
        0.9570092579, 0.9432744844, 0.9379055234, 0.9500351228, 0.9395566853,
        0.9241032598, 0.9306808264, 0.9423198798, 0.9426912298, 0.9405672552,
        0.9391840960, 0.9468720836, 0.9399823209, 0.9484022030, 0.9424259947,
        0.9397163212, 0.9382786204, 0.9468192767, 0.9542361999, 0.9380654403,
        0.9253468688, 0.9377988970, 0.9546552908, 0.9392373321, 0.9454981472,
        0.9455510279, 0.9332024375, 0.9388646161, 0.9476110704, 0.9394502463,
        0.9722996952, 0.9361446349, 0.9360378079, 0.9455510279, 0.9329345049,
        0.9364650427, 0.9444928530, 0.9407267196, 0.9401950664, 0.9476110704,
        0.9447575078, 0.9545505354, 0.9380654403, 0.9537120773, 0.9506664679,
        0.9414705299, 0.9424790477, 0.9620637769, 0.9471888628, 0.9400886996,
        0.9471360736, 0.9426381888, 0.9459211094, 0.9432214762, 0.9461325196,
        0.9516127001, 0.9621677144, 0.9433804920, 0.9422668179, 0.9473999901,
        0.9373722701, 0.9518228461, 0.9443869703, 0.9950209444, 0.9475583047,
        0.9390776148, 0.9535547843, 0.9387580986, 0.9500351228, 0.9481912988,
        0.9323447823, 0.9505086709, 0.9483494814, 0.9516652409, 0.9415767405,
        0.9430624336, 0.9286370714, 0.9414174200, 0.9531876664, 0.9281523651,
        0.9509294046, 0.9462910463, 0.9442281240, 0.9312179112, 0.9468192767,
    ];

    #[rustfmt::skip]
    const K_HATS: [f64; 480] = [
       -0.0104000000,  0.0037000000,  0.0026000000,  0.0202000000, -0.0048000000,
        0.0074000000, -0.0115000000,  0.0043000000,  0.0035000000,  0.0121000000,
       -0.0008000000, -0.0065000000,  0.0047000000, -0.0238000000, -0.0265000000,
       -0.0107000000,  0.0000000000,  0.0053000000, -0.0137000000, -0.0111000000,
        0.0007000000,  0.0145000000, -0.0010000000,  0.0009000000, -0.0114000000,
       -0.0006000000, -0.0091000000, -0.0160000000, -0.0056000000, -0.0032000000,
       -0.0042000000, -0.0185000000,  0.0072000000, -0.0024000000,  0.0011000000,
        0.0015000000, -0.0069000000,  0.0026000000,  0.0009000000, -0.0107000000,
        0.0011000000, -0.0076000000, -0.0155000000,  0.0061000000,  0.0019000000,
       -0.0092000000,  0.0085000000,  0.0150000000, -0.0070000000,  0.0157000000,
        0.0154000000,  0.0067000000,  0.0199000000,  0.0025000000, -0.0055000000,
       -0.0011000000, -0.0022000000,  0.0049000000, -0.0053000000, -0.0058000000,
       -0.0421000000, -0.0014000000,  0.0223000000, -0.0078000000,  0.0046000000,
        0.0528000000, -0.0059000000, -0.0062000000, -0.0052000000,  0.0008000000,
       -0.0034000000, -0.0132000000, -0.0005000000,  0.0009000000, -0.0130000000,
        0.0015000000,  0.0119000000, -0.0115000000, -0.0012000000,  0.0046000000,
        0.0058000000, -0.0008000000, -0.0095000000,  0.0207000000, -0.0009000000,
       -0.0219000000,  0.0338000000,  0.0179000000, -0.0139000000, -0.0113000000,
       -0.0015000000,  0.0103000000,  0.0108000000,  0.0019000000,  0.0102000000,
        0.0046000000, -0.0029000000, -0.0192000000,  0.0438000000,  0.0012000000,
       -0.0169000000, -0.0025000000, -0.0081000000, -0.0102000000, -0.0091000000,
        0.0138000000, -0.0005000000,  0.0085000000, -0.0027000000,  0.0159000000,
       -0.0148000000,  0.0163000000,  0.0186000000,  0.0100000000, -0.0083000000,
       -0.0114000000,  0.0108000000, -0.0051000000, -0.0067000000,  0.0056000000,
       -0.0113000000,  0.0097000000,  0.0003000000,  0.0009000000, -0.0152000000,
       -0.0131000000, -0.0082000000, -0.0037000000, -0.0131000000, -0.0152000000,
        0.0000000000,  0.0019000000,  0.0104000000, -0.0005000000,  0.0247000000,
        0.0176000000,  0.0066000000,  0.0048000000,  0.0070000000,  0.0160000000,
        0.0105000000, -0.0140000000, -0.0217000000, -0.0067000000, -0.0136000000,
        0.0099000000, -0.0032000000,  0.0220000000, -0.0038000000,  0.0002000000,
       -0.0099000000, -0.0005000000,  0.0115000000,  0.0190000000, -0.0099000000,
        0.0097000000, -0.0141000000,  0.0059000000,  0.0117000000,  0.0023000000,
       -0.0113000000,  0.0235000000, -0.0066000000,  0.0004000000, -0.0011000000,
        0.0020000000, -0.0067000000, -0.0164000000, -0.0070000000, -0.0058000000,
       -0.0025000000,  0.0015000000,  0.0068000000,  0.0208000000,  0.0130000000,
        0.0020000000, -0.0047000000,  0.0046000000,  0.0273000000,  0.0010000000,
       -0.0055000000, -0.0038000000, -0.0007000000, -0.0100000000, -0.0060000000,
       -0.0444000000,  0.0157000000, -0.0083000000,  0.0116000000,  0.0028000000,
        0.0207000000, -0.0059000000,  0.0173000000, -0.0075000000, -0.0057000000,
        0.0071000000, -0.0016000000,  0.0218000000, -0.0125000000,  0.0054000000,
       -0.0085000000, -0.0006000000,  0.0030000000, -0.0024000000,  0.0008000000,
        0.0016000000, -0.0041000000, -0.0007000000,  0.0085000000, -0.0121000000,
       -0.0114000000, -0.0054000000, -0.0148000000, -0.0011000000,  0.0112000000,
       -0.0005000000, -0.0085000000, -0.0109000000, -0.0280000000,  0.0110000000,
       -0.0061000000, -0.0107000000,  0.0076000000,  0.0028000000, -0.0144000000,
        0.0044000000,  0.0156000000,  0.0222000000,  0.0136000000, -0.0216000000,
       -0.0013000000,  0.0088000000, -0.0024000000, -0.0043000000,  0.0162000000,
        0.0189000000,  0.0101000000, -0.0199000000,  0.0071000000, -0.0175000000,
        0.0039000000, -0.0119000000, -0.0045000000,  0.0132000000,  0.0091000000,
       -0.0005000000, -0.0089000000, -0.0050000000, -0.0063000000, -0.0035000000,
        0.0059000000,  0.0111000000, -0.0037000000,  0.0215000000,  0.0035000000,
        0.0036000000, -0.0037000000, -0.0016000000,  0.0090000000,  0.0047000000,
       -0.0055000000, -0.0289000000, -0.0131000000,  0.0057000000, -0.0443000000,
        0.0063000000, -0.0087000000, -0.0050000000, -0.0036000000, -0.0144000000,
        0.0012000000, -0.0223000000,  0.0020000000,  0.0048000000, -0.0148000000,
       -0.0236000000, -0.0015000000,  0.0089000000,  0.0076000000, -0.0027000000,
        0.0009000000,  0.0015000000,  0.0050000000,  0.0067000000, -0.0034000000,
       -0.0203000000,  0.0035000000,  0.0108000000, -0.0029000000, -0.0014000000,
        0.0266000000, -0.0056000000,  0.0032000000, -0.0049000000,  0.0141000000,
       -0.0045000000,  0.0090000000,  0.0022000000,  0.0000000000, -0.0059000000,
        0.0021000000,  0.0152000000, -0.0009000000,  0.0056000000, -0.0144000000,
       -0.0263000000,  0.0026000000, -0.0072000000,  0.0036000000, -0.0142000000,
       -0.0167000000,  0.0004000000, -0.0093000000, -0.0122000000,  0.0065000000,
       -0.0198000000,  0.0163000000, -0.0131000000,  0.0136000000,  0.0218000000,
       -0.0187000000,  0.0126000000,  0.0130000000,  0.0048000000,  0.0079000000,
        0.0018000000,  0.0122000000, -0.0037000000, -0.0054000000, -0.0051000000,
        0.0285000000, -0.0288000000, -0.0134000000, -0.0146000000,  0.0028000000,
        0.0092000000,  0.0104000000,  0.0024000000,  0.0429000000,  0.0119000000,
       -0.0102000000, -0.0004000000,  0.0072000000,  0.0101000000, -0.0044000000,
       -0.0081000000,  0.0112000000,  0.0075000000, -0.0009000000, -0.0134000000,
        0.0106000000,  0.0049000000, -0.0020000000,  0.0099000000, -0.0049000000,
        0.0058000000,  0.0055000000,  0.0021000000, -0.0047000000,  0.0039000000,
       -0.0045000000,  0.0050000000, -0.0034000000, -0.0075000000,  0.0111000000,
        0.0010000000,  0.0004000000, -0.0015000000, -0.0076000000, -0.0014000000,
       -0.0066000000,  0.0153000000, -0.0040000000, -0.0210000000, -0.0025000000,
        0.0118000000,  0.0035000000,  0.0075000000, -0.0103000000,  0.0221000000,
       -0.0054000000,  0.0043000000, -0.0089000000,  0.0117000000, -0.0039000000,
        0.0117000000, -0.0108000000,  0.0029000000,  0.0029000000, -0.0180000000,
        0.0130000000, -0.0101000000,  0.0254000000,  0.0035000000,  0.0131000000,
       -0.0114000000, -0.0080000000,  0.0086000000, -0.0004000000,  0.0147000000,
       -0.0247000000,  0.0005000000, -0.0082000000, -0.0006000000,  0.0009000000,
        0.0002000000,  0.0002000000,  0.0179000000,  0.0082000000,  0.0059000000,
        0.0076000000,  0.0025000000,  0.0014000000, -0.0022000000,  0.0043000000,
       -0.0067000000,  0.0157000000,  0.0168000000, -0.0136000000,  0.0058000000,
       -0.0113000000, -0.0114000000,  0.0226000000, -0.0061000000,  0.0124000000,
        0.0413000000, -0.0130000000,  0.0019000000, -0.0091000000,  0.0258000000,
        0.0187000000,  0.0080000000,  0.0053000000, -0.0005000000,  0.0035000000,
       -0.0053000000,  0.0017000000,  0.0078000000,  0.0126000000,  0.0020000000,
        0.0003000000, -0.0101000000,  0.0168000000, -0.0154000000,  0.0141000000,
        0.0095000000, -0.0057000000,  0.0155000000,  0.0046000000, -0.0061000000,
       -0.0147000000,  0.0063000000, -0.0200000000,  0.0002000000, -0.0067000000,
        0.0126000000,  0.0090000000,  0.0211000000,  0.0523000000,  0.0192000000,
        0.0051000000,  0.0124000000, -0.0128000000, -0.0013000000, -0.0080000000,
        0.0086000000, -0.0193000000,  0.0084000000,  0.0067000000,  0.0231000000,
        0.0079000000,  0.0041000000,  0.0140000000, -0.0025000000, -0.0170000000,
       -0.0030000000,  0.0034000000, -0.0111000000,  0.0065000000,  0.0102000000,
    ];

    #[rustfmt::skip]
    const KAPPAS: [f64; 480] = [
       -0.01040000012,  0.00370000004,  0.00260000006,  0.02020000028, -0.00480000006,
        0.00740000013, -0.01150000024,  0.00430000007,  0.00350000004,  0.01210000016,
       -0.00080000001, -0.00650000017,  0.00470000007, -0.02380000031, -0.02650000044,
       -0.01070000032,  0.00000000000,  0.00530000013, -0.01370000025, -0.01110000021,
        0.00070000001,  0.01450000020, -0.00100000001,  0.00090000001, -0.01140000018,
       -0.00060000001, -0.00910000011, -0.01600000024, -0.00560000007, -0.00320000002,
       -0.00420000012, -0.01850000052,  0.00720000012, -0.00240000005,  0.00110000001,
        0.00150000002, -0.00690000008,  0.00260000004,  0.00090000003, -0.01070000018,
        0.00110000001, -0.00760000014, -0.01550000027,  0.00610000021,  0.00190000005,
       -0.00920000018,  0.00850000012,  0.01500000029, -0.00700000022,  0.01570000028,
        0.01540000010,  0.00670000009,  0.01990000038,  0.00250000002, -0.00550000012,
       -0.00110000002, -0.00220000002,  0.00490000006, -0.00530000009, -0.00580000009,
       -0.04210000070, -0.00140000002,  0.02230000030, -0.00780000010,  0.00460000006,
        0.05280000044, -0.00590000014, -0.00620000012, -0.00520000008,  0.00080000001,
       -0.00340000007, -0.01320000017, -0.00050000001,  0.00090000001, -0.01300000032,
        0.00150000003,  0.01190000021, -0.01150000022, -0.00120000002,  0.00460000007,
        0.00580000007, -0.00080000001, -0.00950000011,  0.02070000035, -0.00090000001,
       -0.02190000032,  0.03380000061,  0.01790000030, -0.01390000031, -0.01130000013,
       -0.00150000001,  0.01030000027,  0.01080000012,  0.00190000003,  0.01020000016,
        0.00460000006, -0.00290000006, -0.01920000044,  0.04380000018,  0.00120000002,
       -0.01690000020, -0.00250000003, -0.00810000016, -0.01020000014, -0.00910000020,
        0.01380000022, -0.00050000000,  0.00850000013, -0.00270000005,  0.01590000031,
       -0.01480000031,  0.01630000031,  0.01860000030,  0.01000000013, -0.00830000017,
       -0.01140000013,  0.01080000021, -0.00510000010, -0.00670000009,  0.00560000010,
       -0.01130000042,  0.00970000024,  0.00030000001,  0.00090000002, -0.01520000084,
       -0.01310000067, -0.00820000032, -0.00370000007, -0.01310000035, -0.01520000052,
        0.00000000000,  0.00190000005,  0.01040000044, -0.00050000001,  0.02470000159,
        0.01760000063,  0.00660000023,  0.00480000018,  0.00700000035,  0.01600000076,
        0.01050000030, -0.01400000069, -0.02170000108, -0.00670000026, -0.01360000076,
        0.00990000047, -0.00320000018,  0.02200000128, -0.00380000016,  0.00020000001,
       -0.00990000051, -0.00050000003,  0.01150000061,  0.01900000065, -0.00990000041,
        0.00970000027, -0.01410000067,  0.00590000034,  0.01170000043,  0.00230000012,
       -0.01130000108,  0.02350000177, -0.00660000032,  0.00040000001, -0.00110000005,
        0.00200000011, -0.00670000028, -0.01640000089, -0.00700000027, -0.00580000028,
       -0.00250000013,  0.00150000008,  0.00680000028,  0.02080000067,  0.01300000074,
        0.00200000018, -0.00470000027,  0.00460000014,  0.02730000150,  0.00100000004,
       -0.00550000024, -0.00380000026, -0.00070000003, -0.01000000040, -0.00600000033,
       -0.04440000071,  0.01570000097, -0.00830000051,  0.01160000050,  0.00280000019,
        0.02070000127, -0.00590000027,  0.01730000091, -0.00750000040, -0.00570000023,
        0.00710000032, -0.00160000005,  0.02180000127, -0.01250000040,  0.00540000019,
       -0.00850000043, -0.00060000003,  0.00300000007, -0.00240000009,  0.00080000004,
        0.00160000006, -0.00410000020, -0.00070000003,  0.00850000041, -0.01210000052,
       -0.01140000040, -0.00540000012, -0.01480000070, -0.00110000005,  0.01120000045,
       -0.00050000003, -0.00850000029, -0.01090000049, -0.02800000018,  0.01100000044,
       -0.00610000033, -0.01070000034,  0.00760000042,  0.00280000010, -0.01440000057,
        0.00440000031,  0.01560000056,  0.02220000088,  0.01360000047, -0.02160000109,
       -0.00130000006,  0.00880000071, -0.00240000012, -0.00430000014,  0.01620000134,
        0.01890000068,  0.01010000043, -0.01990000091,  0.00710000052, -0.01750000074,
        0.00390000004, -0.01190000013, -0.00450000010,  0.01320000017,  0.00910000011,
       -0.00050000000, -0.00890000018, -0.00500000008, -0.00630000007, -0.00350000004,
        0.00590000011,  0.01110000028, -0.00370000005,  0.02150000028,  0.00350000006,
        0.00360000011, -0.00370000007, -0.00160000003,  0.00900000016,  0.00470000009,
       -0.00550000015, -0.02890000039, -0.01310000016,  0.00570000009, -0.04430000072,
        0.00630000011, -0.00870000010, -0.00500000008, -0.00360000004, -0.01440000011,
        0.00120000003, -0.02230000064,  0.00200000003,  0.00480000010, -0.01480000020,
       -0.02360000046, -0.00150000002,  0.00890000014,  0.00760000025, -0.00270000004,
        0.00090000001,  0.00150000003,  0.00500000008,  0.00670000024, -0.00340000010,
       -0.02030000041,  0.00350000005,  0.01080000021, -0.00290000009, -0.00140000002,
        0.02660000018, -0.00560000008,  0.00320000006, -0.00490000004,  0.01410000033,
       -0.00450000010,  0.00900000009,  0.00220000002,  0.00000000000, -0.00590000009,
        0.00210000003,  0.01520000024, -0.00090000001,  0.00560000007, -0.01440000019,
       -0.02630000021,  0.00260000006, -0.00720000014,  0.00360000006, -0.01420000022,
       -0.01670000037,  0.00040000000, -0.00930000024, -0.01220000013,  0.00650000016,
       -0.01980000052,  0.01630000030, -0.01310000025,  0.01360000025,  0.02180000035,
       -0.01870000023,  0.01260000025,  0.01300000017,  0.00480000008,  0.00790000012,
        0.00180000002,  0.01220000022, -0.00370000006, -0.00540000012, -0.00510000006,
        0.02850000036, -0.02880000076, -0.01340000015, -0.01460000027,  0.00280000004,
        0.00920000014,  0.01040000021,  0.00240000005,  0.04290000017,  0.01190000022,
       -0.01020000013, -0.00040000000,  0.00720000014,  0.01010000014, -0.00440000010,
       -0.00810000013,  0.01120000020,  0.00750000011, -0.00090000001, -0.01340000027,
        0.01060000022,  0.00490000009, -0.00200000003,  0.00990000014, -0.00490000010,
        0.00580000006,  0.00550000010,  0.00210000004, -0.00470000006,  0.00390000007,
       -0.00450000016,  0.00500000012, -0.00340000012, -0.00750000023,  0.01110000061,
        0.00100000005,  0.00040000001, -0.00150000002, -0.00760000020, -0.00140000004,
       -0.00660000015,  0.01530000040, -0.00400000017, -0.02100000060, -0.00250000016,
        0.01180000042,  0.00350000012,  0.00750000028, -0.01030000051,  0.02210000104,
       -0.00540000016,  0.00430000021, -0.00890000044,  0.01170000046, -0.00390000022,
        0.01170000055, -0.01080000063,  0.00290000017,  0.00290000012, -0.01800000110,
        0.01300000068, -0.01010000067,  0.02540000135,  0.00350000012,  0.01310000055,
       -0.01140000033, -0.00800000038,  0.00860000050, -0.00040000001,  0.01470000081,
       -0.02470000237,  0.00050000003, -0.00820000040, -0.00060000002,  0.00090000004,
        0.00020000001,  0.00020000000,  0.01790000098,  0.00820000032,  0.00590000029,
        0.00760000041,  0.00250000014,  0.00140000006, -0.00220000007,  0.00430000025,
       -0.00670000061,  0.01570000092,  0.01680000053, -0.01360000075,  0.00580000025,
       -0.01130000049, -0.01140000078,  0.02260000126, -0.00610000025,  0.01240000069,
        0.04130000065, -0.01300000081,  0.00190000011, -0.00910000040,  0.02580000182,
        0.01870000114,  0.00800000037,  0.00530000028, -0.00050000002,  0.00350000014,
       -0.00530000024,  0.00170000005,  0.00780000045,  0.01260000041,  0.00200000007,
        0.00030000001, -0.01010000050,  0.01680000041, -0.01540000064,  0.01410000076,
        0.00950000039, -0.00570000027,  0.01550000068,  0.00460000022, -0.00610000026,
       -0.01470000051,  0.00630000015, -0.02000000095,  0.00020000001, -0.00670000027,
        0.01260000075,  0.00900000031,  0.02110000096,  0.05230000034,  0.01920000080,
        0.00510000028,  0.01240000041, -0.01280000072, -0.00130000004, -0.00800000032,
        0.00860000061, -0.01930000071,  0.00840000033,  0.00670000023,  0.02310000117,
        0.00790000038,  0.00410000033,  0.01400000072, -0.00250000008, -0.01700000141,
       -0.00300000010,  0.00340000014, -0.01110000051,  0.00650000048,  0.01020000043,
    ];

    #[test]
    // compare values from pyuvdata
    fn test_van_vleck_crosses_int() {
        let result = van_vleck_crosses_int(&K_HATS, &SIGMAS1, &SIGMAS2);
        println!("{:?}", result);
        for (&r, &s) in result.iter().zip(KAPPAS.iter()) {
            assert_approx_eq!(f64, r, s, epsilon = 1e-10);
        }
    }
}
