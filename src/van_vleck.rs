use crate::{
    ndarray::{parallel::prelude::*, prelude::*},
    BirliError, Jones,
};
use errorfunctions::RealErrorFunctions;
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use itertools::{izip, Itertools};
use log::trace;
use marlu::{
    constants::VEL_C,
    hifitime::{Duration, Epoch},
    io::error::BadArrayShape,
    mwalib::{CorrelatorContext, MWAVersion},
    precession::precess_time,
    Complex, LatLngHeight, RADec, VisSelection, XyzGeodetic, UVW,
};
use std::{
    f32::consts::SQRT_2,
    f64::consts::TAU,
    ops::{Index, Range},
};
use thiserror::Error;

#[derive(Error, Debug)]
/// Error for Passband Corrections
pub enum VanVleckCorrection {
    #[error(transparent)]
    /// Error for bad array shape in provided argument
    BadArrayShape(#[from] BadArrayShape),
}

/// Perform Van Vleck corrections given an observation's
/// [`mwalib::CorrelatorContext`](crate::mwalib::CorrelatorContext) and an
/// [`ndarray::Array3`](crate::ndarray::Array3) of [`marlu::Jones`] visibilities (time, frequency, baseline)
///
/// Credit: https://github.com/RadioAstronomySoftwareGroup/pyuvdata/blob/3b412bc5397dd2a87a3bb6239e9c584b52306603/src/pyuvdata/uvdata/mwa_corr_fits.py#L775
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
/// # Errors:
/// - [`VanVleckCorrection::BadArrayShape`] - If the shape of the provided `ant_pairs` argument is incorrect.
pub fn correct_van_vleck(
    corr_ctx: &CorrelatorContext,
    mut jones_array: ArrayViewMut3<Jones<f32>>,
    // baseline_idxs: &[usize],
    ant_pairs: &[(usize, usize)],
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
    let nsamples = (corr_ctx.metafits_context.corr_fine_chan_width_hz)
        * (corr_ctx.metafits_context.corr_int_time_ms as u32 / 1_000);
    dbg!(&nsamples);

    // TODO: double precision, but ensure this is applied everywhere it needs to be applied
    // jones_array.mapv_inplace(|j| j / (2.0 * nsamples as f32));
    let jones_double = jones_array.mapv(|j| Jones::<f64>::from(j) / (2.0 * nsamples as f64));

    // TODO: figure out antenna flags
    // TODO: ensure not mwax

    // mask for ant_pairs which are autocorrelations
    let auto_mask = ant_pairs
        .iter()
        .enumerate()
        .filter_map(|(i, (ant1, ant2))| if ant1 == ant2 { Some(i) } else { None })
        .collect::<Vec<usize>>();
    let n_autos = auto_mask.len();
    dbg!(&n_autos);

    // new mutable copy of jones_array, only autos
    // TODO: pyuvdata does `[timestep][baseleine][channel]` -> `[baseleine][timestep * channel]`
    // but we do `[timestep][channel][baseleine]`.
    let mut jones_autos = jones_double.select(Axis(2), &auto_mask);

    dbg!(&jones_autos.dim());

    // square root auto xx and yys
    // TODO: what about cross-pols?
    // TODO: assuming everything is real?
    // TODO: this is going to be annoying to unpack.
    let sighat = jones_autos
        .iter()
        .flat_map(|j| [j[0].re.sqrt(), j[3].re.sqrt()])
        .collect_vec();

    // TODO: get max value of jones_array
    let max_sighat = sighat.iter().fold(0.0_f64, |acc, &j| acc.max(j));
    dbg!(&max_sighat);

    // correct autos
    let sigma = van_vleck_autos(sighat);

    Ok(())
}

/// Use Newton's method to solve the inverse of sighat_vector.
/// TODO: error handling for when it doesn't converge
pub fn van_vleck_autos(hat: Vec<f64>) -> Vec<f64> {
    let tol = 1e-10;
    hat.par_iter()
        .enumerate()
        .map(|(i, &s)| {
            // cut off small sigmas that will not converge
            if s < 0.5 {
                return s;
            }
            let mut sigma = s;
            let mut delta = sighat(sigma) - s;
            let mut count = 0;
            while delta.abs() > tol {
                sigma -= delta / sighat_prime(sigma);
                delta = sighat(sigma) - s;
                count += 1;
                if count > 100 {
                    log::warn!(
                        "Van Vleck correction did not converge for sigma={} at index {}",
                        s,
                        i
                    );
                    return s;
                }
            }
            sigma
        })
        .collect()
}

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

/// The derivative of the $\hat\sigma$ function .
/// $$\hat\sigma^\prime(\sigma_i) = \sum_{k=0}^6{(2k+1)(k+0.5)e^{-(k+0.5)^2/2\sigma_i^2}\over\sqrt{2\pi}\sigma_i^2}$$
/// substitute $K=(k+0.5)^2$, $S=\sigma_i^2$:
/// $$ = \sum_{k=0}^6{(2K^2e^{-K^2/2S}\over\sqrt{2\pi}S}$$
fn sighat_prime(sigma: f64) -> f64 {
    let s_ = sigma.powi(2);
    let mut sum = 0.0;
    for k in 0..=6 {
        let k_ = (k as f64 + 0.5).powi(2);
        sum += 2.0 * k_ * (-(k_ / (2.0 * s_))).exp() / (SQRT_TAU * s_);
    }
    sum / sighat(sigma)
}

#[cfg(test)]
mod tests {
    use float_cmp::assert_approx_eq;

    use super::*;

    const sighats: [f64; 20] = [
        1.3732557118031588,
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
    const sigmas: [f64; 20] = [
        1.3425715134733938,
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
    const sighats_prime: [f64; 20] = [
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
        let result = sigmas.iter().map(|&s| sighat(s)).collect_vec();
        for (&r, &s) in result.iter().zip(sighats.iter()) {
            assert_approx_eq!(f64, r, s, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_sighat_prime() {
        let result = sigmas.iter().map(|&s| sighat_prime(s)).collect_vec();
        for (&r, &s) in result.iter().zip(sighats_prime.iter()) {
            assert_approx_eq!(f64, r, s, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_van_vleck_autos() {
        let result = van_vleck_autos(sighats.to_vec());
        for (&r, &s) in result.iter().zip(sigmas.iter()) {
            assert_approx_eq!(f64, r, s, epsilon = 1e-6);
        }
    }
}
