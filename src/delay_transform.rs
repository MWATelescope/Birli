//! Delay transform of real-valued spectra.
//!
//! A delay transform is the Fourier transform of a spectrum along its frequency
//! axis. A ripple with period `P` Hz in a spectrum shows up as a peak at a delay
//! of `1 / P` seconds, so the delay spectrum of a tile's bandpass shows cable
//! reflections and other standing waves as peaks at the delay of the
//! reflection.
//!
//! [`delay_transform`] takes one spectrum per row (for example, the
//! time-averaged autocorrelation of each tile), and for each row it:
//!
//! 1. replaces non-finite samples (flagged channels) with zero,
//! 2. applies a Blackman window, to keep the sidelobes of the strong
//!    zero-delay term out of the delays of interest,
//! 3. zero-pads the spectrum by a whole number of times, so that the delay
//!    bins are no wider than [`DelayTransformConfig::target_delay_res_ns`],
//! 4. takes the amplitude of the forward FFT,
//! 5. keeps only the delays between [`DelayTransformConfig::min_delay_ns`] and
//!    [`DelayTransformConfig::max_delay_ns`].
//!
//! The amplitudes are not normalised: they scale with the number of channels
//! and with the window.
//!
//! The spectra are real, so the delay spectrum is symmetric: for a channel
//! width `df`, delays above `1 / (2 df)` are aliases of the delays below it.
//!
//! # Examples
//!
//! A 2 MHz ripple shows up as a peak at 500 ns.
//!
//! ```rust
//! use birli::{delay_transform, marlu::ndarray::{Array1, Array2}, DelayTransformConfig};
//! use std::f64::consts::TAU;
//!
//! // 1024 channels of 40 kHz
//! let freq_hz = Array1::from_iter((0..1024).map(|i| 150e6 + 40e3 * i as f64));
//! // one tile, whose bandpass has a 2 MHz ripple
//! let spectrum = Array2::from_shape_fn((1, 1024), |(_, i)| {
//!     1.0 + 0.1 * (TAU * freq_hz[i] / 2e6).cos()
//! });
//! let config = DelayTransformConfig {
//!     min_delay_ns: 100.0,
//!     max_delay_ns: 3000.0,
//!     target_delay_res_ns: 1.0,
//! };
//!
//! let result = delay_transform(&spectrum, &freq_hz, &config).unwrap();
//!
//! let (peak_idx, _) = result
//!     .delay_spectrum
//!     .row(0)
//!     .iter()
//!     .enumerate()
//!     .max_by(|a, b| a.1.total_cmp(b.1))
//!     .unwrap();
//! assert!((result.delays_ns[peak_idx] - 500.0).abs() < 2.0);
//! ```

use std::{f64::consts::PI, ops::Range};

use log::trace;
use marlu::ndarray::{Array1, Array2};
use rustfft::{num_complex::Complex64, FftPlanner};
use thiserror::Error;

/// If fewer than this many delay bins fall between
/// [`DelayTransformConfig::min_delay_ns`] and
/// [`DelayTransformConfig::max_delay_ns`], the requested range is too narrow
/// to be useful and [`delay_transform`] falls back to a broader range.
const MIN_DELAY_BINS: usize = 6;

/// Errors that can occur in a delay transform
#[derive(Error, Debug)]
pub enum DelayTransformError {
    /// There are no delays at or below `max_delay_ns`, which only happens when
    /// it is negative.
    #[error("no delays in the requested range, max_delay_ns={max_delay_ns}")]
    NoDelaysInRange {
        /// The maximum delay that was asked for \[ns\]
        max_delay_ns: f64,
    },
}

/// Configuration for the delay transform.
#[derive(Debug, Clone, Copy)]
pub struct DelayTransformConfig {
    /// Minimum delay to include in output \[ns\]
    pub min_delay_ns: f64,
    /// Maximum delay to include in output \[ns\]
    pub max_delay_ns: f64,
    /// Target delay resolution \[ns\]. The spectrum is zero-padded until the
    /// delay bins are at most this wide.
    pub target_delay_res_ns: f64,
}

impl Default for DelayTransformConfig {
    fn default() -> Self {
        Self {
            min_delay_ns: 0.0,
            max_delay_ns: 1000.0,
            target_delay_res_ns: 1.0,
        }
    }
}

/// Result of the delay transform.
#[derive(Debug)]
pub struct DelayTransformResult {
    /// Amplitude of the delay spectrum of each input row `(nants, ndelays)`
    pub delay_spectrum: Array2<f64>,
    /// The delay of each bin in `delay_spectrum` \[ns\]
    pub delays_ns: Array1<f64>,
}

/// The dimensions of a delay transform, see [`calculate_delay_channels`].
#[derive(Debug)]
pub struct DelayChannelInfo {
    /// Number of delay bins in the output
    pub n_delay_channels: usize,
    /// Delay of the first bin in the output \[ns\]
    pub start_delay_ns: f64,
    /// Number of frequency channels after zero-padding, which is also the
    /// number of delay bins before the requested range is selected
    pub nfreqs_padded: usize,
    /// Actual delay resolution achieved \[ns\]
    pub actual_delay_res_ns: f64,
    /// Delay of the last bin before the requested range is selected \[ns\]
    pub max_computable_delay_ns: f64,
}

/// The indices of the padded delay bins that end up in the output.
///
/// - If at least [`MIN_DELAY_BINS`] bins are between `config.min_delay_ns` and
///   `config.max_delay_ns`, it is those bins.
/// - Otherwise the range that was asked for is too narrow, and `min_delay_ns`
///   is ignored: it is all the bins up to `config.max_delay_ns`.
///
/// `None` if there are no such bins, which needs a negative `max_delay_ns`.
fn selected_delay_bins(
    nfreqs_padded: usize,
    delay_res_ns: f64,
    config: &DelayTransformConfig,
) -> Option<Range<usize>> {
    let delay_ns = |i: usize| i as f64 * delay_res_ns;
    let in_range = |i: &usize| {
        let delay_ns = delay_ns(*i);
        delay_ns >= config.min_delay_ns && delay_ns <= config.max_delay_ns
    };

    // delays increase with bin index, so the bins in range are contiguous
    let first = (0..nfreqs_padded).find(in_range);
    let last = (0..nfreqs_padded).rev().find(in_range);
    if let (Some(first), Some(last)) = (first, last) {
        if last + 1 - first >= MIN_DELAY_BINS {
            return Some(first..last + 1);
        }
    }

    let num_below_max = (0..nfreqs_padded)
        .take_while(|&i| delay_ns(i) <= config.max_delay_ns)
        .count();
    (num_below_max > 0).then_some(0..num_below_max)
}

/// Calculate the dimensions of the result of [`delay_transform`] for spectra
/// of `nfreqs` channels at frequencies `freq_hz`, without doing the transform.
///
/// The channel width is taken from the first two frequencies, which should be
/// evenly spaced.
///
/// # Panics
///
/// Panics if `freq_hz` has fewer than two channels.
pub fn calculate_delay_channels(
    nfreqs: usize,
    freq_hz: &Array1<f64>,
    config: &DelayTransformConfig,
) -> DelayChannelInfo {
    let freq_res_hz = freq_hz[1] - freq_hz[0];

    // Zero-padding factor to reach the target delay resolution
    let pad_factor =
        ((1e9 / (config.target_delay_res_ns * freq_res_hz * nfreqs as f64)).ceil() as usize).max(1);
    let nfreqs_padded = nfreqs * pad_factor;

    let actual_delay_res_ns = 1e9 / (nfreqs_padded as f64 * freq_res_hz);
    let max_computable_delay_ns = (nfreqs_padded - 1) as f64 * actual_delay_res_ns;

    let (n_delay_channels, start_delay_ns) =
        selected_delay_bins(nfreqs_padded, actual_delay_res_ns, config).map_or((1, 0.0), |bins| {
            (bins.len(), bins.start as f64 * actual_delay_res_ns)
        });

    DelayChannelInfo {
        n_delay_channels,
        start_delay_ns,
        nfreqs_padded,
        actual_delay_res_ns,
        max_computable_delay_ns,
    }
}

/// Blackman window (same as `numpy.blackman`)
fn blackman_window(n: usize) -> Array1<f64> {
    if n == 1 {
        return Array1::ones(1);
    }
    Array1::from_shape_fn(n, |i| {
        let a = 2.0 * PI * i as f64 / (n as f64 - 1.0);
        0.42 - 0.5 * a.cos() + 0.08 * (2.0 * a).cos()
    })
}

/// Perform a delay transform (the amplitude of the FFT along the frequency
/// axis) of each row of `rdx_ant_freq`. See the [module docs](self) for the
/// details.
///
/// - `rdx_ant_freq` - spectra with dimensions `(nants, nfreqs)`. Non-finite
///   samples are treated as zero. A row with no finite samples gives a row of
///   `NaN`.
/// - `freq_hz` - the centre frequency of each channel \[Hz\], which should be
///   evenly spaced.
///
/// The number of delay bins in the result can be calculated ahead of time with
/// [`calculate_delay_channels`].
///
/// # Errors
///
/// Returns [`DelayTransformError::NoDelaysInRange`] if no delays are at or
/// below `config.max_delay_ns`.
///
/// # Panics
///
/// Panics if `freq_hz` has fewer than two channels.
pub fn delay_transform(
    rdx_ant_freq: &Array2<f64>,
    freq_hz: &Array1<f64>,
    config: &DelayTransformConfig,
) -> Result<DelayTransformResult, DelayTransformError> {
    let (nants, nfreqs) = rdx_ant_freq.dim();

    let delay_info = calculate_delay_channels(nfreqs, freq_hz, config);
    let nfreqs_padded = delay_info.nfreqs_padded;
    let delay_res_ns = delay_info.actual_delay_res_ns;
    trace!(
        "delay transform: zero-padding {nfreqs} -> {nfreqs_padded} channels for {delay_res_ns:.2} ns resolution"
    );

    let delay_bins = selected_delay_bins(nfreqs_padded, delay_res_ns, config).ok_or(
        DelayTransformError::NoDelaysInRange {
            max_delay_ns: config.max_delay_ns,
        },
    )?;

    let fft = FftPlanner::new().plan_fft_forward(nfreqs_padded);
    let window = blackman_window(nfreqs);

    let mut delay_spectrum = Array2::from_elem((nants, delay_bins.len()), f64::NAN);
    let mut buffer = vec![Complex64::new(0.0, 0.0); nfreqs_padded];
    for (spectrum, mut delay_spectrum) in rdx_ant_freq
        .outer_iter()
        .zip(delay_spectrum.outer_iter_mut())
    {
        if !spectrum.iter().any(|x| x.is_finite()) {
            continue;
        }
        buffer.fill(Complex64::new(0.0, 0.0));
        for ((buf, &x), &w) in buffer.iter_mut().zip(spectrum.iter()).zip(window.iter()) {
            if x.is_finite() {
                *buf = Complex64::new(x * w, 0.0);
            }
        }
        fft.process(&mut buffer);
        for (out, transformed) in delay_spectrum.iter_mut().zip(&buffer[delay_bins.clone()]) {
            *out = transformed.norm();
        }
    }

    let delays_ns = delay_bins.map(|i| i as f64 * delay_res_ns).collect();

    Ok(DelayTransformResult {
        delay_spectrum,
        delays_ns,
    })
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;
    use std::f64::consts::TAU;

    use super::*;

    /// 128 channels of 10 kHz, starting at 133 MHz.
    fn fine_freqs_hz() -> Array1<f64> {
        Array1::from_iter((0..128).map(|i| 133e6 + i as f64 * 10e3))
    }

    #[test]
    fn test_blackman_window() {
        // numpy.blackman(5)
        let window = blackman_window(5);
        for (&obs, exp) in window.iter().zip([0.0, 0.34, 1.0, 0.34, 0.0]) {
            assert_abs_diff_eq!(obs, exp, epsilon = 1e-12);
        }
        // numpy.blackman(1)
        assert_eq!(blackman_window(1).to_vec(), vec![1.0]);
        assert_eq!(blackman_window(0).len(), 0);
    }

    #[test]
    fn test_calculate_delay_channels() {
        let nfreqs = 64;
        let freq_hz = Array1::linspace(100e6, 200e6, nfreqs);
        let config = DelayTransformConfig {
            min_delay_ns: 6.0,
            max_delay_ns: 60.0,
            target_delay_res_ns: 0.5,
        };

        let info = calculate_delay_channels(nfreqs, &freq_hz, &config);
        // the channels are 100/63 MHz wide, so 64 channels give 9.84 ns
        // resolution, and they need to be padded 20 times to get under 0.5 ns
        assert_eq!(info.nfreqs_padded, 64 * 20);
        assert_abs_diff_eq!(info.actual_delay_res_ns, 0.4921875, epsilon = 1e-9);
        assert_abs_diff_eq!(
            info.max_computable_delay_ns,
            1279.0 * 0.4921875,
            epsilon = 1e-6
        );
        // bins 13 (6.40 ns) to 121 (59.55 ns)
        assert_eq!(info.n_delay_channels, 109);
        assert_abs_diff_eq!(info.start_delay_ns, 13.0 * 0.4921875, epsilon = 1e-9);
    }

    #[test]
    fn test_no_padding_when_resolution_is_already_fine_enough() {
        // 128 channels of 10 kHz give 781.25 ns resolution
        let freq_hz = fine_freqs_hz();
        let config = DelayTransformConfig {
            min_delay_ns: 0.0,
            max_delay_ns: 1e6,
            target_delay_res_ns: 1000.0,
        };
        let info = calculate_delay_channels(128, &freq_hz, &config);
        assert_eq!(info.nfreqs_padded, 128);
        assert_abs_diff_eq!(info.actual_delay_res_ns, 781.25, epsilon = 1e-9);
        assert_eq!(info.n_delay_channels, 128);
    }

    #[test]
    fn test_delay_transform_shape() {
        let nants = 3;
        let nfreqs = 64;
        let rdx_ant_freq = Array2::ones((nants, nfreqs));
        let freq_hz = Array1::linspace(100e6, 200e6, nfreqs);
        let config = DelayTransformConfig::default();

        // calculate_delay_channels predicts the size of the result
        let info = calculate_delay_channels(nfreqs, &freq_hz, &config);
        let result = delay_transform(&rdx_ant_freq, &freq_hz, &config).unwrap();

        assert_eq!(result.delay_spectrum.dim(), (nants, info.n_delay_channels));
        assert_eq!(result.delays_ns.len(), info.n_delay_channels);
        // the default range starts at zero delay, and the bins are evenly spaced
        assert_abs_diff_eq!(result.delays_ns[0], 0.0);
        assert_abs_diff_eq!(info.start_delay_ns, 0.0);
        assert_abs_diff_eq!(
            result.delays_ns[1] - result.delays_ns[0],
            info.actual_delay_res_ns,
            epsilon = 1e-9
        );
        assert!(*result.delays_ns.last().unwrap() <= config.max_delay_ns);
    }

    /// A ripple in the spectrum gives a peak at the delay of the ripple, with
    /// an amplitude that scales with the amplitude of the ripple.
    #[test]
    fn test_delay_transform_finds_ripple() {
        // 1024 channels of 40 kHz give 24 ns resolution before padding
        let freq_hz = Array1::from_iter((0..1024).map(|i| 150e6 + i as f64 * 40e3));
        let ripples = [(0.1, 250e-9), (0.2, 400e-9), (0.4, 400e-9)];
        let rdx_ant_freq = Array2::from_shape_fn((ripples.len(), freq_hz.len()), |(ant, f)| {
            let (amp, delay_s) = ripples[ant];
            1e6 * (1.0 + amp * (TAU * freq_hz[f] * delay_s).cos())
        });
        let config = DelayTransformConfig {
            min_delay_ns: 100.0,
            max_delay_ns: 3000.0,
            target_delay_res_ns: 1.0,
        };

        let result = delay_transform(&rdx_ant_freq, &freq_hz, &config).unwrap();

        assert!(result.delays_ns[0] >= config.min_delay_ns);
        assert!(result.delays_ns[0] < config.min_delay_ns + 1.0);
        assert!(*result.delays_ns.last().unwrap() <= config.max_delay_ns);
        assert!(*result.delays_ns.last().unwrap() > config.max_delay_ns - 1.0);

        // the delay and amplitude of the peak of each spectrum
        let peaks: Vec<(f64, f64)> = result
            .delay_spectrum
            .outer_iter()
            .map(|row| {
                let (idx, amp) = row
                    .iter()
                    .enumerate()
                    .max_by(|a, b| a.1.total_cmp(b.1))
                    .unwrap();
                (result.delays_ns[idx], *amp)
            })
            .collect();
        assert_abs_diff_eq!(peaks[0].0, 250.0, epsilon = 2.0);
        assert_abs_diff_eq!(peaks[1].0, 400.0, epsilon = 2.0);
        assert_abs_diff_eq!(peaks[2].0, 400.0, epsilon = 2.0);
        // twice the ripple, twice the peak
        assert_abs_diff_eq!(peaks[2].1 / peaks[1].1, 2.0, epsilon = 0.01);
        // a ripple of relative amplitude `a` on a spectrum of amplitude `s`
        // peaks at `a * s / 2` times the sum of the window.
        let window_sum = blackman_window(1024).sum();
        assert_abs_diff_eq!(
            peaks[1].1 / (0.2 * 1e6 / 2.0 * window_sum),
            1.0,
            epsilon = 0.01
        );
    }

    /// The transform is linear in the spectrum
    #[test]
    fn test_delay_transform_is_linear() {
        let freq_hz = fine_freqs_hz();
        let spectrum = Array2::from_shape_fn((1, 128), |(_, f)| ((f * 7) % 13) as f64);
        let config = DelayTransformConfig::default();
        let result = delay_transform(&spectrum, &freq_hz, &config).unwrap();
        let result_x3 = delay_transform(&(&spectrum * 3.0), &freq_hz, &config).unwrap();
        for (a, b) in result
            .delay_spectrum
            .iter()
            .zip(result_x3.delay_spectrum.iter())
        {
            assert_abs_diff_eq!(a * 3.0, *b, epsilon = 1e-9);
        }
        // zero delay is the sum of the windowed spectrum
        let windowed_sum: f64 = spectrum
            .row(0)
            .iter()
            .zip(blackman_window(128).iter())
            .map(|(x, w)| x * w)
            .sum();
        assert_abs_diff_eq!(result.delay_spectrum[[0, 0]], windowed_sum, epsilon = 1e-9);
    }

    /// Flagged (NaN) channels are treated as zero, and a spectrum which is
    /// entirely flagged gives NaN.
    #[test]
    fn test_delay_transform_handles_nan() {
        let freq_hz = fine_freqs_hz();
        let mut spectra = Array2::from_shape_fn((3, 128), |(_, f)| 1.0 + (f % 5) as f64);
        // tile 0 has a flagged channel, tile 1 has the same channel zeroed
        spectra[[0, 17]] = f64::NAN;
        spectra[[1, 17]] = 0.0;
        // tile 2 is flagged
        spectra.row_mut(2).fill(f64::NAN);

        let result = delay_transform(&spectra, &freq_hz, &DelayTransformConfig::default()).unwrap();

        assert!(result.delay_spectrum.row(0).iter().all(|x| x.is_finite()));
        assert_eq!(result.delay_spectrum.row(0), result.delay_spectrum.row(1));
        assert!(result.delay_spectrum.row(2).iter().all(|x| x.is_nan()));
    }

    /// When fewer than `MIN_DELAY_BINS` bins are in the requested range, all
    /// the delays up to the maximum are given instead.
    #[test]
    fn test_narrow_delay_range_falls_back_to_broader_range() {
        // 128 channels of 10 kHz give 781.25 ns resolution without padding
        let freq_hz = fine_freqs_hz();
        let spectra = Array2::ones((2, 128));

        // only 2 bins (1562.5, 2343.75) are between 1000 and 3000 ns, so this
        // gives bins 0..=3
        let config = DelayTransformConfig {
            min_delay_ns: 1000.0,
            max_delay_ns: 3000.0,
            target_delay_res_ns: 1000.0,
        };
        let info = calculate_delay_channels(128, &freq_hz, &config);
        let result = delay_transform(&spectra, &freq_hz, &config).unwrap();
        assert_eq!(info.n_delay_channels, 4);
        assert_eq!(
            result.delays_ns.to_vec(),
            vec![0.0, 781.25, 1562.5, 2343.75]
        );
        assert_eq!(result.delay_spectrum.dim(), (2, 4));

        // nothing is in range, and the maximum is beyond the last bin, so this
        // gives every bin.
        let config = DelayTransformConfig {
            min_delay_ns: 1e6,
            max_delay_ns: 2e6,
            target_delay_res_ns: 1000.0,
        };
        let info = calculate_delay_channels(128, &freq_hz, &config);
        let result = delay_transform(&spectra, &freq_hz, &config).unwrap();
        assert_eq!(info.n_delay_channels, 128);
        assert_eq!(result.delays_ns.len(), 128);
        assert_abs_diff_eq!(
            *result.delays_ns.last().unwrap(),
            info.max_computable_delay_ns
        );
    }

    #[test]
    fn test_no_delays_in_range() {
        let freq_hz = fine_freqs_hz();
        let spectra = Array2::ones((1, 128));
        let config = DelayTransformConfig {
            min_delay_ns: -20.0,
            max_delay_ns: -10.0,
            target_delay_res_ns: 1.0,
        };
        let err = delay_transform(&spectra, &freq_hz, &config).unwrap_err();
        assert!(matches!(
            err,
            DelayTransformError::NoDelaysInRange { max_delay_ns } if max_delay_ns < 0.0
        ));
        assert!(err.to_string().contains("no delays in the requested range"));
        // the size calculation never gives zero
        assert_eq!(
            calculate_delay_channels(128, &freq_hz, &config).n_delay_channels,
            1
        );
    }
}
