//! Delay transform utilities for Birli
//!
//! This module provides a delay transform for per-antenna spectra, using
//! a Blackman window and zero-padding to achieve the requested delay resolution.

use marlu::ndarray::{Array1, Array2};
use rustfft::{num_complex::Complex64, FftPlanner};
use std::f64::consts::PI;

/// Configuration for the delay transform.
#[derive(Debug, Clone, Copy)]
pub struct DelayTransformConfig {
    /// Minimum delay to include in output [ns]
    pub min_delay_ns: f64,
    /// Maximum delay to include in output [ns]
    pub max_delay_ns: f64,
    /// Target delay resolution [ns]
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
    /// Delay spectrum (nants, ndelays)
    pub delay_spectrum: Array2<f64>,
    /// Delay axis [ns]
    pub delays_ns: Array1<f64>,
}

/// Information about delay transform dimensions
#[derive(Debug)]
pub struct DelayChannelInfo {
    /// Number of delay channels in the output
    pub n_delay_channels: usize,
    /// Zero-padded frequency channels
    pub nfreqs_padded: usize,
    /// Actual delay resolution achieved [ns]
    pub actual_delay_res_ns: f64,
    /// Maximum computable delay [ns]
    pub max_computable_delay_ns: f64,
}

/// Calculate the number of delay channels that will be in the transform result
pub fn calculate_delay_channels(
    nfreqs: usize,
    freq_hz: &Array1<f64>,
    config: &DelayTransformConfig,
) -> DelayChannelInfo {
    // Calculate zero-padding factor for better delay resolution
    let freq_res_hz = freq_hz[1] - freq_hz[0];

    // Zero-padding factor to reach target delay resolution
    let pad_factor =
        ((1e9 / (config.target_delay_res_ns * freq_res_hz * nfreqs as f64)).ceil() as usize).max(1);
    let nfreqs_padded = nfreqs * pad_factor;

    // Actual delay resolution and range
    let delay_res_s = 1.0 / (nfreqs_padded as f64 * freq_res_hz);
    let actual_delay_res_ns = delay_res_s * 1e9;
    let max_computable_delay_ns = (nfreqs_padded - 1) as f64 * actual_delay_res_ns;

    // Create delay axis to determine filtering
    let mut delays_ns = Array1::zeros(nfreqs_padded);
    for i in 0..nfreqs_padded {
        delays_ns[i] = i as f64 * actual_delay_res_ns;
    }

    // Determine final number of channels based on filtering logic
    let delay_mask: Vec<bool> = delays_ns
        .iter()
        .map(|&d| d >= config.min_delay_ns && d <= config.max_delay_ns)
        .collect();
    let mask_count = delay_mask.iter().filter(|&&x| x).count();

    let n_delay_channels = if mask_count > 5 {
        // Will use filtered range
        mask_count
    } else {
        // Will use broader range
        if max_computable_delay_ns < config.max_delay_ns {
            // Show full computed range
            nfreqs_padded
        } else {
            // Show subset up to max_delay_ns
            delays_ns
                .iter()
                .take_while(|&&d| d <= config.max_delay_ns)
                .count()
                .max(1) // Ensure at least 1 channel
        }
    };

    DelayChannelInfo {
        n_delay_channels,
        nfreqs_padded,
        actual_delay_res_ns,
        max_computable_delay_ns,
    }
}

/// Blackman window (same as numpy.blackman)
fn blackman_window(n: usize) -> Array1<f64> {
    let mut window = Array1::zeros(n);
    if n == 1 {
        window[0] = 1.0;
        return window;
    }
    for i in 0..n {
        let a = 2.0 * PI * i as f64 / (n as f64 - 1.0);
        window[i] = 0.42 - 0.5 * a.cos() + 0.08 * (2.0 * a).cos();
    }
    window
}

/// Perform a delay transform (FFT along frequency axis) for each antenna.
/// Input: rdx_ant_freq (nants, nfreqs), freq_hz (nfreqs)
pub fn delay_transform(
    rdx_ant_freq: &Array2<f64>,
    freq_hz: &Array1<f64>,
    config: &DelayTransformConfig,
) -> Result<DelayTransformResult, Box<dyn std::error::Error>> {
    let nants = rdx_ant_freq.nrows();
    let nfreqs = rdx_ant_freq.ncols();

    eprintln!(
        "Computing delay transform: {}-{} ns",
        config.min_delay_ns, config.max_delay_ns
    );

    // Get delay channel information
    let delay_info = calculate_delay_channels(nfreqs, freq_hz, config);
    let nfreqs_padded = delay_info.nfreqs_padded;

    eprintln!(
        "Zero-padding: {} -> {} for {:.1} ns resolution",
        nfreqs, nfreqs_padded, delay_info.actual_delay_res_ns
    );

    // Initialize output array
    let mut rdx_ant_delay = Array2::from_elem((nants, nfreqs_padded), f64::NAN);

    // FFT setup
    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(nfreqs_padded);

    // Blackman window
    let window = blackman_window(nfreqs);

    // Process each antenna
    for ant in 0..nants {
        let ant_data = rdx_ant_freq.row(ant);
        let valid_mask: Vec<bool> = ant_data.iter().map(|&x| x.is_finite()).collect();
        if !valid_mask.iter().any(|&x| x) {
            continue;
        }
        let mut data_windowed = Array1::zeros(nfreqs);
        for i in 0..nfreqs {
            data_windowed[i] = if valid_mask[i] {
                ant_data[i] * window[i]
            } else {
                0.0
            };
        }
        let mut data_padded: Vec<Complex64> = vec![Complex64::new(0.0, 0.0); nfreqs_padded];
        for i in 0..nfreqs {
            data_padded[i] = Complex64::new(data_windowed[i], 0.0);
        }
        fft.process(&mut data_padded);
        for i in 0..nfreqs_padded {
            rdx_ant_delay[[ant, i]] = data_padded[i].norm();
        }
    }

    // Delay axis [ns]
    let mut delays_ns = Array1::zeros(nfreqs_padded);
    for i in 0..nfreqs_padded {
        delays_ns[i] = i as f64 * delay_info.actual_delay_res_ns;
    }

    eprintln!(
        "Actual delay range: 0 to {:.1} ns, res={:.2} ns",
        delay_info.max_computable_delay_ns, delay_info.actual_delay_res_ns
    );
    eprintln!(
        "Channel BW = {:.3} MHz -> delay period = {:.1} ns",
        (freq_hz[1] - freq_hz[0]) / 1e6,
        1e9 / (freq_hz[1] - freq_hz[0])
    );

    // Check if we have finite data
    let has_finite = rdx_ant_delay.iter().any(|&x| x.is_finite());
    eprintln!(
        "Delay data shape: {:?}, has finite data: {}",
        rdx_ant_delay.dim(),
        has_finite
    );

    if has_finite {
        let min_val = rdx_ant_delay
            .iter()
            .filter(|&&x| x.is_finite())
            .fold(f64::INFINITY, |a, &b| a.min(b));
        let max_val = rdx_ant_delay
            .iter()
            .filter(|&&x| x.is_finite())
            .fold(f64::NEG_INFINITY, |a, &b| a.max(b));
        eprintln!("Data range: {:.2e} to {:.2e}", min_val, max_val);
    }

    // Filter to requested delay range (same logic as in calculate_delay_channels)
    let delay_mask: Vec<bool> = delays_ns
        .iter()
        .map(|&d| d >= config.min_delay_ns && d <= config.max_delay_ns)
        .collect();
    let mask_count = delay_mask.iter().filter(|&&x| x).count();

    eprintln!(
        "Requested range {}-{} ns has {} points",
        config.min_delay_ns, config.max_delay_ns, mask_count
    );

    let (filtered_delay_spectrum, filtered_delays) = if mask_count > 5 {
        // Apply the mask
        let valid_indices: Vec<usize> = delay_mask
            .iter()
            .enumerate()
            .filter(|(_, &mask)| mask)
            .map(|(i, _)| i)
            .collect();
        let mut filtered_spectrum = Array2::zeros((nants, valid_indices.len()));
        let mut filtered_delays_array = Array1::zeros(valid_indices.len());
        for (new_idx, &old_idx) in valid_indices.iter().enumerate() {
            for ant in 0..nants {
                filtered_spectrum[[ant, new_idx]] = rdx_ant_delay[[ant, old_idx]];
            }
            filtered_delays_array[new_idx] = delays_ns[old_idx];
        }

        eprintln!(
            "Using filtered range: {:.1} to {:.1} ns",
            filtered_delays_array[0],
            filtered_delays_array[filtered_delays_array.len() - 1]
        );

        (filtered_spectrum, filtered_delays_array)
    } else {
        // Show broader range if requested range is too narrow
        if delay_info.max_computable_delay_ns < config.max_delay_ns {
            eprintln!(
                "Showing full computed range: 0 to {:.1} ns",
                delay_info.max_computable_delay_ns
            );
            (rdx_ant_delay, delays_ns)
        } else {
            // Show subset up to max_delay_ns
            let subset_indices: Vec<usize> = delays_ns
                .iter()
                .enumerate()
                .take_while(|(_, &d)| d <= config.max_delay_ns)
                .map(|(i, _)| i)
                .collect();
            if !subset_indices.is_empty() {
                let subset_end = subset_indices.len();
                let mut subset_spectrum = Array2::zeros((nants, subset_end));
                let mut subset_delays = Array1::zeros(subset_end);
                for i in 0..subset_end {
                    for ant in 0..nants {
                        subset_spectrum[[ant, i]] = rdx_ant_delay[[ant, i]];
                    }
                    subset_delays[i] = delays_ns[i];
                }

                eprintln!(
                    "Showing subset: 0 to {:.1} ns",
                    subset_delays[subset_end - 1]
                );
                (subset_spectrum, subset_delays)
            } else {
                eprintln!("No delays <= {} ns found!", config.max_delay_ns);
                return Err("No valid delays in requested range".into());
            }
        }
    };

    eprintln!(
        "Final result: {:.1} to {:.1} ns, {} points",
        filtered_delays[0],
        filtered_delays[filtered_delays.len() - 1],
        filtered_delays.len()
    );

    Ok(DelayTransformResult {
        delay_spectrum: filtered_delay_spectrum,
        delays_ns: filtered_delays,
    })
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use super::*;

    #[test]
    fn test_blackman_window() {
        let window = blackman_window(10);
        assert_eq!(window.len(), 10);
        // Blackman window should start and end near zero
        assert!(window[0] < 0.1);
        assert!(window[9] < 0.1);
        // And peak in the middle
        assert!(window[4] > 0.3);
        assert!(window[5] > 0.3);
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
        assert_eq!(info.n_delay_channels, 109);
        assert!(info.nfreqs_padded >= nfreqs);
        assert_abs_diff_eq!(info.actual_delay_res_ns, 0.5, epsilon = 1e-1);
        assert!(info.max_computable_delay_ns >= 60.0);
    }

    #[test]
    fn test_delay_transform_shape() {
        let nants = 3;
        let nfreqs = 64;
        let rdx_ant_freq = Array2::ones((nants, nfreqs));
        let freq_hz = Array1::linspace(100e6, 200e6, nfreqs);
        let config = DelayTransformConfig::default();

        // Test that calculate_delay_channels predicts the correct size
        let info = calculate_delay_channels(nfreqs, &freq_hz, &config);
        let result = delay_transform(&rdx_ant_freq, &freq_hz, &config).unwrap();

        assert_eq!(result.delay_spectrum.nrows(), nants);
        assert_eq!(result.delay_spectrum.ncols(), info.n_delay_channels);
        assert_eq!(result.delays_ns.len(), info.n_delay_channels);
    }

    #[test]
    fn test_uvfits_like_parameters() {
        // Test with parameters similar to the UVfits file that Python processes
        let nants = 3;
        let nfreqs = 128;

        // MWA-like frequency setup: 133-156 MHz, 10 kHz channels
        let start_freq = 133e6; // 133 MHz
        let freq_res = 10e3; // 10 kHz
        let freq_hz = Array1::from_iter((0..nfreqs).map(|i| start_freq + i as f64 * freq_res));

        // Create some realistic test data with different spectral signatures per antenna
        let mut rdx_ant_freq = Array2::zeros((nants, nfreqs));
        for ant in 0..nants {
            for freq in 0..nfreqs {
                // Simulate different antenna responses
                let normalized_freq = freq as f64 / nfreqs as f64;
                rdx_ant_freq[[ant, freq]] = match ant {
                    0 => 1e6 * (1.0 + 0.1 * (normalized_freq * 4.0 * std::f64::consts::PI).sin()), // Strong sine pattern
                    1 => 5e5 * (1.0 + 0.05 * (normalized_freq * 8.0 * std::f64::consts::PI).sin()), // Weaker, higher freq
                    _ => 2e6 * (1.0 + 0.2 * (normalized_freq * 2.0 * std::f64::consts::PI).cos()), // Cosine pattern
                };
            }
        }

        let config = DelayTransformConfig {
            min_delay_ns: 10.0,
            max_delay_ns: 1000.0,
            target_delay_res_ns: 1.0,
        };

        let result = delay_transform(&rdx_ant_freq, &freq_hz, &config).unwrap();

        eprintln!("UVfits-like test results:");
        eprintln!("Input shape: {:?}", rdx_ant_freq.dim());
        eprintln!("Output shape: {:?}", result.delay_spectrum.dim());
        eprintln!(
            "Delay range: {:.1} to {:.1} ns",
            result.delays_ns[0],
            result.delays_ns[result.delays_ns.len() - 1]
        );
        eprintln!(
            "Input data range: {:.2e} to {:.2e}",
            rdx_ant_freq.iter().fold(f64::INFINITY, |a, &b| a.min(b)),
            rdx_ant_freq
                .iter()
                .fold(f64::NEG_INFINITY, |a, &b| a.max(b))
        );
        eprintln!(
            "Output data range: {:.2e} to {:.2e}",
            result
                .delay_spectrum
                .iter()
                .filter(|&&x| x.is_finite())
                .fold(f64::INFINITY, |a, &b| a.min(b)),
            result
                .delay_spectrum
                .iter()
                .filter(|&&x| x.is_finite())
                .fold(f64::NEG_INFINITY, |a, &b| a.max(b))
        );

        // Basic sanity checks
        assert_eq!(result.delay_spectrum.nrows(), nants);
        assert!(result.delay_spectrum.ncols() > 100); // Should have reasonable number of delay channels
        assert!(result.delays_ns[0] >= 10.0);
        assert!(result.delays_ns[result.delays_ns.len() - 1] <= 1000.0);

        // Check that we have finite, non-zero data
        let has_finite = result
            .delay_spectrum
            .iter()
            .any(|&x| x.is_finite() && x > 0.0);
        assert!(
            has_finite,
            "Delay spectrum should contain finite, positive values"
        );
    }

    #[test]
    #[ignore] // Run with: cargo test test_real_data_debug -- --ignored --nocapture
    fn test_real_data_debug() {
        use crate::cli::BirliContext;
        use std::path::Path;

        let metafits_path = "/Users/dev/Code/mwa-demo/demo/data/1361397000/raw/1361397000.metafits";
        let gpubox_path = "/Users/dev/Code/mwa-demo/demo/data/1361397000/raw/1361397000_20230225214942_ch137_000.fits";
        let output_path = "test_debug_metrics.fits";

        // Check if input files exist
        if !Path::new(metafits_path).exists() || !Path::new(gpubox_path).exists() {
            eprintln!("Skipping test - input files not found");
            eprintln!("Expected: {} and {}", metafits_path, gpubox_path);
            return;
        }

        // Use BirliContext following the pattern from existing tests
        let args = vec![
            "birli",
            "--metrics-out",
            output_path,
            "-m",
            metafits_path,
            "--no-draw-progress",
            gpubox_path,
        ];

        let birli_ctx = BirliContext::from_args(&args).unwrap();
        birli_ctx.run().unwrap();

        // Check if output file was created
        assert!(
            Path::new(output_path).exists(),
            "Output metrics file was not created"
        );

        eprintln!("Successfully created {}", output_path);
        eprintln!("You can view it with: carta {}", output_path);
        eprintln!("Debug output should show spectral analysis above");

        // Clean up
        std::fs::remove_file(output_path).ok();
    }
}
