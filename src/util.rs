//! Utility functions

use crate::cxx_aoflagger::ffi::{CxxFlagMask, CxxImageSet};
use cxx::UniquePtr;
use std::{cmp::min, fmt::Write};

/// Peek into a baseline imgset for debugging
///
/// # Examples
///
/// ```rust
/// use birli::{context_to_baseline_imgsets, cxx_aoflagger_new, util::dump_baseline_imgsets, mwalib};
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
/// let baseline_imgsets = unsafe {
///     let aoflagger = cxx_aoflagger_new();
///     context_to_baseline_imgsets(
///         &aoflagger,
///         &context,
///         &context.common_coarse_chan_indices.clone(),
///         &context.common_timestep_indices.clone(),
///         None,
///     )
/// };
///
/// dump_baseline_imgsets(&baseline_imgsets, Some(2), Some(2), Some(2), Some(16));
/// ```
pub fn dump_baseline_imgsets(
    baseline_imgsets: &[UniquePtr<CxxImageSet>],
    baseline_limit: Option<usize>,
    freq_limit: Option<usize>,
    timestep_limit: Option<usize>,
    radix: Option<usize>,
) -> String {
    let mut out = String::new();
    let baseline_limit = baseline_limit.unwrap_or(baseline_imgsets.len());
    baseline_imgsets
        .iter()
        .take(baseline_limit)
        .enumerate()
        .for_each(|(baseline_idx, imgset)| {
            writeln!(&mut out, "baseline {}", &baseline_idx).unwrap();
            dump_imgset(imgset, freq_limit, timestep_limit, radix);
        });
    out
}

/// Peek into a single imgset for debugging
///
/// # Examples
///
/// ```rust
/// use birli::{context_to_baseline_imgsets, cxx_aoflagger_new, util::dump_imgset, mwalib};
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
/// let baseline_imgsets = unsafe {
///     let aoflagger = cxx_aoflagger_new();
///     context_to_baseline_imgsets(
///         &aoflagger,
///         &context,
///         &context.common_coarse_chan_indices.clone(),
///         &context.common_timestep_indices.clone(),
///         None,
///     )
/// };
/// // if log_enabled!(Level::Debug) {
/// baseline_imgsets.iter().take(128).enumerate().for_each(|(baseline_idx, imgset)| {
///     let baseline = &context.metafits_context.baselines[baseline_idx];
///     if baseline.ant1_index == 0 && baseline.ant2_index == 2 {
///         println!("baseline {:04} ({}, {})\n{}", baseline_idx, baseline.ant1_index, baseline.ant2_index, dump_imgset(imgset, Some(128), Some(128), Some(16)));
///     }
/// });
/// ```
pub fn dump_imgset(
    imgset: &UniquePtr<CxxImageSet>,
    freq_limit: Option<usize>,
    timestep_limit: Option<usize>,
    radix: Option<usize>,
) -> String {
    let mut out = String::new();
    let stride = imgset.HorizontalStride();
    let width = imgset.Width();
    let height = imgset.Height();
    let count = imgset.ImageCount();
    let freq_limit = min(freq_limit.unwrap_or(height), height);
    let timestep_limit = min(timestep_limit.unwrap_or(width), width);
    (0..count).for_each(|pol| {
        let img_buf = imgset.ImageBuffer(pol);
        img_buf
            .chunks(stride)
            .into_iter()
            .take(freq_limit)
            .enumerate()
            .for_each(|(freq_idx, chunk)| {
                write!(&mut out, "pol {:02} freq {:03} | ", &pol, &freq_idx).unwrap();
                chunk.iter().take(timestep_limit).for_each(|&fl| {
                    let radix = radix.unwrap_or(0);
                    let vis = fl as i32;
                    if radix == 2 {
                        write!(&mut out, "{:032b}", vis).unwrap();
                    } else if radix == 16 {
                        write!(
                            &mut out,
                            "{}0x{:08x}",
                            if vis < 0 { "-" } else { "+" },
                            if vis < 0 { -vis } else { vis }
                        )
                        .unwrap();
                    }
                    write!(&mut out, " ").unwrap();
                });
                writeln!(&mut out).unwrap();
            });
    });
    out
}

/// Peek into a baseline flagmask vector for debugging
///
/// # Examples
///
/// ```rust
/// use birli::{context_to_baseline_imgsets, cxx_aoflagger_new, flag_imgsets, util::dump_baseline_flagmasks, mwalib};
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
/// let baseline_flagmasks = unsafe {
///     let aoflagger = cxx_aoflagger_new();
///     let baseline_imgsets = context_to_baseline_imgsets(
///         &aoflagger,
///         &context,
///         &context.common_coarse_chan_indices.clone(),
///         &context.common_timestep_indices.clone(),
///         None,
///     );
///     let strategy_filename = &aoflagger.FindStrategyFileMWA();
///     flag_imgsets(&aoflagger, &strategy_filename, baseline_imgsets)
/// };
///
/// dump_baseline_flagmasks(&baseline_flagmasks, Some(2), Some(2), Some(16));
/// ```
pub fn dump_baseline_flagmasks(
    baseline_flagmasks: &[UniquePtr<CxxFlagMask>],
    baseline_limit: Option<usize>,
    freq_limit: Option<usize>,
    timestep_limit: Option<usize>,
) -> String {
    let mut out = String::new();
    let baseline_limit = baseline_limit.unwrap_or(baseline_flagmasks.len());
    baseline_flagmasks
        .iter()
        .take(baseline_limit)
        .enumerate()
        .for_each(|(baseline_idx, flagmask)| {
            writeln!(&mut out, "baseline {}", &baseline_idx).unwrap();
            dump_flagmask(flagmask, freq_limit, timestep_limit);
        });
    out
}

/// Peek into a flagmask for debugging
///
/// # Examples
///
/// ```rust
/// use birli::{context_to_baseline_imgsets, cxx_aoflagger_new, flag_imgsets, util::dump_flagmask, mwalib};
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
/// let baseline_flagmasks = unsafe {
///     let aoflagger = cxx_aoflagger_new();
///     let baseline_imgsets = context_to_baseline_imgsets(
///         &aoflagger,
///         &context,
///         &context.common_coarse_chan_indices.clone(),
///         &context.common_timestep_indices.clone(),
///         None,
///     );
///     let strategy_filename = &aoflagger.FindStrategyFileMWA();
///     flag_imgsets(&aoflagger, &strategy_filename, baseline_imgsets)
/// };
///
/// baseline_flagmasks.iter().take(128).enumerate().for_each(|(baseline_idx, flagmask)| {
///     let baseline = &context.metafits_context.baselines[baseline_idx];
///     if baseline.ant1_index == 0 && baseline.ant2_index == 2 {
///         println!("baseline {:04} ({}, {})\n{}", baseline_idx, baseline.ant1_index, baseline.ant2_index, dump_flagmask(flagmask, Some(128), Some(128)));
///     }
/// });
/// ```
pub fn dump_flagmask(
    flagmask: &UniquePtr<CxxFlagMask>,
    freq_limit: Option<usize>,
    timestep_limit: Option<usize>,
) -> String {
    let mut out = String::new();
    let flag_buf = flagmask.Buffer();
    let stride = flagmask.HorizontalStride();
    let width = flagmask.Width();
    let height = flagmask.Height();
    let freq_limit = freq_limit.unwrap_or(height);
    let timestep_limit = timestep_limit.unwrap_or(width);
    flag_buf
        .chunks(stride)
        .into_iter()
        .take(freq_limit)
        .enumerate()
        .for_each(|(freq_idx, chunk)| {
            write!(&mut out, "freq {:03} | ", &freq_idx).unwrap();
            chunk.iter().take(timestep_limit).for_each(|&fl| {
                let symbol = if fl as i8 > 0 { "#" } else { "." };
                write!(&mut out, "{} ", symbol).unwrap();
            });
            writeln!(&mut out).unwrap();
        });
    out
}

#[cfg(test)]
mod tests {
    use super::{dump_baseline_flagmasks, dump_baseline_imgsets};
    use crate::{
        context_to_baseline_imgsets, cxx_aoflagger_new, flag_imgsets, get_flaggable_timesteps,
    };
    use mwa_rust_core::mwalib;
    use mwalib::CorrelatorContext;

    // TODO: deduplicate this from lib.rs
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
    fn test_dump_baseline_imgsets_mwa_ord() {
        let aoflagger = unsafe { cxx_aoflagger_new() };
        let context = get_mwa_ord_context();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;

        let baseline_imgsets = context_to_baseline_imgsets(
            &aoflagger,
            &context,
            img_coarse_chan_idxs,
            &img_timestep_idxs,
            None,
        );

        dump_baseline_imgsets(&baseline_imgsets, Some(2), Some(2), Some(2), Some(16));
    }

    #[test]
    fn test_dump_baseline_flagmasks_mwax() {
        let aoflagger = unsafe { cxx_aoflagger_new() };
        let context = get_mwax_context();

        let img_timestep_idxs = get_flaggable_timesteps(&context).unwrap();
        let img_coarse_chan_idxs = &context.common_coarse_chan_indices;

        let baseline_imgsets = context_to_baseline_imgsets(
            &aoflagger,
            &context,
            img_coarse_chan_idxs,
            &img_timestep_idxs,
            None,
        );
        let strategy_filename = &aoflagger.FindStrategyFileMWA();
        let baseline_flagmasks = flag_imgsets(&aoflagger, strategy_filename, baseline_imgsets);

        dump_baseline_flagmasks(&baseline_flagmasks, Some(2), Some(2), Some(16));
    }
}
