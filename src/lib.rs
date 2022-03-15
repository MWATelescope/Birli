#![warn(missing_docs)]
#![warn(rustdoc::missing_doc_code_examples)]
#![warn(clippy::missing_safety_doc)]
#![warn(clippy::missing_errors_doc)]

//! Birli is a library of common preprocessing tasks performed in the data pipeline of the Murchison
//! Widefield Array (MWA) Telescope.
//!
//! # Examples
//!
//! Here's an example of how to flag some visibility files
//!
//! ```rust
//! use birli::{
//!     context_to_jones_array, write_flags,
//!     mwalib::CorrelatorContext, write_uvfits,
//!     add_dimension, get_weight_factor, flag_to_weight_array,
//!     FlagContext, VisSelection
//! };
//! use tempfile::tempdir;
//!
//! // define our input files
//! let metafits_path = "tests/data/1297526432_mwax/1297526432.metafits";
//! let gpufits_paths = vec![
//!     "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_000.fits",
//!     "tests/data/1297526432_mwax/1297526432_20210216160014_ch117_001.fits",
//!     "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_000.fits",
//!     "tests/data/1297526432_mwax/1297526432_20210216160014_ch118_001.fits",
//! ];
//!
//! // define a temporary directory for output files
//! let tmp_dir = tempdir().unwrap();
//!
//! // define our output paths
//! let flag_template = tmp_dir.path().join("Flagfile%%%.mwaf");
//! let uvfits_out = tmp_dir.path().join("1297526432.uvfits");
//!
//! // Create an mwalib::CorrelatorContext for accessing visibilities.
//! let corr_ctx = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
//!
//! // Determine which timesteps and coarse channels we want to use
//! let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
//!
//! // Prepare our flagmasks with known bad antennae
//! let flag_ctx = FlagContext::from_mwalib(&corr_ctx);
//! let flag_array = flag_ctx.to_array(
//!     &vis_sel.timestep_range,
//!     &vis_sel.coarse_chan_range,
//!     vis_sel.get_ant_pairs(&corr_ctx.metafits_context)
//! ).unwrap();
//!
//! // load visibilities into our array of jones matrices
//! let (mut jones_array, flag_array) = context_to_jones_array(
//!     &corr_ctx,
//!     &vis_sel.timestep_range,
//!     &vis_sel.coarse_chan_range,
//!     Some(flag_array),
//!     false,
//! ).unwrap();
//!
//! // write the flags to disk as .mwaf
//! write_flags(&corr_ctx, &flag_array, flag_template.to_str().unwrap(), &vis_sel.coarse_chan_range).unwrap();
//! // write the visibilities to disk as .uvfits
//!
//! let num_pols = corr_ctx.metafits_context.num_visibility_pols;
//! let flag_array = add_dimension(flag_array.view(), num_pols);
//! let weight_factor = get_weight_factor(&corr_ctx);
//! let weight_array = flag_to_weight_array(flag_array.view(), weight_factor);
//! write_uvfits(
//!     uvfits_out.as_path(),
//!     &corr_ctx,
//!     jones_array.view(),
//!     weight_array.view(),
//!     flag_array.view(),
//!     &vis_sel.timestep_range,
//!     &vis_sel.coarse_chan_range,
//!     &vis_sel.baseline_idxs,
//!     None,
//!     None,
//!     1,
//!     1,
//!     false,
//! ).unwrap();
//! ```
//!
//! # Details
//!
//! Birli reads visibilities with [`MWALib`] and uses CXX to bind to the [`AOFlagger`] C++ library.
//! For more details about AOFlagger's interface, check out the [`aoflagger::AOFlagger`]
//! documentation
//!
//! [`MWALib`]: https://github.com/MWATelescope/mwalib
//! [`AOFlagger`]: https://gitlab.com/aroffringa/aoflagger
//! [`aoflagger::AOFlagger`]: http://www.andreoffringa.org/aoflagger/doxygen/classaoflagger_1_1AOFlagger.html

use cfg_if::cfg_if;
use log::{trace, warn};
cfg_if! {
    if #[cfg(feature = "aoflagger")] {
        use ndarray::ArrayView2;
    }
}
use std::ops::Range;

pub mod io;
pub use io::{mwaf::FlagFileSet, uvfits::UvfitsWriter, write_ms, write_uvfits};
pub mod corrections;
pub use corrections::{correct_cable_lengths, correct_geometry};
pub mod calibration;
pub mod flags;
#[cfg(test)]
pub use approx;
#[cfg(test)]
pub(crate) mod types;
pub use flags::{add_dimension, flag_to_weight_array, get_weight_factor, write_flags, FlagContext};
pub mod passband_gains;
pub use marlu;
pub use marlu::{
    mwalib,
    mwalib::{fitsio, fitsio_sys, CorrelatorContext},
    ndarray,
    ndarray::{parallel::prelude::*, Array3, Axis},
    Complex, Jones,
};
#[cfg(test)]
pub(crate) use types::TestJones;

mod error;
pub use error::BirliError;

pub mod preprocessing;
pub use preprocessing::PreprocessContext;

pub mod selection;
pub use selection::VisSelection;

pub mod cli;
pub use cli::BirliContext;

cfg_if! {
    if #[cfg(feature = "aoflagger")] {
        pub use flags::{flag_jones_array, flag_jones_array_existing};
        pub use aoflagger_sys::{cxx_aoflagger_new, CxxAOFlagger, CxxFlagMask, UniquePtr, CxxImageSet};
        use ndarray::{ArrayBase, Dim, ViewRepr};
    }
}

#[cfg(test)]
pub mod test_common;

#[macro_export]
/// Time a statement and increment the timer given by name in the hashmap of durations
macro_rules! with_increment_duration {
    ($durs:expr, $name:literal, $($s:stmt);+ $(;)?) => {
        {
            let _now = std::time::Instant::now();
            let _res = {
                $(
                    $s
                )*
            };
            *$durs.entry($name.into()).or_insert(Duration::default()) += _now.elapsed();
            _res
        }
    };
}

/// generate a 3 dimensional array of Jones matrices from an observation's
/// [`mwalib::CorrelatorContext`], for all baselines, over a given range of
/// mwalib timestep and coarse channel indices.
///
/// The dimensions of the array are:
///  - timestep
///  - channel
///  - baseline
///
/// An equally sized flag array is also returned with flags indicating when reading via mwalib
/// causes a GPUBoxError.
///
/// # Details:
///
/// Observations can sometimes be too large to fit in memory. This method will only load
/// visibilities from the provided timesteps and coarse channels, in order to enable the visibility to
/// be read in user-defined "chunks" of time or frequency.
///
/// The timesteps are indices in the [`mwalib::CorrelatorContext`]'s timestep array, which should be a contiguous
/// superset of times from all provided coarse gpubox files. A similar concept applies to coarse channels.
/// Instead of reading visibilities for all known timesteps / coarse channels, it is recommended to use
/// `common_coarse_chan_indices` and `common_timestep_indices`, as these ignore timesteps and coarse channels
/// which are missing contiguous data. `common_good_timestep_indices` is also a good choice to avoid quack time.
///
/// For more details, see the [documentation](https://docs.rs/mwalib/latest/mwalib/struct.CorrelatorContext.html).
///
/// Note: it doesn't make sense to ask aoflagger to flag non-contiguous timesteps
/// or coarse channels, and so this interface only allows to ranges to be used.
/// For flagging an obeservation with "picket fence" coarse channels or timesteps,
/// contiguous ranges should be flagged separately.
///
/// [`mwalib::CorrelatorContext`]: https://docs.rs/mwalib/latest/mwalib/struct.CorrelatorContext.html
///
/// # Examples
///
/// ```rust
/// use birli::{context_to_jones_array, mwalib::CorrelatorContext, VisSelection};
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
/// let img_timestep_idxs = &corr_ctx.common_timestep_indices;
/// let good_timestep_idxs = &corr_ctx.common_good_timestep_indices;
///
/// let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
/// vis_sel.timestep_range =
///     *img_timestep_idxs.first().unwrap()..(*img_timestep_idxs.last().unwrap() + 1);
///
///
/// // read visibilities out of the gpubox files
/// let (jones_array, _) = context_to_jones_array(
///     &corr_ctx,
///     &vis_sel.timestep_range,
///     &vis_sel.coarse_chan_range,
///     None,
///     false,
/// ).unwrap();
///
/// let dims_common = jones_array.dim();
///
/// // now try only with good timesteps
/// vis_sel.timestep_range =
///     *good_timestep_idxs.first().unwrap()..(*good_timestep_idxs.last().unwrap() + 1);
///
/// // read visibilities out of the gpubox files
/// let (jones_array, _) = context_to_jones_array(
///     &corr_ctx,
///     &vis_sel.timestep_range,
///     &vis_sel.coarse_chan_range,
///     None,
///     false,
/// ).unwrap();
///
/// let dims_good = jones_array.dim();
///
/// assert_ne!(dims_common, dims_good);
/// ```
///
/// # Errors
///
/// can throw BadArrayShape flag_array is provided and its shape does not match
/// that of the timestep and coarse channel ranges.
// #[deprecated(since = "0.6.0", note = "use `VisSelection::read_mwalib` instead")]
pub fn context_to_jones_array(
    corr_ctx: &CorrelatorContext,
    timestep_range: &Range<usize>,
    coarse_chan_range: &Range<usize>,
    flag_array: Option<Array3<bool>>,
    draw_progress: bool,
) -> Result<(Array3<Jones<f32>>, Array3<bool>), BirliError> {
    trace!("start context_to_jones_array");

    let vis_sel = VisSelection {
        timestep_range: timestep_range.clone(),
        coarse_chan_range: coarse_chan_range.clone(),
        ..VisSelection::from_mwalib(corr_ctx)?
    };

    // allocate our result

    let num_timesteps = timestep_range.len();
    let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
    let num_coarse_chans = coarse_chan_range.len();
    let num_chans = num_coarse_chans * fine_chans_per_coarse;
    let num_baselines = corr_ctx.metafits_context.num_baselines;

    let shape = (num_timesteps, num_chans, num_baselines);
    let num_elems = num_timesteps * num_chans * num_baselines;

    // We need this many gibibytes to do processing (visibilities and flags).
    let need_gib = (shape.0 * shape.1 * shape.2)
        * (std::mem::size_of::<Jones<f32>>() + std::mem::size_of::<bool>())
        / 1024_usize.pow(3);

    let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse)?;
    let mut flag_array = if let Some(flag_array_) = flag_array {
        if flag_array_.dim() != shape {
            return Err(BirliError::BadArrayShape {
                argument: "flag_array".to_string(),
                function: "context_to_jones_array".to_string(),
                expected: format!("({}, {}, {}, 4)", shape.0, shape.1, shape.2),
                received: format!("{:?}", shape),
            });
        };
        flag_array_
    } else {
        let mut v = Vec::new();
        match v.try_reserve_exact(shape.0 * shape.1 * shape.2) {
            Ok(()) => {
                v.resize(num_elems, false);
                Array3::from_shape_vec(shape, v).unwrap()
            }
            // Instead of erroring out with how many GiB we need for *this*
            // array, error out with how many we need total.
            Err(_) => return Err(BirliError::InsufficientMemory { need_gib }),
        }
    };

    vis_sel.read_mwalib(corr_ctx, &mut jones_array, &mut flag_array, draw_progress)?;

    Ok((jones_array, flag_array))
}

/// Create an aoflagger [`CxxImageSet`] for a particular baseline from the given jones array
///
/// # Assumptions
///
/// - `baseline_jones_view` is [timestep][channel] for one baseline
/// - imageset is timesteps wide, and channels high
/// - jones matrics are always XX, YY, XY, YX
///
/// # Errors
///
/// Will throw [`BirliError`] if there was an error reading from corr_ctx.
#[cfg(feature = "aoflagger")]
pub fn jones_baseline_view_to_imageset(
    aoflagger: &CxxAOFlagger,
    // jones_array: &Array3<Jones<f32>>,
    baseline_jones_view: &ArrayBase<ViewRepr<&Jones<f32>>, Dim<[usize; 2]>>,
) -> Result<UniquePtr<CxxImageSet>, BirliError> {
    let array_dims = baseline_jones_view.dim();
    let img_count = 8;
    let imgset = unsafe {
        aoflagger.MakeImageSet(
            array_dims.0,
            array_dims.1,
            img_count,
            0 as f32,
            array_dims.0,
        )
    };
    let img_stride = imgset.HorizontalStride();
    let mut img_bufs: Vec<&mut [f32]> = (0..img_count)
        .map(|img_idx| unsafe { imgset.ImageBufferMutUnsafe(img_idx) })
        .collect();

    // TODO: benchmark if iterate over pol first
    for (img_timestep_idx, timestep_jones_view) in baseline_jones_view.outer_iter().enumerate() {
        for (img_chan_idx, singular_jones_view) in timestep_jones_view.outer_iter().enumerate() {
            let jones = singular_jones_view.get(()).unwrap();
            for (img_idx, img_buf) in img_bufs.iter_mut().enumerate() {
                let pol_idx = img_idx / 2;
                img_buf[img_chan_idx * img_stride + img_timestep_idx] = if img_idx % 2 == 0 {
                    jones[pol_idx].re
                } else {
                    jones[pol_idx].im
                };
            }
        }
    }

    Ok(imgset)
}

/// Create an aoflagger [`CxxFlagMask`] for a from the given flag array view
///
/// # Assumptions
///
/// - flag array view is [timestep][channel] for one baseline
/// - flagmask is timesteps wide, and channels high
///
/// # Errors
///
/// Will throw [`BirliError`] if there was an error reading from corr_ctx.
#[cfg(feature = "aoflagger")]
pub fn flag_baseline_view_to_flagmask(
    aoflagger: &CxxAOFlagger,
    baseline_flag_view: &ArrayView2<bool>,
) -> Result<UniquePtr<CxxFlagMask>, BirliError> {
    let array_dims = baseline_flag_view.dim();
    let mut flag_mask = unsafe { aoflagger.MakeFlagMask(array_dims.0, array_dims.1, false) };
    let stride = flag_mask.HorizontalStride();
    let flag_buf = flag_mask.pin_mut().BufferMut();

    // TODO: assign by slice
    // flag_buf.copy_from_slice(baseline_flag_view.as_slice().unwrap());
    for (img_timestep_idx, timestep_flag_view) in baseline_flag_view.outer_iter().enumerate() {
        for (img_chan_idx, singular_flag_view) in timestep_flag_view.outer_iter().enumerate() {
            flag_buf[img_chan_idx * stride + img_timestep_idx] =
                *singular_flag_view.get(()).unwrap();
        }
    }
    Ok(flag_mask)
}
