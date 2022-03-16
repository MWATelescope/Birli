//! Selecting a subset of visibilities from an observation using ranges and vectors of indices.
//!
//! Observations can sometimes be too large to fit in memory. This method will only load
//! visibilities from the selected timesteps, coarse channels and baselines in order to enable
//! processing in "chunks"
//!
//! The timesteps are specified as a range of indices in the [`marlu::mwalib::CorrelatorContext`]'s
//! timestep array, which should be a contiguous superset of times from all provided coarse gpubox
//! files. A similar concept applies to coarse channels. Instead of reading visibilities for all
//! known timesteps / coarse channels, it is recommended to use `common_coarse_chan_indices` and
//! `common_timestep_indices`, as these ignore timesteps and coarse channels which are missing
//! contiguous data. `common_good_timestep_indices` is also a good choice to avoid quack time.
//!
//! For more details, see the [documentation](https://docs.rs/mwalib/latest/mwalib/struct.CorrelatorContext.html).
//!
//! Note: it doesn't make sense to ask aoflagger to flag non-contiguous timesteps
//! or coarse channels, and so this interface only allows to ranges to be used.
//! For flagging an obeservation with "picket fence" coarse channels or timesteps,
//! contiguous ranges should be flagged separately.
//!
//! [`marlu::mwalib::CorrelatorContext`]: https://docs.rs/mwalib/latest/mwalib/struct.CorrelatorContext.html
//!
//! # Examples
//!
//! ```rust
//! use birli::{mwalib::CorrelatorContext, FlagContext, VisSelection};
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
//! // Create an mwalib::CorrelatorContext for accessing visibilities.
//! let corr_ctx = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
//!
//! // Determine which timesteps and coarse channels we want to use
//! let img_timestep_idxs = &corr_ctx.common_timestep_indices;
//! let good_timestep_idxs = &corr_ctx.common_good_timestep_indices;
//!
//! let mut vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
//! vis_sel.timestep_range =
//!     *img_timestep_idxs.first().unwrap()..(*img_timestep_idxs.last().unwrap() + 1);
//!
//! // Create a blank array to store flags and visibilities
//! let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
//! let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
//! let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
//!
//! // read visibilities out of the gpubox files
//! vis_sel
//!     .read_mwalib(&corr_ctx, jones_array.view_mut(), flag_array.view_mut(), false)
//!     .unwrap();
//!
//! let dims_common = jones_array.dim();
//!
//! // now try only with good timesteps
//! vis_sel.timestep_range =
//!     *good_timestep_idxs.first().unwrap()..(*good_timestep_idxs.last().unwrap() + 1);
//!
//! // read visibilities out of the gpubox files
//! let mut flag_array = vis_sel.allocate_flags(fine_chans_per_coarse).unwrap();
//! let mut jones_array = vis_sel.allocate_jones(fine_chans_per_coarse).unwrap();
//! vis_sel
//!     .read_mwalib(&corr_ctx, jones_array.view_mut(), flag_array.view_mut(), false)
//!     .unwrap();
//!
//! let dims_good = jones_array.dim();
//!
//! // different selections have different sized arrays.
//! assert_ne!(dims_common, dims_good);
//! ```

use std::ops::Range;

use crossbeam_channel::unbounded;
use crossbeam_utils::thread;
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};
use itertools::izip;
use log::warn;

use crate::{
    marlu::{
        mwalib::{CorrelatorContext, MetafitsContext},
        ndarray::{Array3, ArrayViewMut3, Axis},
        num_traits::Zero,
        rayon::iter::{IndexedParallelIterator, IntoParallelIterator, ParallelIterator},
        Complex, Jones,
    },
    BirliError,
    BirliError::{NoCommonTimesteps, NoProvidedTimesteps},
};

/// Keep track of which mwalib indices the values in a jones array, its' weights and its' flags
/// came from
///
/// TODO: what about <https://doc.rust-lang.org/std/ops/trait.RangeBounds.html> insetad of Range?
#[derive(Debug, Default, Clone)]
pub struct VisSelection {
    /// selected range of mwalib timestep indices
    pub timestep_range: Range<usize>,
    /// selected range of mwalib coarse channel indices
    pub coarse_chan_range: Range<usize>,
    /// selected mwalib baseline indices
    pub baseline_idxs: Vec<usize>,
}

impl VisSelection {
    /// Produce a [`VisSelection`] from a given [`marlu::mwalib::CorrelatorContext`].
    ///
    /// - timesteps are selected from the first [common](https://docs.rs/mwalib/latest/mwalib/struct.CorrelatorContext.html#structfield.common_timestep_indices) to the last [provided](https://docs.rs/mwalib/latest/mwalib/struct.CorrelatorContext.html#structfield.provided_timestep_indices).
    /// - only [common](https://docs.rs/mwalib/latest/mwalib/struct.CorrelatorContext.html#structfield.common_coarse_chan_indices) coarse channels are selected
    /// - all baselines are selected
    ///
    /// For more details, see [mwalib core concepts](https://github.com/MWATelescope/mwalib/wiki/Key-Concepts#timesteps-and-coarse-channels)
    ///
    /// # Examples
    ///
    /// ```rust
    /// use birli::{mwalib::CorrelatorContext, VisSelection};
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
    /// let corr_ctx = CorrelatorContext::new(&metafits_path, &gpufits_paths).unwrap();
    /// let vis_sel = VisSelection::from_mwalib(&corr_ctx).unwrap();
    ///
    /// assert_eq!(vis_sel.timestep_range.len(), 4);
    /// ```
    ///
    /// # Errors
    /// This will return [`BirliError::NoProvidedTimesteps`] if mwalib has determined that no
    /// timesteps have been provided, [`BirliError::NoCommonTimesteps`] if none of the provided
    /// gpuboxes have timesteps in common
    ///
    /// [`marlu::mwalib::CorrelatorContext.common_timestep_indices`]:
    /// [`marlu::mwalib::CorrelatorContext.provided_timestep_indices`]:
    pub fn from_mwalib(corr_ctx: &CorrelatorContext) -> Result<Self, BirliError> {
        Ok(Self {
            timestep_range: match (
                corr_ctx.common_timestep_indices.first(),
                corr_ctx.provided_timestep_indices.last(),
            ) {
                (Some(&first), Some(&last)) if first <= last => (first..last + 1),
                (.., None) => {
                    return Err(NoProvidedTimesteps {
                        hdu_info: format!("{:?}", &corr_ctx.gpubox_time_map),
                    })
                }
                _ => {
                    return Err(NoCommonTimesteps {
                        hdu_info: format!("{:?}", &corr_ctx.gpubox_time_map),
                    })
                }
            },
            coarse_chan_range: match (
                corr_ctx.common_coarse_chan_indices.first(),
                corr_ctx.common_coarse_chan_indices.last(),
            ) {
                (Some(&first), Some(&last)) if first <= last => (first)..(last + 1),
                _ => {
                    return Err(NoCommonTimesteps {
                        hdu_info: format!("{:?}", &corr_ctx.gpubox_time_map),
                    })
                }
            },
            baseline_idxs: (0..corr_ctx.metafits_context.num_baselines).collect(),
            // ..Default::default()
        })

        // let mut result = Self::default();

        // result.timestep_range = match (
        //     corr_ctx.common_timestep_indices.first(),
        //     corr_ctx.provided_timestep_indices.last(),
        // ) {
        //     (Some(&first), Some(&last)) if first <= last => (first..last + 1),
        //     (.., None) => return Err(NoProvidedTimesteps {
        //         hdu_info: format!("{:?}", &corr_ctx.gpubox_time_map),
        //     }),
        //     _ => return Err(NoCommonTimesteps {
        //         hdu_info: format!("{:?}", &corr_ctx.gpubox_time_map),
        //     }),
        // };

        // result.coarse_chan_range = match (corr_ctx.common_coarse_chan_indices.first(), corr_ctx.common_coarse_chan_indices.last()) {
        //     (Some(&first), Some(&last)) if first <= last => (first)..(last + 1),
        //     _ => return Err(NoCommonTimesteps {
        //         hdu_info: format!("{:?}", &corr_ctx.gpubox_time_map),
        //     }),
        // };

        // Ok(result)
    }

    /// The selected antenna index pairs corresponding to `sel_baselines_idxs`
    pub fn get_ant_pairs(&self, meta_ctx: &MetafitsContext) -> Vec<(usize, usize)> {
        self.baseline_idxs
            .iter()
            .map(|&idx| {
                (
                    meta_ctx.baselines[idx].ant1_index,
                    meta_ctx.baselines[idx].ant2_index,
                )
            })
            .collect()
    }

    /// Get the shape of the jones, flag or weight array for this selection
    pub fn get_shape(&self, fine_chans_per_coarse: usize) -> (usize, usize, usize) {
        let num_chans = self.coarse_chan_range.len() * fine_chans_per_coarse;
        let num_baselines = self.baseline_idxs.len();
        let num_timesteps = self.timestep_range.len();
        (num_timesteps, num_chans, num_baselines)
    }

    /// Estimate the memory size in bytes required to store the given selection with redundant pols.
    ///
    /// TODO: rip out everything that does redundant pols
    pub fn estimate_bytes_worst(&self, fine_chans_per_coarse: usize, num_pols: usize) -> usize {
        let shape = self.get_shape(fine_chans_per_coarse);
        shape.0
            * shape.1
            * shape.2
            * num_pols
            * (std::mem::size_of::<Complex<f32>>()
                + std::mem::size_of::<f32>()
                + std::mem::size_of::<bool>())
    }

    /// Estimate the memory size in bytes required to store the given selection without redundant pols.
    pub fn estimate_bytes_best(&self, fine_chans_per_coarse: usize) -> usize {
        let shape = self.get_shape(fine_chans_per_coarse);
        shape.0
            * shape.1
            * shape.2
            * (std::mem::size_of::<Jones<f32>>()
                + std::mem::size_of::<f32>()
                + std::mem::size_of::<bool>())
    }

    /// Allocate a jones array to store visibilities for the selection
    ///
    /// # Errors
    ///
    /// can raise `BirliError::InsufficientMemory` if not enough memory.
    pub fn allocate_jones(
        &self,
        fine_chans_per_coarse: usize,
    ) -> Result<Array3<Jones<f32>>, BirliError> {
        let shape = self.get_shape(fine_chans_per_coarse);
        let num_elems = shape.0 * shape.1 * shape.2;
        let mut v = Vec::new();

        if v.try_reserve_exact(num_elems) == Ok(()) {
            // Make the vector's length equal to its new capacity.
            v.resize(num_elems, Jones::zero());
            Ok(Array3::from_shape_vec(shape, v).unwrap())
        } else {
            // Instead of erroring out with how many GiB we need for *this*
            // array, error out with how many we need for the whole selection.
            let need_gib = self.estimate_bytes_best(fine_chans_per_coarse) / 1024_usize.pow(3);
            Err(BirliError::InsufficientMemory { need_gib })
        }
    }

    /// Allocate a flag array to store flags for the selection
    ///
    /// # Errors
    ///
    /// can raise `BirliError::InsufficientMemory` if not enough memory.
    pub fn allocate_flags(&self, fine_chans_per_coarse: usize) -> Result<Array3<bool>, BirliError> {
        let shape = self.get_shape(fine_chans_per_coarse);
        let num_elems = shape.0 * shape.1 * shape.2;
        let mut v = Vec::new();

        if v.try_reserve_exact(num_elems) == Ok(()) {
            // Make the vector's length equal to its new capacity.
            v.resize(num_elems, false);
            Ok(Array3::from_shape_vec(shape, v).unwrap())
        } else {
            // Instead of erroring out with how many GiB we need for *this*
            // array, error out with how many we need for the whole selection.
            let need_gib = self.estimate_bytes_best(fine_chans_per_coarse) / 1024_usize.pow(3);
            Err(BirliError::InsufficientMemory { need_gib })
        }
    }

    /// Allocate a weight array to store weights for the selection
    ///
    /// # Errors
    ///
    /// can raise `BirliError::InsufficientMemory` if not enough memory.
    pub fn allocate_weights(
        &self,
        fine_chans_per_coarse: usize,
    ) -> Result<Array3<f32>, BirliError> {
        let shape = self.get_shape(fine_chans_per_coarse);
        let num_elems = shape.0 * shape.1 * shape.2;
        let mut v = Vec::new();

        if v.try_reserve_exact(num_elems) == Ok(()) {
            // Make the vector's length equal to its new capacity.
            v.resize(num_elems, 0.);
            Ok(Array3::from_shape_vec(shape, v).unwrap())
        } else {
            // Instead of erroring out with how many GiB we need for *this*
            // array, error out with how many we need for the whole selection.
            let need_gib = self.estimate_bytes_best(fine_chans_per_coarse) / 1024_usize.pow(3);
            Err(BirliError::InsufficientMemory { need_gib })
        }
    }

    /// Read the visibilities for this selection into the jones array using mwalib,
    /// flag visiblities if they are not provided.
    ///
    /// # Errors
    ///
    /// Can raise [`BirliError::BadArrayShape`] if `jones_array` or `flag_array` does not match the
    /// expected shape of this selection.
    pub fn read_mwalib(
        &self,
        corr_ctx: &CorrelatorContext,
        mut jones_array: ArrayViewMut3<Jones<f32>>,
        mut flag_array: ArrayViewMut3<bool>,
        draw_progress: bool,
    ) -> Result<(), BirliError> {
        let fine_chans_per_coarse = corr_ctx.metafits_context.num_corr_fine_chans_per_coarse;
        let shape = self.get_shape(fine_chans_per_coarse);
        let (num_timesteps, _, _) = shape;
        let num_coarse_chans = self.coarse_chan_range.len();

        if jones_array.dim() != shape {
            return Err(BirliError::BadArrayShape {
                argument: "jones_array".to_string(),
                function: "VisSelection::read_mwalib".to_string(),
                expected: format!("{:?}", shape),
                received: format!("{:?}", jones_array.dim()),
            });
        };

        if flag_array.dim() != shape {
            return Err(BirliError::BadArrayShape {
                argument: "flag_array".to_string(),
                function: "VisSelection::read_mwalib".to_string(),
                expected: format!("{:?}", shape),
                received: format!("{:?}", flag_array.dim()),
            });
        };

        // since we are using read_by_by basline into buffer, the visibilities are read in order:
        // baseline,frequency,pol,r,i

        let floats_per_chan = corr_ctx.metafits_context.num_visibility_pols * 2;
        let floats_per_baseline = floats_per_chan * fine_chans_per_coarse;
        let floats_per_hdu = floats_per_baseline * corr_ctx.metafits_context.num_baselines;

        // A queue of errors
        let (tx_error, rx_error) = unbounded();

        // Progress bar draw target
        let draw_target = if draw_progress {
            ProgressDrawTarget::stderr()
        } else {
            ProgressDrawTarget::hidden()
        };
        // a progress bar containing the progress bars associated with this method
        let multi_progress = MultiProgress::with_draw_target(draw_target);
        // a vector of progress bars for the visibility reading progress of each channel.
        let read_progress: Vec<ProgressBar> = self
            .coarse_chan_range
            .clone()
            .map(|mwalib_coarse_chan_idx| {
                let channel_progress = multi_progress.add(
                    ProgressBar::new(num_timesteps as _)
                        .with_style(
                            ProgressStyle::default_bar()
                                .template("{msg:16}: [{wide_bar:.blue}] {pos:4}/{len:4}")
                                .progress_chars("=> "),
                        )
                        .with_position(0)
                        .with_message(format!("coarse_chan {:03}", mwalib_coarse_chan_idx)),
                );
                channel_progress.set_position(0);
                channel_progress
            })
            .collect();
        // The total reading progress bar.
        let total_progress = multi_progress.add(
            ProgressBar::new((num_timesteps * num_coarse_chans) as _)
                .with_style(
                    ProgressStyle::default_bar()
                        .template(
                            "{msg:16}: [{elapsed_precise}] [{wide_bar:.cyan/blue}] {percent:3}% ({eta:5})",
                        )
                        .progress_chars("=> "),
                )
                .with_position(0)
                .with_message("loading hdus"),
        );

        thread::scope(|scope| {
            // Spawn a thread to draw the progress bars.
            scope.spawn(|_| {
                multi_progress.join().unwrap();
            });

            total_progress.set_position(0);

            // Error handling thread
            scope.spawn(|_| {
                for (mwalib_timestep_idx, mwalib_coarse_chan_idx, err) in rx_error {
                    warn!(
                        "Flagging missing HDU @ ts={}, cc={} {:?}",
                        mwalib_timestep_idx, mwalib_coarse_chan_idx, err
                    );
                }
            });

            // Load HDUs from each coarse channel. arrays: [timestep][chan][baseline]

            //refactor below
            jones_array
                .axis_chunks_iter_mut(Axis(1), fine_chans_per_coarse)
                .into_par_iter()
                .zip(flag_array.axis_chunks_iter_mut(Axis(1), fine_chans_per_coarse))
                .zip(self.coarse_chan_range.clone())
                .zip(read_progress)
                .for_each(
                    |(((mut jones_array, mut flag_array), coarse_chan_idx), progress)| {
                        progress.set_position(0);

                        // buffer: [baseline][chan][pol][complex]
                        let mut hdu_buffer: Vec<f32> = vec![0.0; floats_per_hdu];

                        // arrays: [chan][baseline]
                        for (mut jones_array, mut flag_array, timestep_idx) in izip!(
                            jones_array.outer_iter_mut(),
                            flag_array.outer_iter_mut(),
                            self.timestep_range.clone(),
                        ) {
                            if let Err(err) = corr_ctx.read_by_baseline_into_buffer(
                                timestep_idx,
                                coarse_chan_idx,
                                hdu_buffer.as_mut_slice(),
                            ) {
                                flag_array.fill(true);
                                tx_error.send((timestep_idx, coarse_chan_idx, err)).unwrap();
                                continue;
                            }
                            // arrays: [chan]
                            for (mut jones_array, baseline_idx) in izip!(
                                jones_array.axis_iter_mut(Axis(1)),
                                self.baseline_idxs.clone()
                            ) {
                                let buf_chunk_start = baseline_idx * floats_per_baseline;
                                // buffer: [chan][pol][complex]
                                let hdu_baseline_chunk = &hdu_buffer
                                    [buf_chunk_start..(buf_chunk_start + floats_per_baseline)];
                                for (jones, hdu_chan_chunk) in izip!(
                                    jones_array.iter_mut(),
                                    hdu_baseline_chunk.chunks(floats_per_chan)
                                ) {
                                    *jones = Jones::from([
                                        Complex::new(hdu_chan_chunk[0], hdu_chan_chunk[1]),
                                        Complex::new(hdu_chan_chunk[2], hdu_chan_chunk[3]),
                                        Complex::new(hdu_chan_chunk[4], hdu_chan_chunk[5]),
                                        Complex::new(hdu_chan_chunk[6], hdu_chan_chunk[7]),
                                    ]);
                                }
                            }

                            progress.inc(1);
                            total_progress.inc(1);
                        }
                        progress.finish();
                    },
                );

            drop(tx_error);

            // We're done!
            total_progress.finish();
        })
        .unwrap();

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_abs_diff_eq;

    use crate::{
        compare_jones,
        marlu::{num_traits::Zero, Complex},
        test_common::{get_mwa_ord_context, get_mwa_ord_dodgy_context, get_mwax_context},
        TestJones, VisSelection,
    };

    use super::*;

    #[test]
    /// We expect coarse channel 0 ( fine channels 0,1 ) to be the same as in `get_mwa_ord_context`,
    /// but coarse channel 0 (fine channels 2, 3 ) should be shifted.
    fn test_read_mwalib_mwax_flags_missing_hdus() {
        let corr_ctx = get_mwa_ord_dodgy_context();
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

        // ts 0, chan 0, baseline 0
        compare_jones!(
            jones_array[(0, 0, 0)],
            Jones::from([
                Complex::new(0x10c5be as f32, -0x10c5bf as f32),
                Complex::new(0x10c5ae as f32, 0x10c5af as f32),
                Complex::new(0x10c5ae as f32, -0x10c5af as f32),
                Complex::new(0x10bec6 as f32, -0x10bec7 as f32),
            ])
        );
        // ts 1, chan 0, baseline 0
        compare_jones!(
            jones_array[(1, 0, 0)],
            Jones::from([
                Complex::new(0x14c5be as f32, -0x14c5bf as f32),
                Complex::new(0x14c5ae as f32, 0x14c5af as f32),
                Complex::new(0x14c5ae as f32, -0x14c5af as f32),
                Complex::new(0x14bec6 as f32, -0x14bec7 as f32),
            ])
        );
        // ts 2, chan 0, baseline 0
        compare_jones!(
            jones_array[(2, 0, 0)],
            Jones::from([
                Complex::new(0x18c5be as f32, -0x18c5bf as f32),
                Complex::new(0x18c5ae as f32, 0x18c5af as f32),
                Complex::new(0x18c5ae as f32, -0x18c5af as f32),
                Complex::new(0x18bec6 as f32, -0x18bec7 as f32),
            ])
        );
        // ts 3, chan 0, baseline 0
        compare_jones!(
            jones_array[(3, 0, 0)],
            Jones::from([
                Complex::new(0x1cc5be as f32, -0x1cc5bf as f32),
                Complex::new(0x1cc5ae as f32, 0x1cc5af as f32),
                Complex::new(0x1cc5ae as f32, -0x1cc5af as f32),
                Complex::new(0x1cbec6 as f32, -0x1cbec7 as f32),
            ])
        );

        // Fine channel 2 is in the modified coarse channel, so should have drifted

        // ts 0, chan 2, baseline 0
        compare_jones!(
            jones_array[(0, 2, 0)],
            Jones::from([
                Complex::new(0x04c5be as f32, -0x04c5bf as f32),
                Complex::new(0x04c5ae as f32, 0x04c5af as f32),
                Complex::new(0x04c5ae as f32, -0x04c5af as f32),
                Complex::new(0x04bec6 as f32, -0x04bec7 as f32),
            ])
        );
        assert!(!flag_array[(0, 2, 0)]);
        // ts 1, chan 2, baseline 0
        compare_jones!(
            jones_array[(1, 2, 0)],
            TestJones::from(Jones::<f32>::zero())
        );
        assert!(flag_array[(1, 2, 0)]);
        // ts 2, chan 2, baseline 0 - unchanged
        compare_jones!(
            jones_array[(2, 2, 0)],
            Jones::from([
                Complex::new(0x08c5be as f32, -0x08c5bf as f32),
                Complex::new(0x08c5ae as f32, 0x08c5af as f32),
                Complex::new(0x08c5ae as f32, -0x08c5af as f32),
                Complex::new(0x08bec6 as f32, -0x08bec7 as f32),
            ])
        );
        assert!(!flag_array[(2, 2, 0)]);
        // ts 3, chan 2, baseline 0
        compare_jones!(
            jones_array[(3, 2, 0)],
            TestJones::from(Jones::<f32>::zero())
        );
        assert!(flag_array[(3, 2, 0)]);
    }

    #[test]
    fn test_read_mwalib_mwax() {
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

        // ts 0, chan 0, baseline 0
        compare_jones!(
            jones_array[(0, 0, 0)],
            Jones::from([
                Complex::new(0x410000 as f32, 0x410001 as f32),
                Complex::new(0x410002 as f32, 0x410003 as f32),
                Complex::new(0x410004 as f32, 0x410005 as f32),
                Complex::new(0x410006 as f32, 0x410007 as f32),
            ])
        );

        // ts 0, chan 0 (cc 0, fc 0), baseline 1
        compare_jones!(
            jones_array[(0, 0, 1)],
            Jones::from([
                Complex::new(0x410010 as f32, 0x410011 as f32),
                Complex::new(0x410012 as f32, 0x410013 as f32),
                Complex::new(0x410014 as f32, 0x410015 as f32),
                Complex::new(0x410016 as f32, 0x410017 as f32),
            ])
        );

        // ts 0, chan 0, baseline 2
        compare_jones!(
            jones_array[(0, 0, 2)],
            Jones::from([
                Complex::new(0x410020 as f32, 0x410021 as f32),
                Complex::new(0x410022 as f32, 0x410023 as f32),
                Complex::new(0x410024 as f32, 0x410025 as f32),
                Complex::new(0x410026 as f32, 0x410027 as f32),
            ])
        );

        // ts 0, chan 1, baseline 2
        compare_jones!(
            jones_array[(0, 1, 2)],
            Jones::from([
                Complex::new(0x410028 as f32, 0x410029 as f32),
                Complex::new(0x41002a as f32, 0x41002b as f32),
                Complex::new(0x41002c as f32, 0x41002d as f32),
                Complex::new(0x41002e as f32, 0x41002f as f32),
            ])
        );

        // ts 1, chan 0, baseline 0
        compare_jones!(
            jones_array[(1, 0, 0)],
            Jones::from([
                Complex::new(0x410100 as f32, 0x410101 as f32),
                Complex::new(0x410102 as f32, 0x410103 as f32),
                Complex::new(0x410104 as f32, 0x410105 as f32),
                Complex::new(0x410106 as f32, 0x410107 as f32),
            ])
        );

        // ts 2, chan 0, baseline 0
        compare_jones!(
            jones_array[(2, 0, 0)],
            Jones::from([
                Complex::new(0x410200 as f32, 0x410201 as f32),
                Complex::new(0x410202 as f32, 0x410203 as f32),
                Complex::new(0x410204 as f32, 0x410205 as f32),
                Complex::new(0x410206 as f32, 0x410207 as f32),
            ])
        );

        // ts 3, chan 3 (cc 1, fc 1), baseline 1
        compare_jones!(
            jones_array[(3, 3, 1)],
            Jones::from([
                Complex::new(0x410718 as f32, 0x410719 as f32),
                Complex::new(0x41071a as f32, 0x41071b as f32),
                Complex::new(0x41071c as f32, 0x41071d as f32),
                Complex::new(0x41071e as f32, 0x41071f as f32),
            ])
        );
    }

    #[test]
    fn test_read_mwalib_mwa_ord() {
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

        // ts 0, chan 0, baseline 0
        compare_jones!(
            jones_array[(0, 0, 0)],
            Jones::from([
                Complex::new(0x10c5be as f32, -0x10c5bf as f32),
                Complex::new(0x10c5ae as f32, 0x10c5af as f32),
                Complex::new(0x10c5ae as f32, -0x10c5af as f32),
                Complex::new(0x10bec6 as f32, -0x10bec7 as f32),
            ])
        );
        // ts 1, chan 0, baseline 0
        compare_jones!(
            jones_array[(1, 0, 0)],
            Jones::from([
                Complex::new(0x14c5be as f32, -0x14c5bf as f32),
                Complex::new(0x14c5ae as f32, 0x14c5af as f32),
                Complex::new(0x14c5ae as f32, -0x14c5af as f32),
                Complex::new(0x14bec6 as f32, -0x14bec7 as f32),
            ])
        );
        // ts 2, chan 0, baseline 0
        compare_jones!(
            jones_array[(2, 0, 0)],
            Jones::from([
                Complex::new(0x18c5be as f32, -0x18c5bf as f32),
                Complex::new(0x18c5ae as f32, 0x18c5af as f32),
                Complex::new(0x18c5ae as f32, -0x18c5af as f32),
                Complex::new(0x18bec6 as f32, -0x18bec7 as f32),
            ])
        );
        // ts 3, chan 0, baseline 0
        compare_jones!(
            jones_array[(3, 0, 0)],
            Jones::from([
                Complex::new(0x1cc5be as f32, -0x1cc5bf as f32),
                Complex::new(0x1cc5ae as f32, 0x1cc5af as f32),
                Complex::new(0x1cc5ae as f32, -0x1cc5af as f32),
                Complex::new(0x1cbec6 as f32, -0x1cbec7 as f32),
            ])
        );

        // ts 0, chan 0, baseline 5
        compare_jones!(
            jones_array[(0, 0, 5)],
            Jones::from([
                Complex::new(0x10f1ce as f32, -0x10f1cf as f32),
                Complex::new(0x10ea26 as f32, -0x10ea27 as f32),
                Complex::new(0x10f1be as f32, -0x10f1bf as f32),
                Complex::new(0x10ea16 as f32, -0x10ea17 as f32),
            ])
        );
        // ts 1, chan 0, baseline 5
        compare_jones!(
            jones_array[(1, 0, 5)],
            Jones::from([
                Complex::new(0x14f1ce as f32, -0x14f1cf as f32),
                Complex::new(0x14ea26 as f32, -0x14ea27 as f32),
                Complex::new(0x14f1be as f32, -0x14f1bf as f32),
                Complex::new(0x14ea16 as f32, -0x14ea17 as f32),
            ])
        );
        // ts 2, chan 0, baseline 5
        compare_jones!(
            jones_array[(2, 0, 5)],
            Jones::from([
                Complex::new(0x18f1ce as f32, -0x18f1cf as f32),
                Complex::new(0x18ea26 as f32, -0x18ea27 as f32),
                Complex::new(0x18f1be as f32, -0x18f1bf as f32),
                Complex::new(0x18ea16 as f32, -0x18ea17 as f32),
            ])
        );
        // ts 3, chan 0, baseline 5
        compare_jones!(
            jones_array[(3, 0, 5)],
            Jones::from([
                Complex::new(0x1cf1ce as f32, -0x1cf1cf as f32),
                Complex::new(0x1cea26 as f32, -0x1cea27 as f32),
                Complex::new(0x1cf1be as f32, -0x1cf1bf as f32),
                Complex::new(0x1cea16 as f32, -0x1cea17 as f32),
            ])
        );

        // ts 0, chan 2, baseline 0
        compare_jones!(
            jones_array[(0, 2, 0)],
            Jones::from([
                Complex::new(0x00c5be as f32, -0x00c5bf as f32),
                Complex::new(0x00c5ae as f32, 0x00c5af as f32),
                Complex::new(0x00c5ae as f32, -0x00c5af as f32),
                Complex::new(0x00bec6 as f32, -0x00bec7 as f32),
            ])
        );
        // ts 1, chan 2, baseline 0
        compare_jones!(
            jones_array[(1, 2, 0)],
            Jones::from([
                Complex::new(0x04c5be as f32, -0x04c5bf as f32),
                Complex::new(0x04c5ae as f32, 0x04c5af as f32),
                Complex::new(0x04c5ae as f32, -0x04c5af as f32),
                Complex::new(0x04bec6 as f32, -0x04bec7 as f32),
            ])
        );
        // ts 2, chan 2, baseline 0
        compare_jones!(
            jones_array[(2, 2, 0)],
            Jones::from([
                Complex::new(0x08c5be as f32, -0x08c5bf as f32),
                Complex::new(0x08c5ae as f32, 0x08c5af as f32),
                Complex::new(0x08c5ae as f32, -0x08c5af as f32),
                Complex::new(0x08bec6 as f32, -0x08bec7 as f32),
            ])
        );
        // ts 3, chan 2, baseline 0
        compare_jones!(
            jones_array[(3, 2, 0)],
            Jones::from([
                Complex::new(0x0cc5be as f32, -0x0cc5bf as f32),
                Complex::new(0x0cc5ae as f32, 0x0cc5af as f32),
                Complex::new(0x0cc5ae as f32, -0x0cc5af as f32),
                Complex::new(0x0cbec6 as f32, -0x0cbec7 as f32),
            ])
        );
        // ts 3, chan 3, baseline 5
        compare_jones!(
            jones_array[(3, 3, 5)],
            Jones::from([
                Complex::new(0x0df3ce as f32, -0x0df3cf as f32),
                Complex::new(0x0dec26 as f32, -0x0dec27 as f32),
                Complex::new(0x0df3be as f32, -0x0df3bf as f32),
                Complex::new(0x0dec16 as f32, -0x0dec17 as f32),
            ])
        );
    }
}
