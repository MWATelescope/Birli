//! Selecting a subset of visibilities from an observation using ranges and vectors of indices.

use std::ops::Range;

use marlu::{
    mwalib::{CorrelatorContext, MetafitsContext},
    Complex,
};

use crate::{
    BirliError,
    BirliError::{NoCommonTimesteps, NoProvidedTimesteps},
};

/// Keep track of which mwalib indices the values in a jones array, its' weights and its' flags
/// came from
#[derive(Debug, Default, Clone)]
pub struct VisSelection {
    // TODO: remove sel_ prefix
    /// selected range of mwalib timestep indices
    pub timestep_range: Range<usize>,
    /// selected range of mwalib coarse channel indices
    pub coarse_chan_range: Range<usize>,
    /// selected mwalib baseline indices
    pub baseline_idxs: Vec<usize>,
}

impl VisSelection {
    /// Produce a VisSelection from a given [`mwalib::CorrelatorContext`].
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
    /// [`mwalib::CorrelatorContext.common_timestep_indices`]:
    /// [`mwalib::CorrelatorContext.provided_timestep_indices`]:
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

    /// Estimate the memory size in bytes required to store the given selection.
    pub fn estimate_bytes(&self, fine_chans_per_coarse: usize, num_pols: usize) -> usize {
        let num_sel_chans = self.coarse_chan_range.len() * fine_chans_per_coarse;
        let num_sel_baselines = self.baseline_idxs.len();
        let num_sel_timesteps = self.timestep_range.len();
        num_sel_chans
            * num_sel_baselines
            * num_sel_timesteps
            * num_pols
            * (std::mem::size_of::<Complex<f32>>()
                + std::mem::size_of::<f32>()
                + std::mem::size_of::<bool>())
    }
}
