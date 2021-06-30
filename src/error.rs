//! Errors that can occur in Birli

use thiserror::Error;

#[derive(Error, Debug)]
#[allow(clippy::upper_case_acronyms)]
/// An enum of all the errors possible in Birli
pub enum BirliError {
    #[error("{0}")]
    /// Error derived from [`crate::io::error::IOError`]
    IOError(#[from] crate::io::error::IOError),

    #[error("No common timesteps found. CorrelatorContext timestep info: {timestep_info}")]
    /// Error for when gpuboxes provided have no overlapping visibilities
    NoCommonTimesteps {
        /// display of mwalib::CorrelatorContext::gpubox_time_map
        timestep_info: String,
    },

    #[error("No timesteps were provided. CorrelatorContext timestep info: {timestep_info}")]
    /// Error for when gpuboxes provided have no overlapping visibilities
    NoProvidedTimesteps {
        /// display of mwalib::CorrelatorContext::gpubox_time_map
        timestep_info: String,
    },
}
