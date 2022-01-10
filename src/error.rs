//! Errors that can occur in Birli

use thiserror::Error;

#[derive(Error, Debug)]
#[allow(clippy::upper_case_acronyms)]
/// An enum of all the errors possible in Birli
pub enum BirliError {
    #[error("{0}")]
    /// Error derived from [`crate::io::error::IOError`]
    IOError(#[from] crate::io::error::IOError),

    #[error("No common timesteps found. CorrelatorContext hdu info: {hdu_info}")]
    /// Error for when gpuboxes provided have no overlapping visibilities
    NoCommonTimesteps {
        /// display of mwalib::CorrelatorContext::gpubox_time_map
        hdu_info: String,
    },

    #[error("No timesteps were provided. CorrelatorContext hdu info: {hdu_info}")]
    /// Error for when gpuboxes provided have no overlapping visibilities
    NoProvidedTimesteps {
        /// display of mwalib::CorrelatorContext::gpubox_time_map
        hdu_info: String,
    },

    #[error("No common coarse channels found. CorrelatorContext hdu info: {hdu_info}")]
    /// Error for when gpuboxes provided have no overlapping visibilities
    NoCommonCoarseChans {
        /// display of mwalib::CorrelatorContext::gpubox_time_map
        hdu_info: String,
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
