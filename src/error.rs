//! Errors that can occur in Birli

use marlu::{io::error::BadArrayShape, mwalib};
#[cfg(feature = "cli")]
use shlex::QuoteError;
use thiserror::Error;

use crate::{
    corrections::{DigitalGainCorrection, PassbandCorrection},
    van_vleck::VanVleckCorrection,
};

/// Errors relating to CI
#[cfg(feature = "cli")]
#[derive(Error, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub enum CLIError {
    #[error("Invalid Command Line Argument: option {option} expected {expected} but received {received}")]
    /// When a bad CLI argument is provided
    InvalidCommandLineArgument {
        /// The option for which the argument was provided
        option: String,
        /// The argument that was expected
        expected: String,
        /// The argument that was received instead
        received: String,
    },
    #[error("Invalid range specifier: {reason}")]
    /// When a bad range specifier is provided
    InvalidRangeSpecifier {
        /// Why the range is invalid
        reason: String,
    },
}

/// An enum of all the errors possible in Birli
#[derive(Error, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub enum BirliError {
    #[error(transparent)]
    /// Error derived from [`crate::io::error::IOError`]
    IOError(#[from] crate::io::error::IOError),

    #[error(transparent)]
    /// Error derived from [`crate::calibration::CalibrationError`]
    CalibrationError(#[from] crate::calibration::CalibrationError),

    #[cfg(feature = "cli")]
    #[error(transparent)]
    /// Error derived from [`QuoteError`]
    QuoteError(#[from] QuoteError),

    #[cfg(feature = "cli")]
    #[error(transparent)]
    /// Error derived from [`clap::Error`]
    ClapError(#[from] clap::Error),

    #[cfg(feature = "cli")]
    #[error(transparent)]
    /// Error derived from `CLIError`
    CLIError(#[from] CLIError),

    #[error(transparent)]
    /// Error derived from [`marlu::mwalib::MwalibError`]
    MwalibError(#[from] mwalib::MwalibError),

    #[error(transparent)]
    /// Error derived from [`crate::marlu::selection::SelectionError`]
    SelectionError(#[from] marlu::selection::SelectionError),

    #[error(transparent)]
    /// Error derived from [`crate::corrections::PassbandCorrection`]
    PassbandCorrection(#[from] PassbandCorrection),

    #[error(transparent)]
    /// Error derived from [`crate::corrections::DigitalGainCorrection`]
    DigitalGainCorrection(#[from] DigitalGainCorrection),

    #[error(transparent)]
    /// Error derived from [`crate::van_vleck::VanVleckCorrection`]
    VanVleckCorrection(#[from] VanVleckCorrection),

    #[error("You selected dry run")]
    /// enum variant for when a dry run is selected
    DryRun {},

    #[error(transparent)]
    /// Error for bad array shape in provided argument
    BadArrayShape(#[from] BadArrayShape),

    #[error("Insufficient memory available; need {need_gib} GiB of memory.\nPlease specify the maximum amount of memory to use.")]
    /// Error when we asked for too much memory
    InsufficientMemory {
        /// The amount of memory we think we need
        need_gib: usize,
    },

    #[error("Invalid MWA Version ({version}): {message}")]
    /// When a bad MWA Version is provided
    BadMWAVersion {
        /// The message to display
        message: String,
        /// The version that was provided
        version: String,
    },
}
