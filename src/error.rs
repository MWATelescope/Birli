//! Errors that can occur in Birli

use marlu::mwalib;
use thiserror::Error;

/// Errors relating to CI
#[derive(Error, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub enum CLIError {
    #[error("Invalid Command Line Argument")]
    /// When a bad CLI argument is provided
    InvalidCommandLineArgument {
        /// The option for which the argument was provided
        option: String,
        /// The argument that was expected
        expected: String,
        /// The argument that was received instead
        received: String,
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

    #[error(transparent)]
    /// Error derived from [`clap::Error`]
    ClapError(#[from] clap::Error),

    #[error(transparent)]
    /// Error derived from [`crate::error::CLIError`]
    CLIError(#[from] CLIError),

    #[error(transparent)]
    /// Error derived from [`marlu::mwalib::MwalibError`]
    MwalibError(#[from] mwalib::MwalibError),

    #[error(transparent)]
    /// Error derived from [`marlu::selection::SelectionError`]
    SelectionError(#[from] marlu::selection::SelectionError),

    #[error("You selected dry run")]
    /// enum variant for when a dry run is selected
    DryRun {},

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

    #[error("Insufficient memory available; need {need_gib} GiB of memory.\nPlease specify the maximum amount of memory to use.")]
    /// Error when we asked for too much memory
    InsufficientMemory {
        /// The amount of memory we think we need
        need_gib: usize,
    },

    #[error("Invalid MWA Version")]
    /// When a bad MWA Version is provided
    BadMWAVersion {
        /// The message to display
        message: String,
        /// The version that was provided
        version: String,
    },
}
