//! Errors that can occur in the io module

use thiserror::Error;

use marlu::{fitsio, mwalib, mwalib::FitsError};

#[derive(Error, Debug)]
#[allow(clippy::upper_case_acronyms)]
/// All the errors that can occur in file io operations
pub enum IOError {
    /// An error derived from `FitsError`.
    #[error("{source_file}:{source_line}\nInvalid flag filename template. Must contain \"%%\" (or \"%%%\") for MWAX")]
    InvalidFlagFilenameTemplate {
        /// The file where the error originated (usually `file!()`)
        source_file: &'static str,
        /// The line number where the error originated (usually `line!()`)
        source_line: u32,
        /// The filename templte
        filename_template: String,
    },
    /// Error when opening a fits file.
    #[error("{source_file}:{source_line}\nCouldn't open {fits_filename}: {fits_error}")]
    FitsOpen {
        /// The [`fitsio::errors::Error`]
        fits_error: fitsio::errors::Error,
        /// The filename of the fits file
        fits_filename: String,
        /// The file where the error originated (usually `file!()`)
        source_file: &'static str,
        /// The line number where the error originated (usually `line!()`)
        source_line: u32,
    },
    /// A generic error associated with the fitsio crate.
    #[error("{source_file}:{source_line}\n{fits_filename} HDU {hdu_num}: {fits_error}")]
    #[allow(clippy::upper_case_acronyms)]
    FitsIO {
        /// The [`fitsio::errors::Error`]
        fits_error: fitsio::errors::Error,
        /// The filename of the fits file where the error occurred
        fits_filename: String,
        /// The hdu number in the fits file where the error occurred
        hdu_num: usize,
        /// The file where the error originated (usually `file!()`)
        source_file: &'static str,
        /// The line number where the error originated (usually `line!()`)
        source_line: u32,
    },

    #[error(transparent)]
    /// Error derived from [`marlu::mwalib::FitsError`]
    FitsError(#[from] mwalib::FitsError),

    #[error(transparent)]
    /// Error derived from [`fitsio::errors::Error`]
    FitsioError(#[from] fitsio::errors::Error),

    #[error(transparent)]
    /// Error derived from [`marlu::io::error::IOError`]
    MarluIOError(#[from] marlu::io::error::IOError),

    #[error(transparent)]
    /// Error derived from [`marlu::io::error::UvfitsWriteError`]
    UvfitsWriteError(#[from] marlu::io::error::UvfitsWriteError),

    /// Error to describe some kind of inconsistent state within an mwaf file.
    #[error("Inconsistent mwaf file (file: {file}, expected: {expected}, found: {found})")]
    MwafInconsistent {
        /// The filename of the fits file where the error occurred
        file: String,
        /// The value that was expected
        expected: String,
        /// The unexpected value that was found
        found: String,
    },

    #[error("Invalid GPUBox ID {found}, expected on of {expected}")]
    /// Error for an unexpected gpubox ID
    InvalidGpuBox {
        /// The value that was expected
        expected: String,
        /// The unexpected value that was found
        found: String,
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

#[derive(Error, Debug)]
/// Errors that can occur when reading from a solution file
pub enum ReadSolutionsError {
    #[error("Tried to read calibration solutions file with an unsupported extension '{ext}'!")]
    #[allow(missing_docs)]
    UnsupportedExt { ext: String },

    #[error(
        "When reading {file}, expected MWAOCAL as the first 7 characters, got '{got}' instead!"
    )]
    #[allow(missing_docs)]
    AndreBinaryStr { file: String, got: String },

    #[error(
        "When reading {file}, expected a value {expected} in the header, but got '{got}' instead!"
    )]
    #[allow(missing_docs)]
    AndreBinaryVal {
        file: String,
        expected: &'static str,
        got: String,
    },

    #[error("In file {file} key {key}, could not parse '{got}' as a number!")]
    #[allow(missing_docs)]
    Parse {
        file: String,
        key: &'static str,
        got: String,
    },

    #[error(transparent)]
    #[allow(missing_docs)]
    Fits(#[from] FitsError),

    #[error("IO error: {0}")]
    #[allow(missing_docs)]
    IO(#[from] std::io::Error),
}
