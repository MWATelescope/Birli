//! Errors that can occur in the io module

use thiserror::Error;

use mwa_rust_core::{fitsio, mwalib};

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
    // TODO: address this
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

    #[error("{0}")]
    /// Error derived from [`mwalib::FitsError`]
    FitsError(#[from] mwalib::FitsError),

    #[error("{0}")]
    /// Error derived from [`fitsio::errors::Error`]
    FitsioError(#[from] fitsio::errors::Error),

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
}

///
/// Most of this was blatently stolen (with permission) from [Chris Jordan](https://github.com/cjordan)
#[derive(Error, Debug)]
pub enum UvfitsWriteError {
    /// An error when trying to write to an unexpected row.
    #[error("Tried to write to row number {row_num}, but only {num_rows} rows are expected")]
    BadRowNum {
        /// The row number (0-indexed)
        row_num: usize,
        /// Total number of rows expected.
        num_rows: usize,
    },

    /// An error when less rows were written to an HDU than expected.
    #[error("Expected {total} uvfits rows to be written, but only {current} were written")]
    NotEnoughRowsWritten {
        /// Number of rows written
        current: usize,
        /// Total number of rows expected.
        total: usize,
    },

    /// An error associated with ERFA.
    #[error("{0}")]
    Erfa(#[from] mwa_rust_core::pos::ErfaError),

    /// An error associated with fitsio.
    #[error("{0}")]
    Fitsio(#[from] fitsio::errors::Error),

    /// An error when converting a Rust string to a C string.
    #[error("{0}")]
    BadString(#[from] std::ffi::NulError),

    /// An IO error.
    #[error("{0}")]
    IO(#[from] std::io::Error),
}
