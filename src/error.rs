use thiserror::Error;

#[derive(Error, Debug)]
pub enum BirliError {
    /// An error derived from `FitsError`.
    #[error("{source_file}:{source_line}\nInvalid flag filename template. Must contain \"%%\" (or \"%%%\") for MWAX")]
    InvalidFlagFilenameTemplate {
        source_file: &'static str,
        source_line: u32,
        filename_template: String,
    },
    /// Error when opening a fits file.
    #[error("{source_file}:{source_line}\nCouldn't open {fits_filename}: {fits_error}")]
    FitsOpen {
        fits_error: fitsio::errors::Error,
        fits_filename: String,
        source_file: &'static str,
        source_line: u32,
    },
    /// A generic error associated with the fitsio crate.
    #[error("{source_file}:{source_line}\n{fits_filename} HDU {hdu_num}: {fits_error}")]
    #[allow(clippy::upper_case_acronyms)]
    FitsIO {
        fits_error: fitsio::errors::Error,
        fits_filename: String,
        hdu_num: usize,
        source_file: &'static str,
        source_line: u32,
    },

    #[error("{0}")]
    FitsError(#[from] mwalib::FitsError),

    #[error("{0}")]
    FitsioError(#[from] fitsio::errors::Error),

    /// Error to describe some kind of inconsistent state within an mwaf file.
    #[error("Inconsistent mwaf file (file: {file}, expected: {expected}, found: {found})")]
    MwafInconsistent {
        file: String,
        expected: String,
        found: String,
    },

    #[error("Invalid GPUBox ID {found}, expected on of {expected}")]
    InvalidGpuBox { expected: String, found: String },
}
