// use fitsio;
use thiserror::Error;

// #[derive(Error, Debug)]
// pub enum FitsError {
//     #[error(
//         "{source_file}:{source_line}\nThere was a fitsio error for {fits_filename}: {fits_error}"
//     )]
//     FitsIO {
//         fits_error: fitsio::errors::Error,
//         fits_filename: String,
//         source_file: String,
//         source_line: u32,
//     },
// }

#[derive(Error, Debug)]
pub enum BirliError {
    /// An error derived from `FitsError`.
    // #[error("{0}")]
    // Fits(#[from] FitsError),
    #[error("{source_file}:{source_line}\nInvalid flag filename template. Must contain \"%%\" (or \"%%%\") for MWAX")]
    InvalidFlagFilenameTemplateError {
        source_file: &'static str,
        source_line: u32,
        filename_template: String,
    },
}
