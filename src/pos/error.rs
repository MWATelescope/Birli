//! Error which can occus within the `pos` (positional astronomy) module.

use thiserror::Error;

#[derive(Error, Debug)]
#[error(
    "{source_file}:{source_line} Call to ERFA function {function} returned status code {status}"
)]

/// An error associated with ERFA.
pub struct ErfaError {
    /// The file where the error originated (usually `file!()`)
    pub source_file: &'static str,
    /// The line number where the error originated (usually `line!()`)
    pub source_line: u32,
    /// The status code set by the erfa function
    pub status: i32,
    /// The name of the erfa function called.
    pub function: &'static str,
}
