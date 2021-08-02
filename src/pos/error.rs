
use thiserror::Error;

#[derive(Error, Debug)]
#[error(
    "{source_file}:{source_line} Call to ERFA function {function} returned status code {status}"
)]

/// An error associated with ERFA.
pub struct ErfaError {
    pub source_file: &'static str,
    pub source_line: u32,
    pub status: i32,
    pub function: &'static str,
}