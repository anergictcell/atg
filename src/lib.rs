#![doc = include_str!("../README.md")]
#![doc(
    html_logo_url = "https://raw.githubusercontent.com/anergictcell/atg/main/assets/logo_standard.png",
    html_favicon_url = "https://raw.githubusercontent.com/anergictcell/atg/main/assets/favicon.ico"
)]

pub mod bed;
pub mod fasta;
pub mod genepred;
pub mod genepredext;
pub mod gtf;
pub mod models;
pub mod refgene;
pub mod tests;
pub mod utils;

use crate::models::TranscriptRead;
use crate::models::Transcripts;
use crate::utils::errors::ReadWriteError;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

pub fn read_transcripts<R: TranscriptRead>(
    reader: Result<R, ReadWriteError>,
) -> Result<Transcripts, ReadWriteError> {
    match reader {
        Ok(mut r) => r.transcripts(),
        Err(err) => Err(err),
    }
}
