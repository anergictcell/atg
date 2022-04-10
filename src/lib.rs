#![doc = include_str!("../README.md")]

pub mod bed;
pub mod fasta;
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
