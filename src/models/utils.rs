use std::fmt;
use std::str::FromStr;

use log::debug;
use serde::{Deserialize, Serialize};

use crate::models::{Transcript, Transcripts};
use crate::utils::errors::ReadWriteError;

/// The Status of the start of stop codon of a CDS
///
/// This is used by RefGene and inferred in GTF based
/// on the absence or presence of dedicated Start- & Stop codon records.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum CdsStat {
    None,
    Unknown,
    Incomplete,
    Complete,
}

impl fmt::Display for CdsStat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                CdsStat::None => "none",
                CdsStat::Unknown => "unk",
                CdsStat::Incomplete => "incmpl",
                CdsStat::Complete => "cmpl",
            }
        )
    }
}

impl FromStr for CdsStat {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "incmpl" => Ok(CdsStat::Incomplete),
            "incompl" => Ok(CdsStat::Incomplete),
            "incomplete" => Ok(CdsStat::Incomplete),
            "cmpl" => Ok(CdsStat::Complete),
            "compl" => Ok(CdsStat::Complete),
            "complete" => Ok(CdsStat::Complete),
            "none" => Ok(CdsStat::None),
            "unk" => Ok(CdsStat::Unknown),
            _ => Err(format!("Invaid CdsStat {}.", s)),
        }
    }
}

/// Indicates the strandness of the transcript
///
/// # Examples
/// ```
/// use std::str::FromStr;
/// use atg::models::Strand;
///
/// let strand = Strand::from_str("+").unwrap();
/// assert_eq!(strand, Strand::Plus);
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Strand::Plus => "+",
                Strand::Minus => "-",
                Strand::Unknown => ".",
            }
        )
    }
}

impl Strand {
    pub fn from_string(s: &str) -> Result<Self, String> {
        Strand::from_str(s)
    }
}

impl FromStr for Strand {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Plus),
            "-" => Ok(Strand::Minus),
            "." => Ok(Strand::Unknown),
            _ => Err(format!(
                "invalid strand {}. Strand must be either `+`, `-` or `.`.",
                s
            )),
        }
    }
}

/// Trait to write ['Transcripts'] or ['Transcript'] into a `Writer` instance
pub trait TranscriptWrite {
    fn writeln_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error>;

    fn write_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error>;

    fn write_transcript_vec(&mut self, transcripts: &[Transcript]) -> Result<(), std::io::Error> {
        debug!("Writing transcripts");
        for transcript in transcripts {
            self.writeln_single_transcript(transcript)?;
        }
        debug!("Finished writing transcripts");
        Ok(())
    }

    fn write_transcripts(&mut self, transcripts: &Transcripts) -> Result<(), std::io::Error> {
        self.write_transcript_vec(transcripts.as_vec())
    }
}

/// Trait for parsing the input into [`Transcripts`]
pub trait TranscriptRead {
    /// Consumes the `Reader` and returns `Transcripts`
    fn transcripts(&mut self) -> Result<Transcripts, ReadWriteError>;
}
