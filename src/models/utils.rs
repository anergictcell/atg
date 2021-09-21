use std::fmt;
use std::str::FromStr;

use crate::models::{Transcript, Transcripts};
use crate::utils::errors::ReadWriteError;

#[derive(Clone, Copy, Debug, PartialEq)]
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

#[derive(Clone, Copy, Debug, PartialEq)]
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

pub trait TranscriptWrite {
    fn writeln_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error>;

    fn write_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error>;

    fn write_transcript_vec(&mut self, transcripts: &[Transcript]) -> Result<(), std::io::Error> {
        for transcript in transcripts {
            self.writeln_single_transcript(transcript)?;
        }
        Ok(())
    }

    fn write_transcripts(&mut self, transcripts: &Transcripts) -> Result<(), std::io::Error> {
        self.write_transcript_vec(transcripts.as_vec())
    }
}

pub trait TranscriptRead {
    fn transcripts(&mut self) -> Result<Transcripts, ReadWriteError>;
}

/// Ensures that the chrom values internally will always contain a "chr" prefix
///
/// # Examples
///
/// ```rust
/// use atg::models::parse_chrom;
///
/// assert_eq!(parse_chrom("1"), "chr1");
/// assert_eq!(parse_chrom("chr1"), "chr1");
/// assert_eq!(parse_chrom("CHR1"), "chr1");
/// assert_eq!(parse_chrom("CHr1"), "chr1");
/// assert_eq!(parse_chrom("chrM"), "chrM");
/// assert_eq!(parse_chrom("M"), "chrM");
/// assert_eq!(parse_chrom("chrMT"), "chrMT");
/// assert_eq!(parse_chrom("MT"), "chrMT");
/// ```
pub fn parse_chrom(s: &str) -> String {
    let res = String::from("chr");
    if s.to_ascii_lowercase().starts_with("chr") {
        res + &s[3..]
    } else {
        res + s
    }
}
