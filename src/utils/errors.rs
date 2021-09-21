use crate::gtf::ParseGtfError;
use crate::refgene::ParseRefGeneError;
use std::error::Error;
use std::fmt;

pub struct MissingCDSError;

impl Error for MissingCDSError {}

impl fmt::Display for MissingCDSError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "the exon does not have a coding sequence!")
    }
}

impl fmt::Debug for MissingCDSError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{{ file: {}, line: {} }}", file!(), line!())
    }
}

#[derive(Debug)]
pub struct BuildTranscriptError {
    message: String,
}

impl Error for BuildTranscriptError {}

impl BuildTranscriptError {
    pub fn new<S: fmt::Display>(s: S) -> Self {
        Self {
            message: s.to_string(),
        }
    }
}

impl fmt::Display for BuildTranscriptError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Unable to build the transcript")
    }
}

#[derive(Debug)]
pub struct BuildCodonError {
    message: String,
}

impl Error for BuildCodonError {}

impl BuildCodonError {
    pub fn new<S: fmt::Display>(s: S) -> Self {
        Self {
            message: s.to_string(),
        }
    }
}

impl fmt::Display for BuildCodonError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "unable to build a Codon: {}", self.message)
    }
}

#[derive(Debug)]
pub struct ReadWriteError {
    message: String,
}

impl Error for ReadWriteError {}

impl ReadWriteError {
    pub fn new<S: fmt::Display>(s: S) -> Self {
        Self {
            message: s.to_string(),
        }
    }
}

impl fmt::Display for ReadWriteError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Unable to read or write: {}", self.message)
    }
}

impl From<ParseGtfError> for ReadWriteError {
    fn from(err: ParseGtfError) -> Self {
        Self {
            message: err.to_string(),
        }
    }
}

impl From<ParseRefGeneError> for ReadWriteError {
    fn from(err: ParseRefGeneError) -> Self {
        Self {
            message: err.to_string(),
        }
    }
}
