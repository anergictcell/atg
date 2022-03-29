use core::num::TryFromIntError;
use std::error::Error;
use std::fmt;
use std::num::ParseIntError;

pub struct AtgError {
    message: String,
}

impl Error for AtgError {}

impl AtgError {
    pub fn new<S: fmt::Display>(s: S) -> Self {
        Self {
            message: s.to_string(),
        }
    }
}

impl fmt::Display for AtgError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // user-facing error
        write!(f, "{}", self.message)
    }
}

impl fmt::Debug for AtgError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // developer-facing error
        write!(f, "{{ file: {}, line: {} }}", file!(), line!())
    }
}

impl From<std::io::Error> for AtgError {
    fn from(e: std::io::Error) -> AtgError {
        AtgError {
            message: format!("IO error: {}", e),
        }
    }
}

impl From<ReadWriteError> for AtgError {
    fn from(e: ReadWriteError) -> AtgError {
        AtgError {
            message: format!("ReadWrite error: {}", e),
        }
    }
}

impl From<FastaError> for AtgError {
    fn from(e: FastaError) -> AtgError {
        AtgError {
            message: format!("Fasta error: {}", e),
        }
    }
}

impl From<String> for AtgError {
    fn from(e: String) -> AtgError {
        AtgError { message: e }
    }
}

impl From<Box<bincode::ErrorKind>> for AtgError {
    fn from(e: Box<bincode::ErrorKind>) -> AtgError {
        AtgError {
            message: e.to_string(),
        }
    }
}

pub struct ParseGtfError {
    pub message: String,
}

impl ParseGtfError {
    pub fn new<S: fmt::Display>(s: S) -> Self {
        Self {
            message: s.to_string(),
        }
    }

    pub fn from_chain(err: ParseGtfError, msg: &str) -> Self {
        Self {
            message: format!("{}\nPrevious error: {}", msg, err),
        }
    }
}

impl Error for ParseGtfError {}

impl fmt::Display for ParseGtfError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // user-facing error
        write!(
            f,
            "An error occurred while parsing the GTF input. Please check your input data.\n{}",
            self.message
        )
    }
}

impl fmt::Debug for ParseGtfError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // developer-facing error
        write!(f, "{{ file: {}, line: {} }}", file!(), line!())
    }
}

impl From<BuildTranscriptError> for ParseGtfError {
    fn from(f: BuildTranscriptError) -> ParseGtfError {
        ParseGtfError::new(&f.to_string())
    }
}

pub struct ParseRefGeneError {
    pub message: String,
}

impl ParseRefGeneError {
    pub fn new<S: fmt::Display>(s: S) -> ParseRefGeneError {
        ParseRefGeneError {
            message: s.to_string(),
        }
    }
}

impl Error for ParseRefGeneError {}

impl fmt::Display for ParseRefGeneError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // user-facing error
        write!(
            f,
            "An error occurred while parsing the RefGene input. Please check your input data\n{}",
            self.message
        )
    }
}

impl fmt::Debug for ParseRefGeneError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // developer-facing error
        write!(f, "{{ file: {}, line: {} }}", file!(), line!())
    }
}

impl From<ParseIntError> for ParseRefGeneError {
    fn from(e: ParseIntError) -> ParseRefGeneError {
        ParseRefGeneError {
            message: format!("Unable to parse an integer {}", e),
        }
    }
}

impl From<String> for ParseRefGeneError {
    fn from(e: String) -> ParseRefGeneError {
        ParseRefGeneError { message: e }
    }
}

impl From<&str> for ParseRefGeneError {
    fn from(e: &str) -> ParseRefGeneError {
        ParseRefGeneError {
            message: e.to_string(),
        }
    }
}

pub struct ParseBedError {
    pub message: String,
}

impl ParseBedError {
    pub fn new<S: fmt::Display>(s: S) -> ParseBedError {
        ParseBedError {
            message: s.to_string(),
        }
    }
}

impl Error for ParseBedError {}

impl fmt::Display for ParseBedError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // user-facing error
        write!(
            f,
            "An error occurred while parsing the RefGene input. Please check your input data\n{}",
            self.message
        )
    }
}

impl fmt::Debug for ParseBedError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // developer-facing error
        write!(f, "{{ file: {}, line: {} }}", file!(), line!())
    }
}

impl From<String> for ParseBedError {
    fn from(e: String) -> ParseBedError {
        ParseBedError { message: e }
    }
}

impl From<&str> for ParseBedError {
    fn from(e: &str) -> ParseBedError {
        ParseBedError {
            message: e.to_string(),
        }
    }
}

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
        write!(f, "Unable to build the transcript: {}", self.message)
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
        write!(f, "{}", self.message)
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

impl From<String> for ReadWriteError {
    fn from(e: String) -> ReadWriteError {
        ReadWriteError { message: e }
    }
}

#[derive(Debug)]
pub struct FastaError {
    message: String,
}

impl FastaError {
    pub fn new<S: fmt::Display>(s: S) -> Self {
        FastaError {
            message: s.to_string(),
        }
    }
}

impl Error for FastaError {}

impl fmt::Display for FastaError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl From<String> for FastaError {
    fn from(s: String) -> Self {
        FastaError { message: s }
    }
}

impl From<AtgError> for FastaError {
    fn from(err: AtgError) -> Self {
        FastaError {
            message: err.to_string(),
        }
    }
}

impl From<ParseIntError> for FastaError {
    fn from(err: ParseIntError) -> Self {
        FastaError {
            message: err.to_string(),
        }
    }
}

impl From<std::io::Error> for FastaError {
    fn from(err: std::io::Error) -> Self {
        FastaError {
            message: err.to_string(),
        }
    }
}

impl From<TryFromIntError> for FastaError {
    fn from(err: TryFromIntError) -> Self {
        FastaError {
            message: err.to_string(),
        }
    }
}

impl From<FastaError> for std::io::Error {
    fn from(err: FastaError) -> Self {
        Self::new(std::io::ErrorKind::Other, err)
    }
}
