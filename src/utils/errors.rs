use std::error::Error;
use std::fmt;
use std::num::ParseIntError;

pub struct ParseGtfError {
    pub message: String,
}

impl ParseGtfError {
    pub fn new(m: &str) -> Self {
        Self {
            message: m.to_string(),
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
    pub fn new(message: &str) -> ParseRefGeneError {
        ParseRefGeneError {
            message: message.to_string(),
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
    pub fn new(message: &str) -> ParseBedError {
        ParseBedError {
            message: message.to_string(),
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
