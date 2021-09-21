use std::error::Error;
use std::fmt;
use std::num::ParseIntError;

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
