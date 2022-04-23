use std::convert::TryFrom;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::models::{Transcript, TranscriptRead, Transcripts};
use crate::utils::errors::{ParseRefGeneError, ReadWriteError};

/// Parses GenePredExt data and creates [`Transcript`]s.
///
/// GenePred data can be read from a file, stdin or remote sources
/// All sources are supported that provide a `std::io::Read` implementation.
///
/// # Examples
///
/// ```rust
/// use atg::genepredext::Reader;
/// use atg::models::TranscriptRead;
///
/// // create a reader from the tests GenePred file
/// let reader = Reader::from_file("tests/data/example.genepredext");
/// assert_eq!(reader.is_ok(), true);
///
/// // parse the GenePred file
/// let transcripts = reader
///     .unwrap()
///     .transcripts()
///     .unwrap();
///
/// assert_eq!(transcripts.len(), 27);
/// ```
pub struct Reader<R> {
    inner: std::io::BufReader<R>,
}

impl Reader<File> {
    /// Creates a Reader instance that reads from a File
    ///
    /// Use this method when you want to read from a GenePred file
    /// on your local file system
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg::genepredext::Reader;
    /// use atg::models::TranscriptRead;
    ///
    /// // create a reader from the tests GenePred file
    /// let reader = Reader::from_file("tests/data/example.genepredext");
    /// assert_eq!(reader.is_ok(), true);
    ///
    /// // parse the GenePred file
    /// let transcripts = reader
    ///     .unwrap()
    ///     .transcripts()
    ///     .unwrap();
    ///
    /// assert_eq!(transcripts.len(), 27);
    /// ```
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, ReadWriteError> {
        match File::open(path.as_ref()) {
            Ok(file) => Ok(Self::new(file)),
            Err(err) => Err(ReadWriteError::new(err)),
        }
    }
}

impl<R: std::io::Read> Reader<R> {
    /// creates a new Reader instance from any `std::io::Read` object
    ///
    /// Use this method when you want to read from stdin or from
    /// a remote source, e.g. via HTTP
    pub fn new(reader: R) -> Self {
        Reader {
            inner: BufReader::new(reader),
        }
    }

    /// Creates a new Reader instance with a known capcity
    ///
    /// Use this when you know the size of your GenePredExt source
    pub fn with_capacity(capacity: usize, reader: R) -> Self {
        Reader {
            inner: BufReader::with_capacity(capacity, reader),
        }
    }

    /// Returns one line of a GenePred file as `Transcript`
    ///
    /// Returns `None` if the end of the input is reached
    pub fn line(&mut self) -> Option<Result<Transcript, ParseRefGeneError>> {
        let mut line = String::new();
        match self.inner.read_line(&mut line) {
            Ok(_) => {}
            Err(x) => {
                return Some(Err(ParseRefGeneError {
                    message: x.to_string(),
                }))
            }
        }

        if line.starts_with('#') {
            return self.line();
        }

        if line.is_empty() {
            None
        } else {
            let mut cols: Vec<&str> = line.trim().split('\t').collect();
            // GenePredExt almost identical to RefGene, except
            // for the first column (`bin`). To re-use the parsing
            // functionality from the `refgene` module, we prefix
            // the line with the `bin` column value `0`.
            cols.insert(0, "0");
            Some(Transcript::try_from(cols))
        }
    }
}

impl<R: std::io::Read> TranscriptRead for Reader<R> {
    /// Reads in GenePred data and returns the final list of `Transcripts`
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg::genepredext::Reader;
    /// use atg::models::TranscriptRead;
    ///
    /// // create a reader from the tests GenePred file
    /// let reader = Reader::from_file("tests/data/example.genepredext");
    /// assert_eq!(reader.is_ok(), true);
    ///
    /// // parse the GenePred file
    /// let transcripts = reader
    ///     .unwrap()
    ///     .transcripts()
    ///     .unwrap();
    ///
    /// assert_eq!(transcripts.len(), 27);
    /// ```
    fn transcripts(&mut self) -> Result<Transcripts, ReadWriteError> {
        let mut res = Transcripts::new();
        while let Some(line) = self.line() {
            match line {
                Ok(t) => res.push(t),
                Err(x) => return Err(ReadWriteError::from(x)),
            }
        }
        Ok(res)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::TranscriptRead;
    use crate::tests::transcripts;

    #[test]
    fn test_nm_001365057() {
        let transcripts = Reader::from_file("tests/data/NM_001365057.2.genepredext")
            .unwrap()
            .transcripts()
            .unwrap();

        assert_eq!(
            transcripts.by_name("NM_001365057.2")[0],
            &transcripts::nm_001365057()
        )
    }

    #[test]
    fn test_nm_001365408() {
        let transcripts = Reader::from_file("tests/data/NM_001365408.1.genepredext")
            .unwrap()
            .transcripts()
            .unwrap();

        assert_eq!(
            transcripts.by_name("NM_001365408.1")[0],
            &transcripts::nm_001365408()
        )
    }

    #[test]
    fn test_nm_001371720() {
        let transcripts = Reader::from_file("tests/data/NM_001371720.1.genepredext")
            .unwrap()
            .transcripts()
            .unwrap();

        assert_eq!(
            transcripts.by_name("NM_001371720.1")[0],
            &transcripts::nm_001371720(false)
        )
    }

    #[test]
    fn test_nm_201550() {
        let transcripts = Reader::from_file("tests/data/NM_201550.4.genepredext")
            .unwrap()
            .transcripts()
            .unwrap();

        assert_eq!(
            transcripts.by_name("NM_201550.4")[0],
            &transcripts::nm_201550()
        )
    }
}
