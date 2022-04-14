use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::models::{Transcript, TranscriptWrite};
use crate::utils::errors::ReadWriteError;

/// Writes [`Transcript`]s into a `BufWriter`
///
/// # Examples
///
/// ```rust
/// use std::io;
/// use atg::tests;;
/// use atg::genepredext::Writer;
/// use atg::models::TranscriptWrite;
///
/// let transcripts = vec![tests::transcripts::standard_transcript()];
///
/// let output = Vec::new(); // substitute this with proper IO (io::stdout())
/// let mut writer = Writer::new(output);
/// writer.write_transcript_vec(&transcripts);
///
/// let written_output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
/// assert_eq!(written_output.starts_with("Test-Transcript\tchr1\t"), true);
/// ```
pub struct Writer<W: std::io::Write> {
    inner: BufWriter<W>,
}

impl Writer<File> {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, ReadWriteError> {
        match File::create(path.as_ref()) {
            Ok(file) => Ok(Self::new(file)),
            Err(err) => Err(ReadWriteError::new(err)),
        }
    }
}

impl<W: std::io::Write> Writer<W> {
    /// Creates a new generic Writer for any `std::io::Read`` object
    ///
    /// Use this method when you want to write to stdout or
    /// a remote source, e.g. via HTTP
    pub fn new(writer: W) -> Self {
        Writer {
            inner: BufWriter::new(writer),
        }
    }

    pub fn with_capacity(capacity: usize, writer: W) -> Self {
        Writer {
            inner: BufWriter::with_capacity(capacity, writer),
        }
    }

    pub fn flush(&mut self) -> Result<(), std::io::Error> {
        self.inner.flush()
    }

    pub fn into_inner(self) -> Result<W, std::io::IntoInnerError<BufWriter<W>>> {
        self.inner.into_inner()
    }
}

impl<W: std::io::Write> TranscriptWrite for Writer<W> {
    /// Writes a single transcript formatted as GenePred with an extra newline
    ///
    /// This method adds an extra newline at the end of the row
    /// to allow writing multiple transcripts continuosly
    fn writeln_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error> {
        self.write_single_transcript(transcript)?;
        self.inner.write_all("\n".as_bytes())
    }

    /// Writes a single transcript formatted as GenePred
    ///
    /// Consider [`writeln_single_transcript`](Writer::writeln_single_transcript)
    /// to ensure that an extra newline is added to the output
    fn write_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error> {
        let columns: Vec<String> = Vec::from(transcript);
        // GenePredExt is similar to RefGene, just missing the first column `bin`.
        self.inner.write_all((columns[1..].join("\t")).as_bytes())
    }
}

#[cfg(test)]
mod tests {
    use super::Writer;
    use crate::genepredext::Reader;
    use crate::models::{TranscriptRead, TranscriptWrite};
    use crate::tests::transcripts;

    #[test]
    fn test_nm_001365057() {
        let transcripts = vec![transcripts::nm_001365057()];
        let mut writer = Writer::new(Vec::new());
        let _ = writer.write_transcript_vec(&transcripts);

        let output = writer.into_inner().unwrap();

        assert!(output.len() > 10);

        // Since it's a bit too complicated to compare GenePred file directly
        // this tests re-parses the GenePred data back into a Transcript
        // and compares this
        let mut reader = Reader::new(&*output);
        let read_transcripts = reader.transcripts().unwrap();

        assert_eq!(
            read_transcripts.by_name("NM_001365057.2")[0],
            &transcripts[0]
        );
    }

    #[test]
    fn test_nm_001365408() {
        let transcripts = vec![transcripts::nm_001365408()];
        let mut writer = Writer::new(Vec::new());
        let _ = writer.write_transcript_vec(&transcripts);

        let output = writer.into_inner().unwrap();

        assert!(output.len() > 10);

        // Since it's a bit too complicated to compare GenePred file directly
        // this tests re-parses the GenePred data back into a Transcript
        // and compares this
        let mut reader = Reader::new(&*output);
        let read_transcripts = reader.transcripts().unwrap();

        assert_eq!(
            read_transcripts.by_name("NM_001365408.1")[0],
            &transcripts[0]
        );
    }

    #[test]
    fn test_nm_001371720() {
        let transcripts = vec![transcripts::nm_001371720(false)];
        let mut writer = Writer::new(Vec::new());
        let _ = writer.write_transcript_vec(&transcripts);

        let output = writer.into_inner().unwrap();

        assert!(output.len() > 10);

        // Since it's a bit too complicated to compare GenePred file directly
        // this tests re-parses the GenePred data back into a Transcript
        // and compares this
        let mut reader = Reader::new(&*output);
        let read_transcripts = reader.transcripts().unwrap();

        assert_eq!(
            read_transcripts.by_name("NM_001371720.1")[0],
            &transcripts[0]
        );
    }

    #[test]
    fn test_nm_201550() {
        let transcripts = vec![transcripts::nm_201550()];
        let mut writer = Writer::new(Vec::new());
        let _ = writer.write_transcript_vec(&transcripts);

        let output = writer.into_inner().unwrap();

        assert!(output.len() > 10);

        // Since it's a bit too complicated to compare GenePred file directly
        // this tests re-parses the GenePred data back into a Transcript
        // and compares this
        let mut reader = Reader::new(&*output);
        let read_transcripts = reader.transcripts().unwrap();

        assert_eq!(read_transcripts.by_name("NM_201550.4")[0], &transcripts[0]);
    }
}
