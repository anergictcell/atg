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
/// use atg::genepred::Writer;
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
        // GenePred is similar to RefGene, but without the first `bin` column
        // and some other columns at the end
        self.inner.write_all((columns[1..11].join("\t")).as_bytes())
    }
}

#[cfg(test)]
mod tests {
    use super::Writer;
    use crate::models::TranscriptWrite;
    use crate::tests::transcripts;

    #[test]
    fn test_nm_001365057() {
        let transcripts = vec![transcripts::nm_001365057()];
        let mut writer = Writer::new(Vec::new());
        let _ = writer.write_transcript_vec(&transcripts);

        let output = writer.into_inner().unwrap();
        assert_eq!(
            std::str::from_utf8(&output).unwrap(),
            "NM_001365057.2\tchr9\t+\t74526554\t74600974\t74526650\t74597573\t3\t74526554,74561921,74597572,\t74526752,74562028,74600974,\n"
        );
    }

    #[test]
    fn test_nm_001365408() {
        let transcripts = vec![transcripts::nm_001365408()];
        let mut writer = Writer::new(Vec::new());
        let _ = writer.write_transcript_vec(&transcripts);

        let output = writer.into_inner().unwrap();

        assert_eq!(
            std::str::from_utf8(&output).unwrap(),
            "NM_001365408.1\tchr16\t+\t66969418\t66978999\t66972142\t66977928\t12\t66969418,66971939,66973119,66974124,66974339,66975026,66975408,66975670,66976007,66976550,66977201,66977741,\t66969778,66972144,66973261,66974258,66974598,66975125,66975549,66975751,66976152,66976640,66977274,66978999,\n"
        );
    }

    #[test]
    fn test_nm_001371720() {
        let transcripts = vec![transcripts::nm_001371720(false)];
        let mut writer = Writer::new(Vec::new());
        let _ = writer.write_transcript_vec(&transcripts);

        let output = writer.into_inner().unwrap();
        assert_eq!(
            std::str::from_utf8(&output).unwrap(),
            "NM_001371720.1\tchr1\t-\t155158299\t155162700\t155158610\t155162634\t8\t155158299,155159700,155159930,155160197,155160483,155160638,155161619,155162576,\t155158685,155159850,155160052,155160334,155160539,155161619,155162101,155162700,\n"
        );
    }

    #[test]
    fn test_nm_201550() {
        let transcripts = vec![transcripts::nm_201550()];
        let mut writer = Writer::new(Vec::new());
        let _ = writer.write_transcript_vec(&transcripts);

        let output = writer.into_inner().unwrap();
        assert_eq!(
            std::str::from_utf8(&output).unwrap(),
            "NM_201550.4\tchr12\t-\t70002343\t70004687\t70003784\t70004618\t1\t70002343,\t70004687,\n"
        );
    }
}
