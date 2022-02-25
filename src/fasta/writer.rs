use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::str::FromStr;

use crate::fasta::FastaReader;
use crate::models::{SequenceBuilder, Transcript, TranscriptWrite};
use crate::utils::errors::ReadWriteError;

/// Writes [`Transcript`]s into a `BufWriter`
///
/// # Examples
///
/// ```no_run
/// use atg::tests;
/// use atg::fasta::Writer;
/// use atg::fasta::FastaReader;
/// use atg::models::TranscriptWrite;
///
/// let transcripts = vec![tests::transcripts::standard_transcript()];
///
/// let output = Vec::new(); // substitute this with proper IO (io::stdout())
/// let mut writer = Writer::new(output);
/// // specify the reference genome fasta file
/// writer.fasta_reader(FastaReader::from_file("/path/to/hg19.fasta").unwrap());
/// writer.write_transcript_vec(&transcripts);
/// ```
pub struct Writer<W: std::io::Write> {
    inner: BufWriter<W>,
    seq_builder: SequenceBuilder,
    fasta_reader: Option<FastaReader<File>>,
}

impl Writer<File> {
    /// Creates a new Writer to write into a file
    pub fn from_file<P: AsRef<Path> + std::fmt::Display>(path: P) -> Result<Self, ReadWriteError> {
        match File::create(path.as_ref()) {
            Ok(file) => Ok(Self::new(file)),
            Err(err) => Err(ReadWriteError::new(format!(
                "unable to open file {} for writing: {}",
                path, err
            ))),
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
            seq_builder: SequenceBuilder::Cds,
            fasta_reader: None,
        }
    }

    /// Constructs a new, empty Writer with the specified capacity.
    pub fn with_capacity(capacity: usize, writer: W) -> Self {
        Writer {
            inner: BufWriter::with_capacity(capacity, writer),
            seq_builder: SequenceBuilder::Cds,
            fasta_reader: None,
        }
    }

    /// Specify a [`FastaReader'](`crate::fasta::FastaReader`) to retrieve
    /// the reference genome sequence.
    ///
    /// You must set a `fasta_reader`, since the `Writer` does not have any
    /// information about the reference genome to use.
    ///
    /// ```rust
    /// use atg::fasta::Writer;
    /// use atg::fasta::FastaReader;
    ///
    /// let output = Vec::new(); // substitute this with proper IO (io::stdout())
    /// let mut writer = Writer::new(output);
    /// // specify the reference genome fasta file
    /// writer.fasta_reader(FastaReader::from_file("tests/data/small.fasta").unwrap());
    /// ```
    pub fn fasta_reader(&mut self, r: FastaReader<File>) {
        self.fasta_reader = Some(r)
    }

    /// Specify how to write the Transcript sequence as Fasta
    ///
    /// You can write either the full length transcript's sequence (exons + introns)
    /// `'transcript'`, write the sequence of all exons
    /// (including non coding, e.g. 5'UTR) `'exons'` or only write the
    /// coding sequence (from `ATG` to the Stop codon.) `'cds'`.
    pub fn fasta_format(&mut self, b: &str) {
        self.seq_builder = SequenceBuilder::from_str(b).unwrap()
    }

    /// Flush this output stream, ensuring that all intermediately buffered contents reach their destination.
    pub fn flush(&mut self) -> Result<(), ReadWriteError> {
        match self.inner.flush() {
            Ok(res) => Ok(res),
            Err(err) => Err(ReadWriteError::from(err.to_string())),
        }
    }

    /// Unwraps this Writer<W>, returning the underlying writer.
    pub fn into_inner(self) -> Result<W, ReadWriteError> {
        match self.inner.into_inner() {
            Ok(res) => Ok(res),
            Err(err) => Err(ReadWriteError::from(err.to_string())),
        }
    }
}

impl<W: std::io::Write> TranscriptWrite for Writer<W> {
    /// Writes a single transcript FASTA sequence with an extra newline
    ///
    /// This method adds an extra newline at the end of the row
    /// to allow writing multiple transcripts continuosly
    fn writeln_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error> {
        self.write_single_transcript(transcript)?;
        self.inner.write_all("\n".as_bytes())
    }

    /// Writes a single transcript FASTA sequence _without_ an extra newline
    ///
    /// This is almost never what you want to do. You most likely want to use
    /// [`writeln_single_transcript`](`Writer::writeln_single_transcript`) instead.
    fn write_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error> {
        if let Some(fasta_reader) = &mut self.fasta_reader {
            self.inner
                .write_all(format!(">{}:{}", transcript.gene(), transcript.name()).as_bytes())?;

            let sequence = self.seq_builder.build(transcript, fasta_reader).to_string();
            let b = sequence.as_bytes();

            // ensure line breaks after 80 nucleotides, as per FASTA specs
            // the last line will _not_ end in a line-break
            for line in b.chunks(80) {
                self.inner.write_all("\n".as_bytes())?;
                self.inner.write_all(line)?;
            }
            Ok(())
        } else {
            Err(std::io::Error::new(
                std::io::ErrorKind::Other,
                "no fasta reader specified",
            ))
        }
    }
}
