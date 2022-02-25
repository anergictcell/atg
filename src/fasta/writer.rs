use std::str::FromStr;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::models::{SequenceBuilder, Transcript, TranscriptWrite};
use crate::utils::errors::ReadWriteError;
use crate::utils::fastareader::FastaReader;

pub struct Writer<W: std::io::Write> {
    inner: BufWriter<W>,
    seq_builder: SequenceBuilder,
    fasta_reader: Option<FastaReader<File>>,
}

impl Writer<File> {
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

    pub fn with_capacity(capacity: usize, writer: W) -> Self {
        Writer {
            inner: BufWriter::with_capacity(capacity, writer),
            seq_builder: SequenceBuilder::Cds,
            fasta_reader: None,
        }
    }

    pub fn fasta_reader(&mut self, r: FastaReader<File>) {
        self.fasta_reader = Some(r)
    }

    pub fn fasta_format(&mut self, b: &str) {
        self.seq_builder = SequenceBuilder::from_str(b).unwrap()
    }

    pub fn flush(&mut self) -> Result<(), ReadWriteError> {
        match self.inner.flush() {
            Ok(res) => Ok(res),
            Err(err) => Err(ReadWriteError::from(err.to_string())),
        }
    }

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

    fn write_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error> {
        // TODO: Line-breaks after 80 characters
        if let Some(fasta_reader) = &mut self.fasta_reader {
            self.inner
                .write_all(format!(">{}:{}", transcript.gene(), transcript.name()).as_bytes())?;

            let sequence = self.seq_builder.build(transcript, fasta_reader).to_string();
            let b = sequence.as_bytes();
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