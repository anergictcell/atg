use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::models::{Transcript, TranscriptWrite};
use crate::refgene::constants::*;
use crate::utils::errors::ParseRefGeneError;
use crate::utils::errors::ReadWriteError;

/// Writes [`Transcript`]s into a `BufWriter`
///
/// # Examples
///
/// ```rust
/// use std::io;
/// use atg::models::helper_functions;
/// use atg::refgene::Writer;
/// use atg::models::TranscriptWrite;
///
/// let transcripts = vec![helper_functions::standard_transcript()];
///
/// let output = Vec::new(); // substitute this with proper IO (io::stdout())
/// let mut writer = Writer::new(output);
/// writer.write_transcript_vec(&transcripts);
///
/// let written_output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
/// assert_eq!(written_output.starts_with("0\tTest-Transcript\tchr1\t"), true);
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

    pub fn flush(&mut self) -> Result<(), ParseRefGeneError> {
        match self.inner.flush() {
            Ok(res) => Ok(res),
            Err(err) => Err(ParseRefGeneError::from(err.to_string())),
        }
    }

    pub fn into_inner(self) -> Result<W, ParseRefGeneError> {
        match self.inner.into_inner() {
            Ok(res) => Ok(res),
            Err(err) => Err(ParseRefGeneError::from(err.to_string())),
        }
    }
}

impl<W: std::io::Write> TranscriptWrite for Writer<W> {
    /// Writes a single transcript formatted as RefGene with an extra newline
    ///
    /// This method adds an extra newline at the end of the row
    /// to allow writing multiple transcripts continuosly
    fn writeln_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error> {
        self.write_single_transcript(transcript)?;
        self.inner.write_all("\n".as_bytes())
    }

    /// Writes a single transcript formatted as RefGene
    ///
    /// Consider [`writeln_single_transcript`](Writer::writeln_single_transcript)
    /// to ensure that an extra newline is added to the output
    fn write_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error> {
        // RefGene is a tab-delimited format, each line has 16 columns
        // Defines a Vector that contains the strings for each column
        let mut columns: Vec<String> = vec!["".to_string(); N_REFGENE_COLUMNS];

        columns[TRANSCRIPT_COL] = transcript.name().to_string();
        columns[CHROMOSOME_COL] = transcript.chrom().to_string();
        columns[STRAND_COL] = transcript.strand().to_string();
        // RefGene handled start coordinates differently, so we substract 1.
        // See the comments in `instantiate_exons`.
        columns[TX_START_COL] = (transcript.tx_start() - 1).to_string();
        columns[TX_END_COL] = transcript.tx_end().to_string();
        columns[GENE_SYMBOL_COL] = transcript.gene().to_string();
        columns[CDS_START_STAT_COL] = transcript.cds_start_stat().to_string();
        columns[CDS_END_STAT_COL] = transcript.cds_end_stat().to_string();
        columns[EXON_COUNT_COL] = transcript.exon_count().to_string();

        // The bin value is not always present, defaulting to 0
        columns[BIN_COL] = match transcript.bin() {
            Some(x) => x.to_string(),
            _ => "0".to_string(),
        };

        // The score value is not always present, default to 0
        columns[SCORE_COL] = match transcript.score() {
            Some(x) => x.to_string(),
            _ => "0".to_string(),
        };

        // If the transcript does not have a CDS
        // the CDS-start and CDS-end values are set to txEnd
        columns[CDS_START_COL] = match transcript.cds_start() {
            // RefGene handled start coordinates differently, so we substract 1.
            // See the comments in `instantiate_exons`.
            Some(x) => (x - 1).to_string(),
            _ => transcript.tx_end().to_string(),
        };
        columns[CDS_END_COL] = match transcript.cds_end() {
            Some(x) => x.to_string(),
            _ => transcript.tx_end().to_string(),
        };

        // Don't ask, but some reason the refGene specs indicate that
        // there is also a trailing comma after the list of exons
        // e.g.: `12227,12721,14409,`
        columns[EXON_STARTS_COL] = transcript
            .exons()
            .iter()
            // RefGene handled start coordinates differently, so we substract 1.
            // See the comments in `instantiate_exons`.
            .map(|exon| format!("{},", (exon.start() - 1)))
            .collect();
        columns[EXON_ENDS_COL] = transcript
            .exons()
            .iter()
            .map(|exon| format!("{},", exon.end()))
            .collect();
        columns[EXON_FRAMES_COL] = transcript
            .exons()
            .iter()
            .map(|exon| format!("{},", exon.frame_offset().to_string()))
            .collect();

        self.inner.write_all((columns.join("\t")).as_bytes())
    }
}
