use std::fmt;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::gtf::record::{GtfFeature, GtfRecord, GtfRecordBuilder};

use crate::models::TranscriptWrite;
use crate::models::{CdsStat, Exon, Frame, Strand, Transcript};
use crate::utils::errors::{ParseGtfError, ReadWriteError};
use crate::utils::subtract;

/// Writes [`Transcript`]s into a `BufWriter`
///
/// # Examples
///
/// ```rust
/// use std::io;
/// use atg::tests;
/// use atg::gtf::Writer;
/// use atg::models::TranscriptWrite;
///
/// let transcripts = vec![tests::transcripts::standard_transcript()];
///
/// let output = Vec::new(); // substitute this with proper IO (io::stdout())
/// let mut writer = Writer::new(output);
/// writer.write_transcript_vec(&transcripts);
///
/// let written_output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
/// assert_eq!(written_output.starts_with("chr1\t"), true);
/// ```
pub struct Writer<W: std::io::Write> {
    inner: BufWriter<W>,
    gtf_source: String,
}

impl Writer<File> {
    /// Creates a new Writer to write into a file
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
            gtf_source: "atg".to_string(),
        }
    }

    pub fn with_capacity(capacity: usize, writer: W) -> Self {
        Writer {
            inner: BufWriter::with_capacity(capacity, writer),
            gtf_source: "atg".to_string(),
        }
    }

    /// Changes the source column of the output GTF data
    ///
    /// GTF contains a source column. Use this method to
    /// globally set the source manually
    pub fn set_source(&mut self, source: &str) {
        self.gtf_source = source.to_string();
    }

    pub fn flush(&mut self) -> Result<(), ParseGtfError> {
        match self.inner.flush() {
            Ok(res) => Ok(res),
            Err(err) => Err(ParseGtfError {
                message: err.to_string(),
            }),
        }
    }

    pub fn into_inner(self) -> Result<W, ParseGtfError> {
        match self.inner.into_inner() {
            Ok(res) => Ok(res),
            Err(err) => Err(ParseGtfError {
                message: err.to_string(),
            }),
        }
    }
}

impl<W: std::io::Write> TranscriptWrite for Writer<W> {
    /// Writes a single transcript formatted as GTF with an extra newline
    ///
    /// One record can consist of multiple rows
    ///
    /// This method adds an extra newline at the end of the last GTF row
    /// to allow writing multiple transcripts continuosly
    fn writeln_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error> {
        self.write_single_transcript(transcript)?;
        self.inner.write_all("\n".as_bytes())
    }

    /// Writes a single transcript formatted as GTF
    ///
    /// One record can consist of multiple rows
    ///
    /// Consider [`writeln_single_transcript`](Writer::writeln_single_transcript)
    /// to ensure that an extra newline is added to the output
    fn write_single_transcript(&mut self, transcript: &Transcript) -> Result<(), std::io::Error> {
        let mut lines: Vec<String> = vec![];
        for gtf_record in compose_lines(transcript, &self.gtf_source) {
            lines.push(gtf_record.to_string());
        }
        self.inner.write_all(lines.join("\n").as_bytes())
    }
}

fn compose_lines(transcript: &Transcript, source: &str) -> Vec<GtfRecord> {
    let mut composer = Composer::new(transcript);
    composer.set_source(source);

    composer.build()
}

enum UtrLocation {
    Left,
    Right,
}

struct Composer<'a> {
    transcript: &'a Transcript,
    source: &'a str,
}

impl<'a> Composer<'a> {
    pub fn new(transcript: &'a Transcript) -> Self {
        Self {
            transcript,
            source: "atg",
        }
    }

    pub fn set_source(&mut self, source: &'a str) {
        self.source = source
    }

    pub fn build(&self) -> Vec<GtfRecord> {
        let mut lines = vec![self.transcript()];

        let start_codon = self.transcript.start_codon();
        let stop_codon = self.transcript.stop_codon();

        if !start_codon.is_empty() && self.transcript.cds_start_codon_stat() == CdsStat::Complete {
            lines.append(&mut self.start_codon(&start_codon));
        }
        if !stop_codon.is_empty() && self.transcript.cds_stop_codon_stat() == CdsStat::Complete {
            lines.append(&mut self.stop_codon(&stop_codon));
        }

        // One record for every exon
        for (idx, exon) in self.transcript.exons().iter().enumerate() {
            let exon_number = idx + 1;
            lines.push(self.exon(exon, &exon_number));

            if exon.is_coding() {
                // One record for the CDS
                if let Some(x) = self.cds(
                    exon,
                    // provide the first and last position of stop codon
                    // see `models::transcript::stop_codon`
                    (&stop_codon[0].0, &stop_codon.last().unwrap().1), // cannot fail, the exon is coding
                    &exon_number,
                ) {
                    lines.push(x);
                }
            }

            if self.transcript.is_coding() {
                // Add UTR where needed

                // upstream UTR
                if exon.start() < self.transcript.cds_start().unwrap() {
                    // cannot fail, the transcript is coding
                    lines.push(self.utr(exon, &UtrLocation::Left, &exon_number))
                }

                // downstream UTR
                if exon.end() > self.transcript.cds_end().unwrap() {
                    // cannot fail, the transcript is coding
                    lines.push(self.utr(exon, &UtrLocation::Right, &exon_number));
                }
            }
        }
        lines
    }

    fn transcript(&self) -> GtfRecord {
        GtfRecordBuilder::new()
            .chrom(self.transcript.chrom())
            .source(self.source)
            .feature(GtfFeature::Transcript)
            .start(self.transcript.tx_start())
            .end(self.transcript.tx_end())
            .score_option(self.transcript.score())
            .strand(self.transcript.strand())
            .frame_offset(Frame::None)
            .gene(self.transcript.gene())
            .transcript(self.transcript.name())
            .build()
            .unwrap() // cannot fail, since we control all fields for the build
    }

    fn exon(&self, exon: &Exon, exon_number: &usize) -> GtfRecord {
        GtfRecordBuilder::new()
            .chrom(self.transcript.chrom())
            .source(self.source)
            .feature(GtfFeature::Exon)
            .start(exon.start())
            .end(exon.end())
            .score_option(self.transcript.score())
            .strand(self.transcript.strand())
            .frame_offset(Frame::None)
            .gene(self.transcript.gene())
            .transcript(self.transcript.name())
            .exon_number(exon_number)
            .build()
            .unwrap() // cannot fail, since we control all fields for the build
    }

    fn cds(&self, exon: &Exon, stop_codon: (&u32, &u32), exon_number: &usize) -> Option<GtfRecord> {
        // GTF specs requires that the stop codon is excluded from the CDS
        // But this only makes sense if the CDS-end-stat is "Complete".
        // For all other types of CDS-end-stat, the full CDS is returned
        let (start, end) = match self.transcript.cds_stop_codon_stat() {
            CdsStat::Complete => {
                match subtract(
                    (&exon.cds_start().unwrap(), &exon.cds_end().unwrap()),
                    (stop_codon.0, stop_codon.1),
                ) {
                    // Stop codon is not part of this CDS
                    fragments if fragments.is_empty() => return None,

                    // Stop codon is part of CDS, transcript is on negative strand
                    fragments if self.transcript.strand() == Strand::Minus => {
                        (fragments.last().unwrap().0, fragments.last().unwrap().1)
                    }

                    // Stop codon is part of CDS on Plus strand
                    fragments => (fragments[0].0, fragments[0].1),
                }
            }
            _ => (exon.cds_start().unwrap(), exon.cds_end().unwrap()),
        };
        Some(
            GtfRecordBuilder::new()
                .chrom(self.transcript.chrom())
                .source(self.source)
                .feature(GtfFeature::CDS)
                .start(start)
                .end(end)
                .score_option(self.transcript.score())
                .strand(self.transcript.strand())
                .frame_offset(*exon.frame_offset())
                .gene(self.transcript.gene())
                .transcript(self.transcript.name())
                .exon_number(exon_number)
                .build()
                .unwrap(), // cannot fail, since we control all fields for the build
        )
    }

    /// Builds the Start Codon GTF line(s)
    fn start_codon(&self, start_codon: &[(u32, u32, Frame)]) -> Vec<GtfRecord> {
        start_codon
            .iter()
            .map(|exon| {
                GtfRecordBuilder::new()
                    .chrom(self.transcript.chrom())
                    .source(self.source)
                    .feature(GtfFeature::StartCodon)
                    .start(exon.0)
                    .end(exon.1)
                    .score_option(self.transcript.score())
                    .strand(self.transcript.strand())
                    .frame_offset(exon.2)
                    .gene(self.transcript.gene())
                    .transcript(self.transcript.name())
                    .build()
                    .unwrap() // cannot fail, since we control all fields for the build
            })
            .collect()
    }

    /// Builds the Stop Codon GTF line(s)
    fn stop_codon(&self, stop_codon: &[(u32, u32, Frame)]) -> Vec<GtfRecord> {
        stop_codon
            .iter()
            .map(|exon| {
                GtfRecordBuilder::new()
                    .chrom(self.transcript.chrom())
                    .source(self.source)
                    .feature(GtfFeature::StopCodon)
                    .start(exon.0)
                    .end(exon.1)
                    .score_option(self.transcript.score())
                    .strand(self.transcript.strand())
                    .frame_offset(exon.2)
                    .gene(self.transcript.gene())
                    .transcript(self.transcript.name())
                    .build()
                    .unwrap() // cannot fail, since we control all fields for the build
            })
            .collect()
    }

    /// Builds a UTR GTF line
    fn utr(&self, exon: &Exon, utr_location: &UtrLocation, exon_number: &usize) -> GtfRecord {
        GtfRecordBuilder::new()
            .chrom(self.transcript.chrom())
            .source(self.source)
            .feature(match (utr_location, self.transcript.strand()) {
                (UtrLocation::Left, Strand::Plus) => GtfFeature::UTR5,
                (UtrLocation::Right, Strand::Plus) => GtfFeature::UTR3,
                (UtrLocation::Left, Strand::Minus) => GtfFeature::UTR3,
                (UtrLocation::Right, Strand::Minus) => GtfFeature::UTR5,
                _ => GtfFeature::UTR,
            })
            .start(match (utr_location, exon.cds_end()) {
                // UTR start will only be different from exon start
                // when the UTR is on the right side and exon is partly coding
                (UtrLocation::Right, Some(x)) => x + 1,
                _ => exon.start(),
            })
            .end(match (utr_location, exon.cds_start()) {
                // UTR end will only be different from exon end
                // when the UTR is on the left side and exon is partly coding
                (UtrLocation::Left, Some(x)) => x - 1,
                _ => exon.end(),
            })
            .score_option(self.transcript.score())
            .strand(self.transcript.strand())
            .frame_offset(Frame::None)
            .gene(self.transcript.gene())
            .transcript(self.transcript.name())
            .exon_number(exon_number)
            .build()
            .unwrap() // cannot fail, since we control all fields for the build
    }
}

impl<'a> fmt::Display for Composer<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            self.build()
                .iter()
                .map(|line| line.to_string())
                .collect::<Vec<String>>()
                .join("\n")
        )
    }
}

#[cfg(test)]
mod tests {
    use super::Writer;
    use crate::gtf::Reader;
    use crate::models::{TranscriptRead, TranscriptWrite};
    use crate::tests::transcripts;

    #[test]
    fn test_nm_001365057() {
        let transcripts = vec![transcripts::nm_001365057()];
        let mut writer = Writer::new(Vec::new());
        let _ = writer.write_transcript_vec(&transcripts);

        let output = writer.into_inner().unwrap();

        assert!(output.len() > 10);

        // Since it's a bit too complicated to compare GTF file directly
        // this tests re-parses the GTF data back into a Transcript
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

        // Since it's a bit too complicated to compare GTF file directly
        // this tests re-parses the GTF data back into a Transcript
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
        let transcripts = vec![transcripts::nm_001371720(true)];
        let mut writer = Writer::new(Vec::new());
        let _ = writer.write_transcript_vec(&transcripts);

        let output = writer.into_inner().unwrap();

        assert!(output.len() > 10);

        // Since it's a bit too complicated to compare GTF file directly
        // this tests re-parses the GTF data back into a Transcript
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

        // Since it's a bit too complicated to compare GTF file directly
        // this tests re-parses the GTF data back into a Transcript
        // and compares this
        let mut reader = Reader::new(&*output);
        let read_transcripts = reader.transcripts().unwrap();

        assert_eq!(read_transcripts.by_name("NM_201550.4")[0], &transcripts[0]);
    }
}

/*

Missing test cases:

Split stop codon
--XXXU---AA--- +

--AA---UXXX--- -


*/
