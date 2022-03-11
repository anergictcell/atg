use std::convert::TryInto;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::str::FromStr;

use crate::fasta::FastaReader;
use crate::models::{Sequence, Transcript, TranscriptWrite};
use crate::utils::errors::ReadWriteError;

/// Writes [`Transcript`]s into a `BufWriter`
///
/// The `Writer` can write the nucleotide sequences of transcripts in 3 different ways:
/// * `cds` Write the actual coding sequence of each transcript
/// * `exons` Write the processed cDNA sequence of each transcript (including 5' and 3' UTR)
/// * `transcript` Write the full genomic DNA sequence of each transcripts (including UTR and introns)
///
/// The nucleotide sequence is stranded, i.e. for `-`-strand transcripts, the sequence will always
/// be the reverse complementary of the reference DNA sequence.
///
/// All sequences are based on DNA nucleotides.
///
/// Sequences will be split into lines with 50 bp each.
///
/// # Examples
///
/// This test uses the [`Standard Transcript`](`crate::tests::transcripts::standard_transcript`)
/// for demonstration purposes.
///
/// ## Write the CDS
/// ```rust
/// use atg::fasta::Writer;
/// use atg::fasta::FastaReader;
/// use atg::models::TranscriptWrite;
/// use atg::tests::transcripts::standard_transcript;
///
/// let transcript = standard_transcript();
///
/// let mut writer = Writer::new(Vec::new());
/// writer.fasta_reader(FastaReader::from_file("tests/data/small.fasta").unwrap());
/// writer.fasta_format("cds");
/// writer.writeln_single_transcript(&transcript).unwrap();
///
/// # let output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
/// # assert_eq!(output.split('\n').collect::<Vec<&str>>()[1], "AGGCCCACTCA");
/// ```
/// This will write the Fasta data like this:
/// ```text
/// > Test-Transcript Test-Gene
/// AGGCCCACTCA
/// ```
///
/// ## Write the full transcript
/// ```rust
/// use atg::fasta::Writer;
/// use atg::fasta::FastaReader;
/// use atg::models::TranscriptWrite;
/// use atg::tests::transcripts::standard_transcript;
///
/// let transcript = standard_transcript();
///
/// let mut writer = Writer::new(Vec::new());
/// writer.fasta_reader(FastaReader::from_file("tests/data/small.fasta").unwrap());
/// writer.fasta_format("transcript");
/// writer.writeln_single_transcript(&transcript).unwrap();
///
/// # let output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
/// # assert_eq!(output.split('\n').collect::<Vec<&str>>()[1], "CACGGGGAAATGGAGGGACTGCCCAGTAGCCTCAGGACACAGGGG");
/// ```
/// This will write the Fasta data like this:
/// ```text
/// > Test-Transcript Test-Gene
/// CACGGGGAAATGGAGGGACTGCCCAGTAGCCTCAGGACACAGGGG
/// ```
pub struct Writer<W: std::io::Write> {
    inner: BufWriter<W>,
    seq_builder: SequenceBuilder,
    fasta_reader: Option<FastaReader<File>>,
    line_length: usize,
    header_template: fn(&Transcript) -> String,
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
        Writer::from_buf_writer(BufWriter::new(writer))
    }

    /// Constructs a new, empty Writer with the specified capacity.
    pub fn with_capacity(capacity: usize, writer: W) -> Self {
        Writer::from_buf_writer(BufWriter::with_capacity(capacity, writer))
    }

    /// Private constructer method to set default values
    fn from_buf_writer(writer: BufWriter<W>) -> Self {
        Writer {
            inner: writer,
            seq_builder: SequenceBuilder::Cds,
            fasta_reader: None,
            line_length: 50,
            header_template: |tx| format!("{} {}", tx.name(), tx.gene()),
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
    ///
    /// The default value is `cds`
    pub fn fasta_format(&mut self, b: &str) {
        self.seq_builder = SequenceBuilder::from_str(b).unwrap()
    }

    /// Specify how many nucleotides can be written per line
    ///
    /// By default, `ATG` will line-break nucleotide sequences
    /// after 50 nucleotides.
    ///
    /// The line length does not apply to the header line.
    pub fn line_length(&mut self, l: usize) {
        self.line_length = l
    }

    /// Specify the format of the header
    ///
    /// By default, `ATG` will use the transcript ID and gene symbol:
    /// ```text
    /// >NM_12345.1 EZH2
    /// ```
    ///
    ///
    ///
    /// # Example
    ///s
    /// ```rust
    /// use atg::fasta::Writer;
    /// use atg::fasta::FastaReader;
    /// use atg::models::{Transcript, TranscriptWrite};
    /// use atg::tests::transcripts::standard_transcript;
    ///
    /// let transcript = standard_transcript();
    ///
    /// let format_func = |tx: &Transcript| format!("{}|{} {}:{}-{}", tx.gene(), tx.name(), tx.chrom(), tx.tx_start(), tx.tx_end());
    ///
    /// let mut writer = Writer::new(Vec::new());
    /// writer.fasta_reader(FastaReader::from_file("tests/data/small.fasta").unwrap());
    /// writer.header_template(format_func);
    /// writer.writeln_single_transcript(&transcript).unwrap();
    ///
    /// # let output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
    /// # assert_eq!(output.split('\n').collect::<Vec<&str>>()[0], ">Test-Gene|Test-Transcript chr1:11-55");
    /// ```
    /// This will write the Fasta data like this:
    /// ```text
    /// >Test-Gene|Test-Transcript chr1:11-55
    /// CACGGGGAAATGGAGGGACTGCCCAGTAGCCTCAGGACACAGGGG
    /// ```
    pub fn header_template(&mut self, func: fn(&Transcript) -> String) {
        self.header_template = func
    }

    /// Flush this output stream, ensuring that all intermediately buffered contents reach their destination.
    pub fn flush(&mut self) -> Result<(), ReadWriteError> {
        match self.inner.flush() {
            Ok(res) => Ok(res),
            Err(err) => Err(ReadWriteError::new(err)),
        }
    }

    /// Unwraps this Writer<W>, returning the underlying writer.
    pub fn into_inner(self) -> Result<W, ReadWriteError> {
        match self.inner.into_inner() {
            Ok(res) => Ok(res),
            Err(err) => Err(ReadWriteError::new(err)),
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
                .write_all(format!(">{}", (self.header_template)(transcript)).as_bytes())?;

            let sequence = self.seq_builder.build(transcript, fasta_reader).to_bytes();
            // ensure line breaks after x nucleotides, as per FASTA specs
            // the last line will _not_ end in a line-break
            for line in sequence.chunks(self.line_length) {
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

/// Creates a `Sequence` from a `Transcript` and a Fasta file
///
/// SequenceBuilder reads the Fasta file and creates a Sequence
/// based on the genomic location of a transcript and the desired target regions.
///
/// Available variants:
/// * Cds - Build only the coding sequence
/// * Exons - Build the sequence of all exons - including 5' and 3' UTR or non-coding transcripts
/// * Transcript - Build the sequence of the full transcript, also including introns
enum SequenceBuilder {
    Cds,
    Exons,
    Transcript,
}

impl SequenceBuilder {
    /// Builds the actual Sequence
    pub fn build(&self, transcript: &Transcript, fasta_reader: &mut FastaReader<File>) -> Sequence {
        let segments = match self {
            SequenceBuilder::Cds => transcript.cds_coordinates(),
            SequenceBuilder::Exons => transcript.exon_coordinates(),
            SequenceBuilder::Transcript => vec![(
                transcript.chrom(),
                transcript.tx_start(),
                transcript.tx_end(),
            )],
        };

        let capacity: u32 = segments.iter().map(|x| x.2 - x.1 + 1).sum();
        let mut seq = Sequence::with_capacity(capacity.try_into().unwrap()); // unwrap can't fail

        for segment in segments {
            seq.append(
                fasta_reader
                    .read_sequence(segment.0, segment.1.into(), segment.2.into())
                    .unwrap(),
            )
        }
        if !transcript.forward() {
            seq.reverse_complement()
        }
        seq
    }
}

impl FromStr for SequenceBuilder {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "cds" => Ok(Self::Cds),
            "exons" => Ok(Self::Exons),
            "transcript" => Ok(Self::Transcript),
            _ => Err(format!("invalid fasta-format {}", s)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::FastaReader;
    use crate::tests::transcripts::standard_transcript;

    #[test]
    fn test_creating_writer() {
        let transcripts = vec![standard_transcript()];

        let mut writer = Writer::new(Vec::new());
        writer.fasta_reader(FastaReader::from_file("tests/data/small.fasta").unwrap());
        writer.fasta_format("exons");
        writer.write_transcript_vec(&transcripts).unwrap();
        let output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        assert_eq!(
            output.split('\n').collect::<Vec<&str>>()[0],
            ">Test-Transcript Test-Gene"
        );
        assert_eq!(
            output.split('\n').collect::<Vec<&str>>()[1],
            "CACGGTGGAGGCCCACTCAGAGGGG"
        );

        let mut writer = Writer::new(Vec::new());
        writer.fasta_reader(FastaReader::from_file("tests/data/small.fasta").unwrap());
        writer.fasta_format("cds");
        writer.write_transcript_vec(&transcripts).unwrap();
        let output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        assert_eq!(
            output.split('\n').collect::<Vec<&str>>()[0],
            ">Test-Transcript Test-Gene"
        );
        assert_eq!(output.split('\n').collect::<Vec<&str>>()[1], "AGGCCCACTCA");

        let mut writer = Writer::new(Vec::new());
        writer.fasta_reader(FastaReader::from_file("tests/data/small.fasta").unwrap());
        writer.fasta_format("transcript");
        writer.write_transcript_vec(&transcripts).unwrap();
        let output = String::from_utf8(writer.into_inner().unwrap()).unwrap();
        assert_eq!(
            output.split('\n').collect::<Vec<&str>>()[0],
            ">Test-Transcript Test-Gene"
        );
        assert_eq!(
            output.split('\n').collect::<Vec<&str>>()[1],
            "CACGGGGAAATGGAGGGACTGCCCAGTAGCCTCAGGACACAGGGG"
        );
    }

    #[test]
    fn test_sequence_building() {
        let transcript = standard_transcript();
        let mut reader = FastaReader::from_file("tests/data/small.fasta").unwrap();
        let seq = SequenceBuilder::Cds.build(&transcript, &mut reader);
        assert_eq!(seq.to_string(), "AGGCCCACTCA".to_string());
    }
}
