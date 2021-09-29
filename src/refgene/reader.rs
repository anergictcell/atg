use std::convert::TryFrom;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::str::FromStr;

use crate::models::Transcripts;
use crate::refgene::constants::*;
use crate::utils::errors::ParseRefGeneError;

use crate::models;
use crate::models::{CdsStat, Exon, Frame, Strand, TranscriptBuilder};
use crate::models::{Transcript, TranscriptRead};
use crate::utils::errors::ReadWriteError;
use crate::utils::exon_cds_overlap;

/// Parses RefGene data and creates [`Transcript`]s.
///
/// RefGene data can be read from a file, stdin or remote sources
/// All sources are supported that provide a `std::io::Read` implementation.
///
/// # Examples
///
/// ```rust
/// use atg::refgene::Reader;
/// use atg::models::TranscriptRead;
///
/// // create a reader from the tests RefGene file
/// let reader = Reader::from_file("tests/data/test.refgene");
/// assert_eq!(reader.is_ok(), true);
///
/// // parse the RefGene file
/// let transcripts = reader
///     .unwrap()
///     .transcripts()
///     .unwrap();
///
/// assert_eq!(transcripts.len(), 10);
/// ```
pub struct Reader<R> {
    inner: std::io::BufReader<R>,
}

impl Reader<File> {
    /// Creates a Reader instance that reads from a File
    ///
    /// Use this method when you want to read from a RefGene file
    /// on your local file system
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg::refgene::Reader;
    /// use atg::models::TranscriptRead;
    ///
    /// // create a reader from the tests RefGene file
    /// let reader = Reader::from_file("tests/data/test.refgene");
    /// assert_eq!(reader.is_ok(), true);
    ///
    /// // parse the RefGene file
    /// let transcripts = reader
    ///     .unwrap()
    ///     .transcripts()
    ///     .unwrap();
    ///
    /// assert_eq!(transcripts.len(), 10);
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
    /// Use this when you know the size of your RefGene source
    pub fn with_capacity(capacity: usize, reader: R) -> Self {
        Reader {
            inner: BufReader::with_capacity(capacity, reader),
        }
    }

    /// Returns one line of a RefGene file as `Transcript`
    ///
    /// This method should rarely be used. GTF files can contain unordered
    /// records and handling lines individually is rarely desired.
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
            let cols: Vec<&str> = line.trim().split('\t').collect();
            Some(Transcript::try_from(cols))
        }
    }
}

impl<R: std::io::Read> TranscriptRead for Reader<R> {
    /// Reads in RefGene data and returns the final list of `Transcripts`
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg::refgene::Reader;
    /// use atg::models::TranscriptRead;
    ///
    /// // create a reader from the tests RefGene file
    /// let reader = Reader::from_file("tests/data/test.refgene");
    /// assert_eq!(reader.is_ok(), true);
    ///
    /// // parse the RefGene file
    /// let transcripts = reader
    ///     .unwrap()
    ///     .transcripts()
    ///     .unwrap();
    ///
    /// assert_eq!(transcripts.len(), 10);
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

impl<R: std::io::Read> IntoIterator for Reader<R> {
    type Item = Result<Transcript, ParseRefGeneError>;
    type IntoIter = TranscriptIterator<R>;

    fn into_iter(self) -> TranscriptIterator<R> {
        TranscriptIterator::new(self)
    }
}

pub struct TranscriptIterator<R> {
    inner: Reader<R>,
}

impl<R> TranscriptIterator<R> {
    fn new(reader: Reader<R>) -> TranscriptIterator<R> {
        TranscriptIterator { inner: reader }
    }
}

impl<R: std::io::Read> Iterator for TranscriptIterator<R> {
    type Item = Result<Transcript, ParseRefGeneError>;

    fn next(&mut self) -> Option<Result<Transcript, ParseRefGeneError>> {
        self.inner.line()
    }
}

impl TryFrom<Vec<&str>> for Transcript {
    type Error = ParseRefGeneError;

    /// Returns a ```Transcript``` based on features of the GTF file,
    /// belonging to one transcript
    fn try_from(cols: Vec<&str>) -> Result<Self, ParseRefGeneError> {
        if cols.len() != N_REFGENE_COLUMNS {
            return Err(ParseRefGeneError {
                message: format!(
                    "Invalid number of columns in line\nvv\n{}\n^^",
                    cols.join("\t")
                ),
            });
        }

        let bin = match cols[BIN_COL].parse::<u16>() {
            Ok(x) => Some(x),
            _ => None,
        };
        let strand = match Strand::from_str(cols[STRAND_COL]) {
            Ok(x) => x,
            Err(message) => return Err(ParseRefGeneError { message }),
        };
        let mut exons = instantiate_exons(&cols)?;
        let cds_start_stat = match CdsStat::from_str(cols[CDS_START_STAT_COL]) {
            Ok(x) => x,
            Err(message) => return Err(ParseRefGeneError { message }),
        };

        let cds_end_stat = match CdsStat::from_str(cols[CDS_END_STAT_COL]) {
            Ok(x) => x,
            Err(message) => return Err(ParseRefGeneError { message }),
        };
        let score = match cols[SCORE_COL].parse::<f32>() {
            Ok(x) => Some(x),
            _ => None,
        };

        let mut transcript = TranscriptBuilder::new()
            .bin(bin)
            .name(cols[TRANSCRIPT_COL])
            .chrom(&models::parse_chrom(cols[CHROMOSOME_COL]))
            .strand(strand)
            .gene(cols[GENE_SYMBOL_COL])
            .cds_start_stat(cds_start_stat)
            .cds_end_stat(cds_end_stat)
            .score(score)
            .build()
            .unwrap();
        transcript.append_exons(&mut exons);
        Ok(transcript)
    }
}

/// Returns a Vector of `Exon`s
///
/// It parses the columns that specify the start and end positions
/// of an exon, as well as frame etc.
/// It also check for CDS overlap.
fn instantiate_exons(cols: &[&str]) -> Result<Vec<Exon>, ParseRefGeneError> {
    // RefGene specifies the number of exons.
    // This assumes this number to be correct.
    let exon_count = cols[EXON_COUNT_COL].parse::<usize>().unwrap();

    // Create an empty vector to hold all exons
    let mut exons: Vec<Exon> = Vec::with_capacity(exon_count);

    // create temporary vectors for the start, end coordinates
    // and the frame-offsets for every exon
    let starts: Vec<&str> = cols[EXON_STARTS_COL]
        .trim_end_matches(',')
        .split(',')
        .collect();
    let ends: Vec<&str> = cols[EXON_ENDS_COL]
        .trim_end_matches(',')
        .split(',')
        .collect();
    let frame_offsets: Vec<&str> = cols[EXON_FRAMES_COL]
        .trim_end_matches(',')
        .split(',')
        .collect();

    // Parse the cds-start and end positions
    // Also check if the transcript is coding at all
    let (coding, cds_start, cds_end) = match (
        cols[CDS_START_COL].parse::<u32>(),
        cols[CDS_END_COL].parse::<u32>(),
    ) {
        (Ok(start), Ok(end)) => (true, Some(start), Some(end)),
        _ => (false, None, None),
    };

    // Iterate through all exons and create `Exon` instances
    for i in 0..exon_count {
        let start = starts
            .get(i)
            .ok_or("Too few exon starts in input")?
            // RefGene start positions are excluded
            // but atg includes them. So we change the start coordinate
            .parse::<u32>()?
            + 1;
        let end = ends
            .get(i)
            .ok_or("Too few exon ends in input")?
            .parse::<u32>()?;

        // Check of the exon is coding
        let exon_cds = if coding {
            exon_cds_overlap(
                &start,
                &end,
                // RefGene start positions are excluded
                // but atg includes them. So we change the start coordinate
                // unwrapping is safe here, this is checked previously
                &(cds_start.unwrap() + 1),
                &cds_end.unwrap(),
            )
        } else {
            (None, None)
        };

        exons.push(Exon::new(
            start,
            end,
            exon_cds.0,
            exon_cds.1,
            Frame::from_refgene(frame_offsets.get(i).ok_or("Too few exon Frame offsets")?)?,
        ));
    }
    Ok(exons)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_exons_no_cds() {
        let cols = vec![
            "585",
            "NR_046018.2",
            "chr1",
            "+",
            "11873",
            "14409",
            "14409",
            "14409",
            "3",
            "11873,12612,13220,",
            "12227,12721,14409,",
            "0",
            "DDX11L1",
            "unk",
            "unk",
            "-1,-1,-1,",
        ];
        let exons = instantiate_exons(&cols).unwrap();

        assert_eq!(exons.len(), 3);

        assert_eq!(exons[0].start(), 11874);
        assert_eq!(exons[0].end(), 12227);
        assert_eq!(*exons[0].cds_start(), None);
        assert_eq!(*exons[0].cds_end(), None);

        assert_eq!(exons[1].start(), 12613);
        assert_eq!(exons[1].end(), 12721);
        assert_eq!(*exons[1].cds_start(), None);
        assert_eq!(*exons[1].cds_end(), None);

        assert_eq!(exons[2].start(), 13221);
        assert_eq!(exons[2].end(), 14409);
        assert_eq!(*exons[2].cds_start(), None);
        assert_eq!(*exons[2].cds_end(), None);
    }

    #[test]
    fn test_missing_exon_stop() {
        let cols = vec![
            "585",
            "NR_046018.2",
            "chr1",
            "+",
            "11873",
            "14409",
            "14409",
            "14409",
            "3",
            "11873,12612,13220",
            "12227,12721,",
            "0",
            "DDX11L1",
            "unk",
            "unk",
            "-1,-1,-1,",
        ];
        let exons = instantiate_exons(&cols);

        assert_eq!(exons.is_err(), true);
        assert_eq!(
            exons.unwrap_err().message,
            "Too few exon ends in input".to_string()
        );
    }

    #[test]
    fn test_missing_exon_start() {
        let cols = vec![
            "585",
            "NR_046018.2",
            "chr1",
            "+",
            "11873",
            "14409",
            "14409",
            "14409",
            "3",
            "11873,12612,",
            "12227,12721,14409,",
            "0",
            "DDX11L1",
            "unk",
            "unk",
            "-1,-1,-1,",
        ];
        let exons = instantiate_exons(&cols);

        assert_eq!(exons.is_err(), true);
        assert_eq!(
            exons.unwrap_err().message,
            "Too few exon starts in input".to_string()
        );
    }

    #[test]
    fn test_missing_exon_frame() {
        let cols = vec![
            "585",
            "NR_046018.2",
            "chr1",
            "+",
            "11873",
            "14409",
            "14409",
            "14409",
            "3",
            "11873,12612,13220",
            "12227,12721,14409,",
            "0",
            "DDX11L1",
            "unk",
            "unk",
            "-1,-1,",
        ];
        let exons = instantiate_exons(&cols);

        assert_eq!(exons.is_err(), true);
        assert_eq!(
            exons.unwrap_err().message,
            "Too few exon Frame offsets".to_string()
        );
    }

    #[test]
    fn test_wrong_exon_frame() {
        let cols = vec![
            "585",
            "NR_046018.2",
            "chr1",
            "+",
            "11873",
            "14409",
            "14409",
            "14409",
            "3",
            "11873,12612,13220",
            "12227,12721,14409,",
            "0",
            "DDX11L1",
            "unk",
            "unk",
            "-1,/,-1",
        ];
        let exons = instantiate_exons(&cols);

        assert_eq!(exons.is_err(), true);
        assert_eq!(
            exons.unwrap_err().message,
            "invalid frame indicator /".to_string()
        );
    }
}
