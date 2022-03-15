use std::collections::HashMap;
use std::convert::TryFrom;
use std::fs::{read_to_string, File};
use std::io::{BufReader, Read, Seek, SeekFrom};

use std::path::Path;
use std::string::ToString;

use crate::models::Sequence;
use crate::utils::errors::FastaError;

type FastaResult<T> = Result<T, FastaError>;

/// ChromosomeIndex contains the index information for one Fasta record
///
/// It is not publicly exposed and is only called from `FastaIndex`
///
/// Is is only used to calculate the byte-offset from a nucleotide position
struct ChromosomeIndex {
    name: String,
    bases: u64,
    start: u64,
    line_bases: u64,
    line_bytes: u64,
}

impl ChromosomeIndex {
    /// Creates a new [`ChromosomeIndex`] that has only purpose:
    /// Calculate the [byte offset](ChromosomeIndex::offset) from a genomic location
    pub fn new(line: &str) -> FastaResult<Self> {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() != 5 {
            return Err(FastaError::new(format!(
                "expected 5 columns but received {}",
                cols.len()
            )));
        }
        Ok(Self {
            name: cols[0].to_string(),
            bases: cols[1].parse::<u64>()?,
            start: cols[2].parse::<u64>()?,
            line_bases: cols[3].parse::<u64>()?,
            line_bytes: cols[4].parse::<u64>()?,
        })
    }

    /// Returns the chromsome/contig name
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Returns the byte-offset in the fasta file
    /// to the specified genomic position.
    pub fn offset(&self, pos: u64) -> FastaResult<u64> {
        if pos > self.bases {
            return Err(FastaError::new(format!(
                "position {} is greater than chromome length {}",
                pos, self.bases
            )));
        }
        let line_offset = (pos - 1) % self.line_bases;
        let line_start = (pos - 1) / self.line_bases * self.line_bytes;
        Ok(self.start + line_start + line_offset)
    }
}

/// FastaIndex holds the index information for all chromosomes in the Fai file
///
/// It is not publicly exposed
///
/// Is is used by the FastaReader to access the [`ChromosomeIndex`] and look
/// up the byte offset in the Fasta file
struct FastaIndex {
    chromosomes: HashMap<String, ChromosomeIndex>,
}

impl FastaIndex {
    /// Crates a new [`FastaIndex`] by parsing the fai file
    pub fn new<P: AsRef<Path> + std::fmt::Display>(filename: P) -> FastaResult<Self> {
        let mut idx = Self {
            chromosomes: HashMap::new(),
        };
        let content = match read_to_string(&filename) {
            Ok(x) => x,
            Err(err) => {
                return Err(FastaError::new(format!(
                    "Unable to read from fasta index file {}: {}",
                    filename, err
                )))
            }
        };
        for line in content.lines() {
            let chrom = ChromosomeIndex::new(line)?;
            idx.chromosomes.insert(chrom.name().to_string(), chrom);
        }
        Ok(idx)
    }

    /// Returns the byte-offset of the fasta file for the
    /// given genomic location
    pub fn offset(&self, chrom: &str, pos: u64) -> FastaResult<u64> {
        match self.chromosomes.get(chrom) {
            Some(idx) => Ok(idx.offset(pos)?),
            None => Err(FastaError::new(format!(
                "index for {} does not exist",
                chrom
            ))),
        }
    }

    /// Returns the byte-offsets of the fasta file for both
    /// genomic start and end positions.
    ///
    /// This method is helpful to determine how many bytes to
    /// read from a fasta file for a given genomic range.
    pub fn offset_range(&self, chrom: &str, start: u64, end: u64) -> FastaResult<(u64, u64)> {
        let offset_start = self.offset(chrom, start)?;

        // seek to the position _before_ the last nucleotide, then add 1
        // to include the last nucleotide as well
        let offset_end = self.offset(chrom, end)? + 1;

        Ok((offset_start, offset_end))
    }
}

/// Provides random access to nucleotide sequences in Fasta files
///
/// It parses the Fasta-index (`*.fai`) to calculate byte offsets
///
/// Please note that this reader is different from [`GtfReader`](`crate::gtf::Reader`)
/// and [`RefGeneReader`](`crate::refgene::Reader`) in that it does _not_ parse transcripts
/// and cannot be used to build [`Transcripts`](`crate::models::Transcript`).
/// Instead, FastaReader can only read nucleotide sequences from Fasta files.
///
/// # Examples
///
/// ```rust
/// use atg;
/// use atg::fasta::FastaReader;
/// let mut reader = FastaReader::from_file("tests/data/small.fasta").unwrap();
/// let seq = reader.read_sequence("chr1", 1, 10).unwrap();
/// assert_eq!(&seq.to_string(), "GCCTCAGAGG");
/// ```
pub struct FastaReader<R> {
    inner: std::io::BufReader<R>,
    idx: FastaIndex,
}

impl FastaReader<File> {
    /// Creates a [`FastaReader`] from a fasta file location
    ///
    /// It assumes that the fasta index (fai) file has the same
    /// filename as the fasta file with an appended `.fai`
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg;
    /// use atg::fasta::FastaReader;
    /// let mut reader = FastaReader::from_file("tests/data/small.fasta").unwrap();
    /// let seq = reader.read_sequence("chr1", 1, 10).unwrap();
    /// assert_eq!(&seq.to_string(), "GCCTCAGAGG");
    /// ```
    pub fn from_file<P: AsRef<Path> + std::fmt::Display>(path: P) -> FastaResult<Self> {
        let fai_path = format!("{}.fai", path);
        Self::new(path, fai_path)
    }

    /// Returns the raw-bytes of the Fasta file for the genomic range
    ///
    /// Reads from the FastaReader and returns the raw bytes
    /// from the fasta file. The raw bytes include newline and
    /// return carriage characters from the Fasta file.
    ///
    /// # Info
    /// There are almost no use-cases to use this method. In most cases
    /// you want to use [`read_sequence`](`FastaReader::read_sequence`) instead.
    ///
    /// This method is using [`BufReader::read_exact`](`std::io::BufReader::read_exact`) internally
    /// to read the exact required amount of bytes from the Fasta file.
    pub fn read_range(&mut self, chrom: &str, start: u64, end: u64) -> FastaResult<Vec<u8>> {
        let (byte_start, byte_end) = self.idx.offset_range(chrom, start, end)?;
        self.inner.seek(SeekFrom::Start(byte_start))?;
        let capacity = usize::try_from(byte_end - byte_start)?;
        let mut buffer: Vec<u8> = vec![0; capacity];
        self.inner.read_exact(&mut buffer)?;
        Ok(buffer)
    }

    /// Returns the Nucleotide [`Sequence`] of the specified region
    ///
    /// The `Sequence` includes both `start` and `end` positions and is 1-based.
    ///
    /// This method is using [`BufReader::read_exact`](`std::io::BufReader::read_exact`) internally
    /// to read the exact required amount of bytes from the Fasta file.
    ///
    /// ```rust
    /// use atg;
    /// use atg::fasta::FastaReader;
    /// let mut reader = FastaReader::from_file("tests/data/small.fasta").unwrap();
    ///
    /// // read the nucleotide at position 150 of chromosome 5
    /// let seq = reader.read_sequence("chr5", 150, 150).unwrap();
    /// assert_eq!(&seq.to_string(), "G");
    ///
    /// // read the first 10 nucleotides of chromosome 1
    /// let seq = reader.read_sequence("chr1", 1, 10).unwrap();
    /// assert_eq!(&seq.to_string(), "GCCTCAGAGG");
    /// ```
    pub fn read_sequence(&mut self, chrom: &str, start: u64, end: u64) -> FastaResult<Sequence> {
        let raw_bytes = self.read_range(chrom, start, end)?;
        // the length of the final sequence is different to the length
        // of the raw bytes, since the raw bytes contain LF and CR characters
        // which are removed from the `Sequence`
        let length = usize::try_from(end - start)?;
        Ok(Sequence::from_raw_bytes(&raw_bytes, length)?)
    }

    /// Creates a `FastaReader` by specifying both fasta and fai file
    ///
    /// Use this method if the fasta-index file (fai) does not follow standard
    /// naming conventions. In most cases, you want to use
    /// [`from_file`](`FastaReader::from_file`) instead.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg;
    /// use atg::fasta::FastaReader;
    /// let mut reader = FastaReader::new("tests/data/small.fasta", "tests/data/small.fasta.fai").unwrap();
    /// let seq = reader.read_sequence("chr1", 1, 10).unwrap();
    /// assert_eq!(&seq.to_string(), "GCCTCAGAGG");
    /// ```
    pub fn new<P: AsRef<Path> + std::fmt::Display, P2: AsRef<Path> + std::fmt::Display>(
        fasta_path: P,
        fai_path: P2,
    ) -> FastaResult<Self> {
        let reader = match File::open(fasta_path.as_ref()) {
            Ok(x) => x,
            Err(err) => {
                return Err(FastaError::new(format!(
                    "unable to open fasta file {}: {}",
                    fasta_path, err
                )))
            }
        };
        Ok(FastaReader {
            inner: BufReader::new(reader),
            idx: FastaIndex::new(fai_path)?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_fai_reading() {
        let fai = FastaIndex::new("tests/data/small.fasta.fai").unwrap();
        assert_eq!(fai.offset("chr1", 1).unwrap(), 6);
        assert_eq!(fai.offset("chr1", 50).unwrap(), 55);
        assert_eq!(fai.offset("chr1", 51).unwrap(), 57);
        assert_eq!(fai.offset("chr1", 100).unwrap(), 106);
        assert_eq!(fai.offset("chr1", 101).unwrap(), 108);
        assert_eq!(fai.offset("chr1", 150).unwrap(), 157);
        assert_eq!(fai.offset("chr1", 151).unwrap(), 159);

        assert_eq!(fai.offset("chr2", 1).unwrap(), 218);
    }

    #[test]
    fn test_fai_errors() {
        let fai = FastaIndex::new("tests/data/small.fasta.fai").unwrap();
        assert_eq!(
            fai.offset("chr6", 1).unwrap_err().to_string(),
            "index for chr6 does not exist".to_string()
        );

        assert_eq!(
            fai.offset("chr1", 202).unwrap_err().to_string(),
            "position 202 is greater than chromome length 201".to_string()
        );

        assert_eq!(
            fai.offset("chr2", 235).unwrap_err().to_string(),
            "position 235 is greater than chromome length 234".to_string()
        );

        assert_eq!(
            fai.offset("chr3", 193).unwrap_err().to_string(),
            "position 193 is greater than chromome length 192".to_string()
        );

        assert_eq!(
            fai.offset("chr4", 150).unwrap_err().to_string(),
            "position 150 is greater than chromome length 149".to_string()
        );

        assert_eq!(
            fai.offset("chr5", 151).unwrap_err().to_string(),
            "position 151 is greater than chromome length 150".to_string()
        );
    }

    #[test]
    fn test_fasta_reading() {
        let mut fasta = FastaReader::from_file("tests/data/small.fasta").unwrap();
        let seq = fasta.read_sequence("chr1", 1, 10).unwrap();
        assert_eq!(&seq.to_string(), "GCCTCAGAGG");

        let seq = fasta.read_sequence("chr1", 51, 53).unwrap();
        assert_eq!(&seq.to_string(), "AGG");

        let seq = fasta.read_sequence("chr2", 55, 60).unwrap();
        assert_eq!(&seq.to_string(), "TCTCAT");

        let seq = fasta.read_sequence("chr4", 75, 110).unwrap();
        assert_eq!(&seq.to_string(), "GCACACCTCCTGCTTCTAACAGCAGAGCTGCCAGGC");

        let seq = fasta.read_sequence("chr1", 201, 201).unwrap();
        assert_eq!(&seq.to_string(), "G");

        let seq = fasta.read_sequence("chr1", 198, 201).unwrap();
        assert_eq!(&seq.to_string(), "GATG");

        let seq = fasta.read_sequence("chr4", 148, 149).unwrap();
        assert_eq!(&seq.to_string(), "TA");

        let seq = fasta.read_sequence("chr5", 101, 150).unwrap();
        assert_eq!(
            &seq.to_string(),
            "TGACCTGCAGGGTCGAGGAGTTGACGGTGCTGAGTTCCCTGCACTCTCAG"
        );

        let seq = fasta.read_sequence("chr5", 150, 150).unwrap();
        assert_eq!(&seq.to_string(), "G");

        let seq = fasta.read_sequence("chr5", 1, 150).unwrap();
        assert_eq!(seq.len(), 150);
    }
}
