use crate::utils::errors::FastaError;
use std::collections::HashMap;
use std::convert::TryFrom;


use std::fs::{read_to_string, File};
use std::io::{BufReader, Read, Seek, SeekFrom};

use std::path::Path;
use std::string::ToString;

use crate::models::Sequence;

type FastaResult<T> = Result<T, FastaError>;

struct ChromosomeIndex {
    name: String,
    bases: u64,
    start: u64,
    line_bases: u64,
    line_bytes: u64,
}

impl ChromosomeIndex {
    /// Creates a new `ChromosomeIndex` that has only purpose:
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
            FastaError::new(format!(
                "start position {} is greater than chromome length {}",
                pos, self.bases
            ));
        }
        let line_offset = (pos - 1) % self.line_bases;
        let line_start = (pos - 1) / self.line_bases * self.line_bytes;
        Ok(self.start + line_start + line_offset)
    }
}

/// Holds the fasta index (fai) information
/// to help with direct access within the fasta file
struct FastaIndex {
    chromosomes: HashMap<String, ChromosomeIndex>,
}

impl FastaIndex {
    /// Crates a new `FastaIndex` by parsing the fai file
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
        let offset_end = self.offset(chrom, end + 1)?;

        Ok((offset_start, offset_end))
    }
}

pub struct FastaReader<R> {
    inner: std::io::BufReader<R>,
    idx: FastaIndex,
}

impl FastaReader<File> {
    /// Creates a `FastaReader` from a fasta file location
    ///
    /// It assumes that the fasta index (fai) file has the same
    /// filename as the fasta file with an appended `.fai`
    pub fn from_file<P: AsRef<Path> + std::fmt::Display>(path: P) -> FastaResult<Self> {
        match File::open(path.as_ref()) {
            Ok(file) => {
                let fai_path = format!("{}.fai", path);
                Ok(Self::new(file, fai_path)?)
            }
            Err(err) => Err(FastaError::new(format!(
                "unable to open fasta file {}: {}",
                path, err
            ))),
        }
    }

    /// Reads from the FastaReader and returns the raw bytes
    /// from the fasta file. The raw bytes include newline and
    /// return carriage characters from the Fasta file.
    pub fn read_range(&mut self, chrom: &str, start: u64, end: u64) -> FastaResult<Vec<u8>> {
        let (byte_start, byte_end) = self.idx.offset_range(chrom, start, end)?;
        self.inner.seek(SeekFrom::Start(byte_start))?;
        let capacity = usize::try_from(byte_end - byte_start)?;
        let mut buffer: Vec<u8> = vec![0; capacity];
        self.inner.read_exact(&mut buffer)?;
        Ok(buffer)
    }

    /// Reads from the FastaReader and returns a `Sequence` of the
    /// specified region. The `Sequence` includes both `start` and `end`
    /// position and is 1-based.
    pub fn read_sequence(&mut self, chrom: &str, start: u64, end: u64) -> FastaResult<Sequence> {
        let raw_bytes = self.read_range(chrom, start, end)?;
        let length = usize::try_from(end - start)?;
        Ok(Sequence::from_raw_bytes(&raw_bytes, length)?)
    }
}

impl<R: std::io::Read> FastaReader<R> {
    /// creates a new FastaReader instance from any `std::io::Read` object
    pub fn new<P: AsRef<Path> + std::fmt::Display>(reader: R, fai_path: P) -> FastaResult<Self> {
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
    }
}
