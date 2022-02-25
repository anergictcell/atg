/// https://www.biostars.org/p/98885/
use core::slice::Chunks;
use core::str::FromStr;
use std::convert::TryFrom;
use std::convert::TryInto;
use std::fmt;
use std::fs::File;

use crate::fasta::FastaReader;
use crate::models::Transcript;

// UTF-8 encoding of all nucleotides
const UPPERCASE_A: u8 = 0x41;
const UPPERCASE_C: u8 = 0x43;
const UPPERCASE_G: u8 = 0x47;
const UPPERCASE_T: u8 = 0x54;
const UPPERCASE_N: u8 = 0x4e;
const LOWERCASE_A: u8 = 0x61;
const LOWERCASE_C: u8 = 0x63;
const LOWERCASE_G: u8 = 0x67;
const LOWERCASE_T: u8 = 0x74;
const LOWERCASE_N: u8 = 0x64;

const LF: u8 = 0xa;
const CR: u8 = 0xd;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Nucleotide {
    A,
    C,
    G,
    T,
    N,
}

impl Nucleotide {
    pub fn new(c: &char) -> Result<Self, String> {
        match c {
            'a' | 'A' => Ok(Self::A),
            'c' | 'C' => Ok(Self::C),
            'g' | 'G' => Ok(Self::G),
            't' | 'T' => Ok(Self::T),
            'n' | 'N' => Ok(Self::N),
            _ => Err("Invalid nucleotide".to_string()),
        }
    }

    pub fn complement(&self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
            Self::N => Self::N,
        }
    }
}

impl FromStr for Nucleotide {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "a" | "A" => Ok(Self::A),
            "c" | "C" => Ok(Self::C),
            "g" | "G" => Ok(Self::G),
            "t" | "T" => Ok(Self::T),
            "n" | "N" => Ok(Self::N),
            _ => Err("Invalid nucleotide".to_string()),
        }
    }
}

impl fmt::Display for Nucleotide {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::A => "A",
                Self::C => "C",
                Self::G => "G",
                Self::T => "T",
                Self::N => "N",
            }
        )
    }
}

impl TryFrom<&char> for Nucleotide {
    type Error = String;
    fn try_from(c: &char) -> Result<Self, Self::Error> {
        match c {
            'a' | 'A' => Ok(Self::A),
            'c' | 'C' => Ok(Self::C),
            'g' | 'G' => Ok(Self::G),
            't' | 'T' => Ok(Self::T),
            'n' | 'N' => Ok(Self::N),
            '\n' | '\r' => Err("newline".to_string()),
            _ => panic!("invalid nucleotide {}", c),
        }
    }
}

impl TryFrom<&u8> for Nucleotide {
    type Error = String;
    fn try_from(b: &u8) -> Result<Nucleotide, String> {
        match b {
            &LOWERCASE_A | &UPPERCASE_A => Ok(Self::A),
            &LOWERCASE_C | &UPPERCASE_C => Ok(Self::C),
            &LOWERCASE_G | &UPPERCASE_G => Ok(Self::G),
            &LOWERCASE_T | &UPPERCASE_T => Ok(Self::T),
            &LOWERCASE_N | &UPPERCASE_N => Ok(Self::N),
            &LF | &CR => Err("newline".to_string()),
            _ => panic!("invalid nucleotide {}", b),
        }
    }
}

impl From<&Nucleotide> for char {
    fn from(n: &Nucleotide) -> Self {
        match n {
            Nucleotide::A => 'A',
            Nucleotide::C => 'C',
            Nucleotide::G => 'G',
            Nucleotide::T => 'T',
            Nucleotide::N => 'N',
        }
    }
}

pub struct Sequence {
    sequence: Vec<Nucleotide>,
}

impl FromStr for Sequence {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let sequence: Vec<Nucleotide> = s.chars().map(|c| Nucleotide::new(&c).unwrap()).collect();
        Ok(Self { sequence })
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut s = String::with_capacity(self.len());
        for n in &self.sequence {
            s.push(n.into())
        }
        write!(f, "{}", s)
    }
}

impl Sequence {
    pub fn with_capacity(capacity: usize) -> Self {
        Sequence {
            sequence: Vec::with_capacity(capacity),
        }
    }
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    pub fn push(&mut self, c: &char) -> Result<(), String> {
        self.sequence.push(Nucleotide::try_from(c)?);
        Ok(())
    }

    pub fn push_nucleotide(&mut self, n: Nucleotide) -> Result<(), String> {
        self.sequence.push(n);
        Ok(())
    }

    pub fn append(&mut self, s: Sequence) {
        self.sequence.append(&mut s.into_inner())
    }

    pub fn into_inner(self) -> Vec<Nucleotide> {
        self.sequence
    }

    pub fn from_raw_bytes(bytes: &[u8], len: usize) -> Result<Self, String> {
        let mut seq = Self::with_capacity(len);
        for b in bytes {
            if let Ok(n) = Nucleotide::try_from(b) {
                seq.push_nucleotide(n)?
            }
        }
        Ok(seq)
    }

    pub fn complement(&mut self) {
        for n in &mut self.sequence {
            *n = n.complement();
        }
    }

    pub fn reverse(&mut self) {
        self.sequence.reverse()
    }

    pub fn reverse_complement(&mut self) {
        self.reverse();
        self.complement();
    }

    pub fn format(self, line_length: usize) -> Vec<Sequence> {
        println!("{}", line_length);
        vec![self]
    }

    pub fn chunks(&self, chunk_size: usize) -> Chunks<'_, Nucleotide> {
        self.sequence.chunks(chunk_size)
    }
}

pub enum SequenceBuilder {
    Cds,
    Exons,
    Transcript,
}

impl SequenceBuilder {
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
        let mut seq = Sequence::with_capacity(capacity.try_into().unwrap());

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

    #[test]
    fn test_create_sequence() {
        let s = "ATCGACGATCGATCGATGAGCGATCGACGATCGCGCTATCGCTA";
        let seq = Sequence::from_str(&s).unwrap();

        assert_eq!(seq.len(), 44);
        assert_eq!(seq.to_string(), s.to_string())
    }
}
