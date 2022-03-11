use core::str::FromStr;
use std::convert::TryFrom;
use std::fmt;

use crate::utils::errors::AtgError;

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

/// Nucleotide is a single DNA nucleotide (A C G T N)
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Nucleotide {
    A,
    C,
    G,
    T,
    N,
}

impl Nucleotide {
    /// Crates a `Nucleotide` from a character
    pub fn new(c: &char) -> Result<Self, AtgError> {
        match c {
            'a' | 'A' => Ok(Self::A),
            'c' | 'C' => Ok(Self::C),
            'g' | 'G' => Ok(Self::G),
            't' | 'T' => Ok(Self::T),
            'n' | 'N' => Ok(Self::N),
            _ => Err(AtgError::new("Invalid nucleotide")),
        }
    }

    /// Returns the complementary nucleotide
    pub fn complement(&self) -> Self {
        match self {
            Self::A => Self::T,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::T => Self::A,
            Self::N => Self::N,
        }
    }

    // Returns the UTF-8 encoding of the Nucleotide string representation
    pub fn to_bytes(self) -> u8 {
        match self {
            Self::A => UPPERCASE_A,
            Self::C => UPPERCASE_C,
            Self::G => UPPERCASE_G,
            Self::T => UPPERCASE_T,
            Self::N => UPPERCASE_N,
        }
    }
}

impl FromStr for Nucleotide {
    type Err = AtgError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "a" | "A" => Ok(Self::A),
            "c" | "C" => Ok(Self::C),
            "g" | "G" => Ok(Self::G),
            "t" | "T" => Ok(Self::T),
            "n" | "N" => Ok(Self::N),
            _ => Err(AtgError::new("Invalid nucleotide")),
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
    type Error = AtgError;
    fn try_from(c: &char) -> Result<Self, Self::Error> {
        match c {
            'a' | 'A' => Ok(Self::A),
            'c' | 'C' => Ok(Self::C),
            'g' | 'G' => Ok(Self::G),
            't' | 'T' => Ok(Self::T),
            'n' | 'N' => Ok(Self::N),
            '\n' | '\r' => Err(AtgError::new("newline")),
            _ => panic!("invalid nucleotide {}", c),
        }
    }
}

impl TryFrom<&u8> for Nucleotide {
    type Error = AtgError;
    fn try_from(b: &u8) -> Result<Nucleotide, AtgError> {
        match b {
            &LOWERCASE_A | &UPPERCASE_A => Ok(Self::A),
            &LOWERCASE_C | &UPPERCASE_C => Ok(Self::C),
            &LOWERCASE_G | &UPPERCASE_G => Ok(Self::G),
            &LOWERCASE_T | &UPPERCASE_T => Ok(Self::T),
            &LOWERCASE_N | &UPPERCASE_N => Ok(Self::N),
            &LF | &CR => Err(AtgError::new("newline")),
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

/// A DNA sequence consisting of Nucleotides.
///
/// It provides some utility methods, like
///[`reverse_complement`](`Sequence::reverse_complement`)
pub struct Sequence {
    sequence: Vec<Nucleotide>,
}

impl FromStr for Sequence {
    type Err = AtgError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut sequence: Vec<Nucleotide> = vec![];
        for c in s.chars(){
            sequence.push(Nucleotide::new(&c)?)
        }
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

impl Default for Sequence {
    fn default() -> Self {
        Self::new()
    }
}

impl Sequence {
    /// Creates a new sequence
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg::models::Sequence;
    ///
    /// let seq = Sequence::new();
    /// assert_eq!(seq.len(), 0)
    /// ```
    pub fn new() -> Self {
        Sequence {
            sequence: Vec::new(),
        }
    }

    /// Creates a new sequence with the specified capacity
    ///
    /// Use this method if you know in advance the final size of the Sequence.
    /// It creates an empty Sequence, but one with an initial buffer that can
    /// hold capacity Nucleotides.
    ///
    /// It is important to note that although the returned Sequence has the capacity specified,
    /// the Sequence will have a zero length
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg::models::Sequence;
    ///
    /// let mut seq = Sequence::with_capacity(5);
    /// assert_eq!(seq.len(), 0);
    ///
    /// // this will not re-allocate memory, since all nucleotides fit into capacity
    /// for c in vec!['A', 'C', 'G', 'T', 'N'] {
    ///     seq.push_char(&c).unwrap();
    /// }
    /// assert_eq!(seq.len(), 5);
    /// ```
    ///
    pub fn with_capacity(capacity: usize) -> Self {
        Sequence {
            sequence: Vec::with_capacity(capacity),
        }
    }

    /// Creates a new `Sequence` from a raw bytes nucleotide sequence, ignoring newlines
    ///
    /// The `len` value is not required to be correct, it helps with allocating the right
    /// amount of memory for the Sequence.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg::models::Sequence;
    ///
    /// let seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.len(), 2);
    /// let seq = Sequence::from_raw_bytes("A\nC\r\nGT".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.len(), 4);
    /// ```
    ///
    pub fn from_raw_bytes(bytes: &[u8], len: usize) -> Result<Self, AtgError> {
        let mut seq = Self::with_capacity(len);
        for b in bytes {
            if let Ok(n) = Nucleotide::try_from(b) {
                seq.push(n)?
            }
        }
        Ok(seq)
    }

    /// Returns the length of the Sequence
    /// # Examples
    ///
    /// ```rust
    /// use atg::models::Sequence;
    ///
    /// let seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.len(), 2);
    /// ```
    ///
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Returns true if the Sequence contains no Nucleotides.
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Appends a `char` as Nucleotide to the back of a collection.
    ///
    /// # Examples
    /// ```rust
    /// use atg::models::{Nucleotide, Sequence};
    ///
    /// let mut seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.to_string(), "AC".to_string());
    ///
    /// seq.push_char(&'T').unwrap();
    /// assert_eq!(seq.to_string(), "ACT".to_string());
    pub fn push_char(&mut self, c: &char) -> Result<(), AtgError> {
        self.sequence.push(Nucleotide::try_from(c)?);
        Ok(())
    }

    /// Appends a Nucleotide to the back of a collection.
    ///
    /// # Examples
    /// ```rust
    /// use atg::models::{Nucleotide, Sequence};
    ///
    /// let mut seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.to_string(), "AC".to_string());
    ///
    /// seq.push(Nucleotide::new(&'T').unwrap()).unwrap();
    /// assert_eq!(seq.to_string(), "ACT".to_string());
    /// ```
    pub fn push(&mut self, n: Nucleotide) -> Result<(), AtgError> {
        self.sequence.push(n);
        Ok(())
    }

    /// Moves all the elements of `other` into `Self`, leaving `other` empty.
    ///
    /// # Examples
    /// ```rust
    /// use atg::models::{Nucleotide, Sequence};
    ///
    /// let mut seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.to_string(), "AC".to_string());
    ///
    /// let seq_2 = Sequence::from_raw_bytes("GT".as_bytes(), 2).unwrap();
    /// seq.append(seq_2);
    /// assert_eq!(seq.to_string(), "ACGT".to_string());
    pub fn append(&mut self, other: Sequence) {
        self.sequence.append(&mut other.into_inner())
    }

    /// Unwraps the Sequence, returning the underlying Vector of [`Nucleotide`]s
    fn into_inner(self) -> Vec<Nucleotide> {
        self.sequence
    }

    /// Changes `Self` to the complementary sequence
    ///
    /// # Examples
    /// ```rust
    /// use atg::models::Sequence;
    ///
    /// let mut seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.to_string(), "AC".to_string());
    ///
    /// seq.complement();
    /// assert_eq!(seq.to_string(), "TG".to_string());
    /// ```
    pub fn complement(&mut self) {
        for n in &mut self.sequence {
            *n = n.complement();
        }
    }

    /// Reverses the `Sequence`, in place
    ///
    /// # Examples
    /// ```rust
    /// use atg::models::Sequence;
    ///
    /// let mut seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.to_string(), "AC".to_string());
    ///
    /// seq.reverse();
    /// assert_eq!(seq.to_string(), "CA".to_string());
    /// ```
    pub fn reverse(&mut self) {
        self.sequence.reverse()
    }

    /// Changes `Self` into the reverse complement sequence
    ///
    /// # Examples
    /// ```rust
    /// use atg::models::Sequence;
    ///
    /// let mut seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.to_string(), "AC".to_string());
    ///
    /// seq.reverse_complement();
    /// assert_eq!(seq.to_string(), "GT".to_string());
    /// ```
    pub fn reverse_complement(&mut self) {
        self.reverse();
        self.complement();
    }

    /// Returns the Sequence as a byte array of UTF-8 encoded nucleotides
    ///
    /// # Examples
    /// ```rust
    /// use atg::models::Sequence;
    ///
    /// let mut seq = Sequence::from_raw_bytes("AC".as_bytes(), 2).unwrap();
    /// assert_eq!(seq.to_bytes(), [0x41, 0x43]);
    /// ```
    pub fn to_bytes(&self) -> Vec<u8> {
        self.sequence.iter().map(|n| n.to_bytes()).collect()
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
