use std::fmt;
use std::ops::Add;
use std::str::FromStr;

/// Frame indicates the reading frame offset of an Exon
///
/// It is based on GTF nomenclature:
/// - 0 indicates that the feature begins with a whole codon at the 5' most base.
/// - 1 means that there is one extra base (the third base of a codon) before the first whole codon
/// - 2 means that there are two extra bases (the second and third bases of the codon) before the first codon.
///
/// # Examples
/// ```rust
/// use std::str::FromStr;
/// use atg::models::Frame;
///
/// let frame_str = Frame::from_str("1").unwrap();
/// let frame_int = Frame::from_int(1).unwrap();
///
/// assert_eq!(frame_str, frame_int);
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Frame {
    Dot,  // used for GTF file writing. Converts to .
    None, // used for non-coding exons. Converts to -1
    Zero, // e.g. --ATG....
    One,  // e.g. --XATG....
    Two,  // e.g. --XXATG....
}

impl fmt::Display for Frame {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Frame::Dot => ".",
                Frame::None => "-1",
                Frame::Zero => "0",
                Frame::One => "1",
                Frame::Two => "2",
            }
        )
    }
}

impl Frame {
    /// Returns a Frame based on the integer frame offset
    pub fn from_int(s: u32) -> Result<Self, String> {
        match s % 3 {
            0 => Ok(Frame::Zero),
            1 => Ok(Frame::One),
            2 => Ok(Frame::Two),
            _ => Err(format!("invalid frame indicator {}", s)),
        }
    }

    fn to_int(self) -> Result<u32, String> {
        match self {
            Frame::Zero => Ok(0),
            Frame::One => Ok(1),
            Frame::Two => Ok(2),
            _ => Err("unspecified frame cannot be converted to int".to_string()),
        }
    }
}

impl Add for Frame {
    type Output = Result<Self, String>;

    /// Calculates the next resulting Frame offset
    fn add(self, other: Self) -> Result<Self, String> {
        match (self.to_int(), other.to_int()) {
            (Ok(x), Ok(y)) => Self::from_int((x + y) % 3),
            (Ok(_), _) => Ok(self),
            (_, Ok(_)) => Ok(other),
            _ => Err("unable to add two unspecified frames".to_string()),
        }
    }
}

impl FromStr for Frame {
    type Err = String;
    /// Creates a [`Frame`] from a string
    ///
    /// Only accepts the following options:
    /// - `-1`
    /// - `0`
    /// - `1`
    /// - `2`
    /// - `.`
    fn from_str(s: &str) -> Result<Self, String> {
        match s {
            "-1" => Ok(Frame::None),
            "0" => Ok(Frame::Zero),
            "1" => Ok(Frame::One),
            "2" => Ok(Frame::Two),
            "." => Ok(Frame::Dot),
            _ => Err(format!("invalid frame indicator {}", s)),
        }
    }
}
