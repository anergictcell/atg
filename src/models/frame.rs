use std::fmt;
use std::ops::Add;
use std::str::FromStr;

use serde::{Serialize, Deserialize};

/// Frame indicates the reading frame offset of an Exon
///
/// It is based on GTF nomenclature:
/// - 0 indicates that the feature begins with a whole codon at the 5' most base.
/// - 1 means that there is one extra base (the third base of a codon) before the first whole codon
/// - 2 means that there are two extra bases (the second and third bases of the codon) before the first codon.
///
/// *Important*: The reading-frame offset is handled differently in RefGene
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
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Frame {
    None, // used for non-coding exons. Converts to `-1` or `.`.
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
                Frame::None => ".",
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

    /// Returns `Frame` from a GTF value
    ///
    /// This method is the same as [`Frame::from_str`] but is present
    /// for consistency reasons
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg::models::Frame;
    /// let frame = Frame::from_gtf("1").unwrap();
    ///
    /// assert_eq!(frame, Frame::One);
    /// ```
    pub fn from_gtf(s: &str) -> Result<Self, String> {
        Frame::from_str(s)
    }

    /// Returns `Frame` from a RefGene value
    ///
    /// RefGene uses a different specification than GTF when it comes
    /// to specifying 1 or 2.
    /// RefGene specifies how many nucleotides of the first codon are
    /// present on the previous exon.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg::models::Frame;
    /// let frame = Frame::from_refgene("1").unwrap();
    ///
    /// assert_eq!(frame, Frame::Two);
    /// ```
    pub fn from_refgene(s: &str) -> Result<Self, String> {
        match s {
            "-1" => Ok(Frame::None),
            "." => Ok(Frame::None), // TODO: Not sure if its really needed
            "0" => Ok(Frame::Zero),
            "1" => Ok(Frame::Two), // yes, that's correct
            "2" => Ok(Frame::One), // this as well
            _ => Err(format!("invalid frame indicator {}", s)),
        }
    }

    /// Returns the GTF String for `Frame`
    ///
    /// This method is the same as the fmt::Display trait but is present
    /// for consistency reasons
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg::models::Frame;
    /// let frame = Frame::One.to_gtf();
    ///
    /// assert_eq!(frame, "1".to_string());
    /// ```
    pub fn to_gtf(&self) -> String {
        self.to_string()
    }

    /// Returns the RefGene String for `Frame`
    ///
    /// RefGene uses a different specification than GTF when it comes
    /// to specifying 1 or 2.
    /// RefGene specifies how many nucleotides of the first codon are
    /// present on the previous exon.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg::models::Frame;
    /// let frame = Frame::One.to_refgene();
    ///
    /// assert_eq!(frame, "2".to_string());
    /// ```
    pub fn to_refgene(&self) -> String {
        match self {
            Frame::Zero => "0",
            Frame::One => "2", // yes, that's correct
            Frame::Two => "1", // this as well
            _ => "-1",
        }
        .to_string()
    }

    /// Returns `true` if the Frame is 0, 1 or 2
    pub fn is_known(&self) -> bool {
        matches!(self, Frame::Zero | Frame::One | Frame::Two)
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
    /// Creates a [`Frame`] from a string (as in GTF format)
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
            "." => Ok(Frame::None),
            _ => Err(format!("invalid frame indicator {}", s)),
        }
    }
}
