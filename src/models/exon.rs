use std::fmt;

use crate::models::Frame;

/// Represents a genomic exon.
///
/// Exons can be coding and non-coding.
/// Coding exons have CDS start and end position and a [frame-offset](crate::models::Frame).
///
/// ```rust
/// use atg::models::{Exon, Frame};
///
/// let start = 1;
/// let end = 10;
/// let non_coding_exon = Exon::new(start, end, None, None, Frame::None);
///
/// assert_eq!(non_coding_exon.is_coding(), false);
/// 
/// let coding_exon = Exon::new(start, end, Some(start), Some(end), Frame::Zero);
///
/// assert_eq!(coding_exon.is_coding(), true);
/// ```
#[derive(Debug, PartialEq, Eq)]
pub struct Exon {
    // u32 max value is 4,294,967,295 => This is sufficient for every human chromosome.
    // If you are working with species with chromsomes with more than 4 Mb per chromosome
    // this library will not work
    start: u32,
    end: u32,
    cds_start: Option<u32>,
    cds_end: Option<u32>,
    frame_offset: Frame,
}

impl Exon {
    /// create a new Exon
    ///
    /// ```rust
    /// use atg::models::{Exon, Frame};
    ///
    /// let start = 1;
    /// let end = 10;
    /// let non_coding_exon = Exon::new(start, end, None, None, Frame::None);
    ///
    /// assert_eq!(non_coding_exon.is_coding(), false);
    /// 
    /// let coding_exon = Exon::new(start, end, Some(start), Some(end), Frame::Zero);
    ///
    /// assert_eq!(coding_exon.is_coding(), true);
    /// ```
    pub fn new(
        start: u32,
        end: u32,
        cds_start: Option<u32>,
        cds_end: Option<u32>,
        frame_offset: Frame,
    ) -> Exon {
        Exon {
            start,
            end,
            cds_start,
            cds_end,
            frame_offset,
        }
    }

    /// Genomic start (leftmost) position of the exon
    pub fn start(&self) -> u32 {
        self.start
    }

    /// modify the [`start`](Exon::start)
    pub fn start_mut(&mut self) -> &mut u32 {
        &mut self.start
    }

    /// Genomic end (rightmost) position of the exon
    pub fn end(&self) -> u32 {
        self.end
    }

    /// modify the [`end`](Exon::end)
    pub fn end_mut(&mut self) -> &mut u32 {
        &mut self.end
    }

    /// If the exon is coding, it contains the leftmost genomic
    /// coding nucleotide position
    pub fn cds_start(&self) -> &Option<u32> {
        &self.cds_start
    }

    /// modify the [`cds_start`](Exon::cds_start)
    pub fn cds_start_mut(&mut self) -> &mut Option<u32> {
        &mut self.cds_start
    }

    /// If the exon is coding, it contains the rightmost genomic
    /// coding nucleotide position
    pub fn cds_end(&self) -> &Option<u32> {
        &self.cds_end
    }

    /// modify the [`cds_end`](Exon::cds_end)
    pub fn cds_end_mut(&mut self) -> &mut Option<u32> {
        &mut self.cds_end
    }

    /// If the exon is coding, the [frame offset](crate::models::Frame)
    /// specifies the offset of the reading frame
    pub fn frame_offset(&self) -> &Frame {
        &self.frame_offset
    }

    /// modify the [frame offset](Exon::frame_offset)
    pub fn frame_offset_mut(&mut self) -> &mut Frame {
        &mut self.frame_offset
    }

    /// Returns true if the exon contains a coding sequence (CDS)
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg::models::{Exon, Frame};
    ///
    /// let start = 1;
    /// let end = 2;
    /// let mut a = Exon::new(start, end, None, None, Frame::None);
    /// assert_eq!(a.is_coding(), false);
    /// *a.cds_start_mut() = Some(1);
    /// *a.cds_end_mut() = Some(2);
    /// assert_eq!(a.is_coding(), true);
    /// ```
    pub fn is_coding(&self) -> bool {
        self.cds_start.is_some()
    }

    /// Returns the number of bp of the exon
    pub fn len(&self) -> u32 {
        // counting the first base as part of the exon
        self.end - self.start + 1
    }

    /// Only implemented to satisfy clippy... Exons cannot be empty.
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Returns the number of bp of the exon's coding sequence
    /// Non-coding exons have 0 bp coding sequence
    pub fn coding_len(&self) -> u32 {
        if !self.is_coding() {
            return 0;
        }
        // counting the first base as part of the ORF
        // using unwrap here is safe, because the exon is coding
        self.cds_end.unwrap() - self.cds_start.unwrap() + 1
    }

    /// Returns the coding frame of the next coding exon
    pub fn downstream_frame(&self) -> Option<Frame> {
        if !self.is_coding() {
            return None;
        }
        let frame = (3 - (self.coding_len() % 3)) % 3;
        // using unwrap here is safe, because we ensure to only have 0,1,2 frame
        match self.frame_offset + Frame::from_int(frame).unwrap() {
            Ok(x) => Some(x),
            Err(_) => None,
        }
    }

    pub fn set_frame(&mut self, frame: Frame) {
        self.frame_offset = frame;
    }
}

impl fmt::Display for Exon {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Exon ({}-{}) [{}-{}]",
            self.start,
            self.end,
            self.cds_start.unwrap_or(0),
            self.cds_end.unwrap_or(0)
        )
    }
}
