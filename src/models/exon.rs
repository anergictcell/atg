use std::fmt;

use crate::models::Frame;

#[derive(Debug, PartialEq, Eq)]
pub struct Exon {
    // u32 max value is 4,294,967,295 => This is sufficient for every human chromosome.
    // If you are working with species with chromsomes with more than 4 Mb per chromosome
    // this library will not work
    pub start: u32,
    pub end: u32,
    pub cds_start: Option<u32>,
    pub cds_end: Option<u32>,
    pub frame_offset: Frame,
}

impl Exon {
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

    pub fn end(&self) -> u32 {
        self.end
    }

    pub fn cds_start(&self) -> &Option<u32> {
        &self.cds_start
    }

    pub fn cds_end(&self) -> &Option<u32> {
        &self.cds_end
    }

    /// Returns true if the exon contains a coding sequence (CDS)
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg::models::{Exon, Frame};
    /// let mut a = Exon {start: 1, end: 2, cds_start: None, cds_end: None, frame_offset: Frame::None};
    /// assert_eq!(a.is_coding(), false);
    /// a.cds_start = Some(1);
    /// a.cds_end = Some(2);
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
