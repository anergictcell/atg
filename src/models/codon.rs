use std::cmp;

use crate::models::transcript::Transcript;
use crate::models::Frame;
use crate::models::Strand;
use crate::utils::errors::BuildCodonError;
use crate::utils::intersect;

/// A single codon, consisting of three nucleotides
///
/// The codon can be split across multiple exons, each section is
/// represented by a `CodonFragment`.
#[derive(Debug)]
pub struct Codon {
    fragments: Vec<CodonFragment>,
}

impl Codon {
    /// Creates and returns a new `Codon` from `CodonFragments`
    ///
    /// The `CodonFragments` must be exactly 3 bp in total length.
    pub fn new(fragments: Vec<CodonFragment>) -> Result<Self, BuildCodonError> {
        let len: u32 = fragments.iter().fold(0, |sum, fgt| sum + fgt.len());

        if len != 3 {
            return Err(BuildCodonError::new("length != 3"));
        }
        Ok(Self { fragments })
    }

    /// Returns one codon from the provided transcript, starting at the given start position
    ///
    /// If the transcript lies on the Minus strand, the codon will be left of start position
    /// ```text
    /// e.g.:   123456789
    ///         ----x---- ==> start = 5
    /// Plus:   ----XXX   ==> 5-7
    /// Minus:  --XXX---- ==> 3-5
    /// ```
    ///
    /// ```rust
    /// use atg;
    /// use atg::models::Codon;
    /// use atg::tests::transcripts;
    /// let transcript = transcripts::standard_transcript();
    /// let start = transcript.cds_start().unwrap();
    /// let codon = Codon::from_transcript(&transcript, &start).unwrap();
    /// assert_eq!(codon.fragments().len(), 2);
    /// assert_eq!(codon.start(), &24);
    /// assert_eq!(codon.end(), &31);
    /// ```
    pub fn from_transcript(transcript: &Transcript, start: &u32) -> Result<Self, BuildCodonError> {
        match transcript.strand() {
            Strand::Plus => Codon::downstream(transcript, start),
            Strand::Minus => Codon::upstream(transcript, start),
            _ => Err(BuildCodonError::new("transcript with unknown Strand")),
        }
    }

    /// Returns one codon downstream (genomic left to right)
    pub fn downstream(transcript: &Transcript, start: &u32) -> Result<Self, BuildCodonError> {
        Codon::sanity_check(transcript, start)?;

        let mut len: u32 = 0;
        let mut fragments: Vec<CodonFragment> = vec![];

        for exon in transcript.exons() {
            if !exon.is_coding() {
                continue;
            }

            let cds_start = exon.cds_start().unwrap(); // cannot fail, exon is coding
            match intersect(
                (
                    &cds_start,
                    &exon.cds_end().unwrap(), // cannot fail, exon is coding
                ),
                // len == 0 => start + 2  ==> 3 bp
                // len == 1 => start + 1  ==> 2 bp
                // len == 2 => start + 0  ==> 1 bp
                (start, &(cmp::max(start, &cds_start) + (2 - len))),
            ) {
                None => continue,
                Some((start, end)) => {
                    fragments.push(CodonFragment::new(
                        transcript.chrom(),
                        start,
                        end,
                        Frame::from_int((3 - len) % 3).unwrap(), // cannot fail
                        transcript.strand(),
                    ));
                    len += end - start + 1
                }
            }
            if len >= 3 {
                break;
            }
        }
        Codon::new(fragments)
    }

    /// Returns one codon upstream (genomic right to left)
    pub fn upstream(transcript: &Transcript, start: &u32) -> Result<Self, BuildCodonError> {
        Codon::sanity_check(transcript, start)?;

        let mut len: u32 = 0;
        let mut fragments: Vec<CodonFragment> = vec![];

        for exon in transcript.exons().iter().rev() {
            if !exon.is_coding() {
                continue;
            }

            let cds_end = exon.cds_end().unwrap();
            match intersect(
                (&exon.cds_start().unwrap(), &cds_end),
                // len == 0 => end-2 - end  ==> 3 bp
                // len == 1 => end-1 - end  ==> 2 bp
                // len == 2 => end-0 - end  ==> 1 bp
                (&(cmp::min(start, &cds_end) - (3 - len - 1)), start),
            ) {
                None => continue,
                Some((start, end)) => {
                    len += end - start + 1;
                    fragments.push(CodonFragment::new(
                        transcript.chrom(),
                        start,
                        end,
                        Frame::from_int(len % 3).unwrap(),
                        transcript.strand(),
                    ));
                }
            }
            if len >= 3 {
                break;
            }
        }
        fragments.reverse();
        Codon::new(fragments)
    }

    /// Returns the leftmost position of the codon
    pub fn start(&self) -> &u32 {
        self.fragments[0].start()
    }

    /// Returns the rightmost position of the codon
    pub fn end(&self) -> &u32 {
        self.fragments.last().unwrap().end()
    }

    /// Returns a vector of `CodonFragment`s
    pub fn fragments(&self) -> &Vec<CodonFragment> {
        &self.fragments
    }

    /// Returns a tuple of `CodonFragment`-like information with
    /// start and end position, as well as Reading `Frame` offset.
    pub fn to_tuple(self) -> Vec<(u32, u32, Frame)> {
        let mut res = vec![];
        for frag in self.fragments() {
            res.push((*frag.start(), *frag.end(), frag.frame_offset()));
        }
        res
    }

    /// Performs basic sanity and QC checks on the transcript that is used
    /// for building the `Codon`.
    fn sanity_check(transcript: &Transcript, start: &u32) -> Result<(), BuildCodonError> {
        if !transcript.is_coding() {
            return Err(BuildCodonError::new("transcript is non-coding"));
        }
        if start > &transcript.cds_end().unwrap() {
            // cannot fail, transcript is coding
            return Err(BuildCodonError::new("start is downstream of the CDS"));
        }
        if start < &transcript.cds_start().unwrap() {
            // cannot fail, transcript is coding
            return Err(BuildCodonError::new("start is upstream of the CDS"));
        }
        Ok(())
    }
}

/// A `Codon` can be split across multiple exons. A `CodonFragment` represents
/// each such exon-specific fragment.
#[derive(Debug)]
pub struct CodonFragment {
    chrom: String,
    start: u32,
    end: u32,
    frame_offset: Frame,
    strand: Strand,
}

impl CodonFragment {
    /// Returns a new `CodonFragment`.
    ///
    /// CodonFragments rarely must be initiated manually. They are created
    /// when creating a `Codon`.
    pub fn new(chrom: &str, start: u32, end: u32, frame_offset: Frame, strand: Strand) -> Self {
        CodonFragment {
            chrom: chrom.to_string(),
            start,
            end,
            frame_offset,
            strand,
        }
    }

    /// Returns the length (in bp) of the `CodonFragment`.
    ///
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> u32 {
        self.end - self.start + 1
    }

    /// Returns the chromomosome of the `CodonFragment`.
    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    /// Returns the start position of the `CodonFragment`.
    pub fn start(&self) -> &u32 {
        &self.start
    }

    /// Returns the end position of the `CodonFragment`.
    pub fn end(&self) -> &u32 {
        &self.end
    }

    /// Returns the frame-offset position of the `CodonFragment`.
    pub fn frame_offset(&self) -> Frame {
        self.frame_offset
    }

    /// Returns the strand position of the `CodonFragment`.
    pub fn strand(&self) -> Strand {
        self.strand
    }
}
