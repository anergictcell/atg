use std::fmt;

use crate::models::codon::Codon;
use crate::models::utils::{CdsStat, Strand};
use crate::models::{Exon, Frame};
use crate::utils::errors::BuildTranscriptError;

/// Transcript is the central data structure of `atg`
///
/// It holds the genomic representation of transcript. The coordinates are 1-based
/// and both start and end coordinate are included.
///
/// A transcript contains exons, some of which may be coding.
/// Transcripts are directional and the direction is encoded through the strand.

/// `Transcript`s should be created using `TranscriptBuilder`
#[derive(Debug)]
pub struct Transcript {
    bin: Option<u16>,
    name: String,
    chrom: String,
    strand: Strand,
    cds_start_stat: CdsStat,
    cds_end_stat: CdsStat,
    exons: Vec<Exon>,
    gene_symbol: String,
    score: Option<f32>,
}

impl Transcript {
    pub fn bin(&self) -> &Option<u16> {
        &self.bin
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn gene(&self) -> &str {
        &self.gene_symbol
    }

    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    pub fn strand(&self) -> Strand {
        self.strand
    }

    pub fn cds_start_stat(&self) -> CdsStat {
        self.cds_start_stat
    }

    pub fn cds_end_stat(&self) -> CdsStat {
        self.cds_end_stat
    }

    pub fn cds_start_codon_stat(&self) -> CdsStat {
        match self.strand {
            Strand::Minus => self.cds_end_stat(),
            _ => self.cds_start_stat(),
        }
    }

    pub fn cds_stop_codon_stat(&self) -> CdsStat {
        match self.strand {
            Strand::Minus => self.cds_start_stat(),
            _ => self.cds_end_stat(),
        }
    }

    pub fn exons(&self) -> &Vec<Exon> {
        &self.exons
    }

    pub fn exons_mut(&mut self) -> &mut Vec<Exon> {
        &mut self.exons
    }

    pub fn score(&self) -> Option<f32> {
        self.score
    }

    pub fn push_exon(&mut self, exon: Exon) {
        self.exons.push(exon)
    }

    pub fn append_exons(&mut self, exons: &mut Vec<Exon>) {
        self.exons.append(exons)
    }

    pub fn set_cds_start_stat(&mut self, stat: CdsStat) {
        self.cds_start_stat = stat;
    }

    pub fn set_cds_end_stat(&mut self, stat: CdsStat) {
        self.cds_end_stat = stat;
    }

    /// Returns if the Transcript is annotated on the + strand
    ///
    /// Note: If the strand is unknown, this method will also return `true`.
    pub fn forward(&self) -> bool {
        match self.strand {
            Strand::Plus => true,
            // most non-stranded features are displayed as if
            // they are on the + strand
            Strand::Unknown => true,
            Strand::Minus => false,
        }
    }

    pub fn exon_count(&self) -> usize {
        self.exons.len()
    }

    pub fn tx_start(&self) -> u32 {
        self.exons[0].start
    }

    pub fn tx_end(&self) -> u32 {
        self.exons[self.exons.len() - 1].end
    }

    pub fn cds_start(&self) -> Option<u32> {
        for exon in &self.exons {
            if let Some(x) = exon.cds_start {
                return Some(x);
            };
        }
        None
    }

    pub fn cds_end(&self) -> Option<u32> {
        for exon in self.exons.iter().rev() {
            if let Some(x) = exon.cds_end {
                return Some(x);
            };
        }
        None
    }

    pub fn is_coding(&self) -> bool {
        for exon in &self.exons {
            if exon.is_coding() {
                return true;
            }
        }
        false
    }

    /// Get a vector of all exons that span the start exon
    ///
    /// The start codon can be split across multiple exons
    /// and thus the coordinates of it can't be easily calculated
    /// Taking into account the following options, by splitting the
    /// start codon into multiple exons
    /// ```text
    ///       1....   2....   3....   4....   5....
    ///       12345   12345   12345   12345   12345
    ///    ---=====---===XX---XXXXX---XXXX=---=====---
    /// 1. ---=====---ATGXX---XXXXX---XXXX=---=====--- >> all in one exon
    /// 2. ---=====---=ATGX---XXXXX---XXXX=---=====--- >> all in one exon
    /// 3. ---=====---==ATG---XXXXX---XXXX=---=====--- >> all in one exon
    /// 4. ---=====---===AT---GXXXX---XXXX=---=====--- >> split
    /// 5. ---=====---====A---TGXXX---XXXX=---=====--- >> split
    /// 6. ---=====---====A-----T-----GXXX=---=====--- >> not really possible, but for the sake of it, let's consider it as well
    /// ```
    pub fn start_codon(&self) -> Vec<(u32, u32, Frame)> {
        if !self.is_coding() {
            return vec![];
        }
        let codon = match self.strand {
            Strand::Minus => Codon::upstream(self, &self.cds_end().unwrap()),
            Strand::Plus => Codon::downstream(self, &self.cds_start().unwrap()),
            _ => return vec![],
        };
        if let Ok(res) = codon {
            res.to_tuple()
        } else {
            vec![]
        }
        // if codon.is_ok() {
        //     codon.unwrap().to_tuple()
        // } else {
        //     vec![]
        // }
    }

    /// Returns the stop codons coordinates, split across exons
    ///
    /// The stop codon can be split across multiple exons
    /// and thus the coordinates of it can't be easily calculated
    /// Taking into account the following options, by splitting the
    /// stop codon into multiple exons
    /// ```text
    ///       1....   2....   3....   4....   5....
    ///       12345   12345   12345   12345   12345
    ///    ---=====---===XX---XXXXX---XXXX=---=====---
    /// 1. ---=====---===XX---XXXXX---XXUAG---=====--- >> all in one exon
    /// 2. ---=====---=XXXX---XXXXX---XUAG=---=====--- >> all in one exon
    /// 3. ---=====---==XXX---XXXXX---UAGX=---=====--- >> all in one exon
    /// 4. ---=====---===XX---XXXXU---AGXX=---=====--- >> split
    /// 5. ---=====---====X---XXXUA---GXXX=---=====--- >> split
    /// 6. ---=====---===XU-----A-----G====---=====--- >> not really possible, but for the sake of it, let's consider it as well
    /// ```
    pub fn stop_codon(&self) -> Vec<(u32, u32, Frame)> {
        if !self.is_coding() {
            return vec![];
        }
        let codon = match self.strand {
            Strand::Minus => Codon::downstream(self, &self.cds_start().unwrap()),
            Strand::Plus => Codon::upstream(self, &self.cds_end().unwrap()),
            _ => return vec![],
        };
        if let Ok(res) = codon {
            res.to_tuple()
        } else {
            vec![]
        }
    }
}

impl fmt::Display for Transcript {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "[{}] {} ({}:{}-{})",
            self.gene_symbol,
            self.name,
            self.chrom,
            self.tx_start(),
            self.tx_end()
        )
    }
}

/// Builds a `Transcript`
pub struct TranscriptBuilder<'a> {
    bin: Option<u16>,
    name: Option<&'a str>,
    chrom: Option<&'a str>,
    strand: Strand,
    cds_start_stat: CdsStat,
    cds_end_stat: CdsStat,
    gene_symbol: Option<&'a str>,
    score: Option<f32>,
}

impl<'a> Default for TranscriptBuilder<'a> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> TranscriptBuilder<'a> {
    /// Starts the Builder process to build a `Transcript`
    ///
    /// TranscriptBuilder methods can be chained for easier creation
    ///
    /// # Examples
    ///
    /// ```rust
    /// use atg;
    /// use atg::models::TranscriptBuilder;
    /// let transcript = TranscriptBuilder::new()
    ///     .name("NM_001203247.2")
    ///     .chrom("chr7")
    ///     .gene("EZH2")
    ///     .strand(atg::models::Strand::Minus)
    ///     .build()
    ///     .unwrap();
    /// assert_eq!(transcript.name(), "NM_001203247.2");
    /// ```
    pub fn new() -> Self {
        Self {
            bin: None,
            name: None,
            chrom: None,
            strand: Strand::Unknown,
            cds_start_stat: CdsStat::None,
            cds_end_stat: CdsStat::None,
            gene_symbol: None,
            score: None,
        }
    }

    /// specify the RefGene bin field
    pub fn bin(&mut self, bin: Option<u16>) -> &mut Self {
        self.bin = bin;
        self
    }

    /// specify the name/ID of the transcript, e.g.: NM_001203247.2
    pub fn name(&mut self, name: &'a str) -> &mut Self {
        self.name = Some(name);
        self
    }

    /// specify the chromosome of the transcript
    pub fn chrom(&mut self, chrom: &'a str) -> &mut Self {
        self.chrom = Some(chrom);
        self
    }

    /// specify the gene symbol
    pub fn gene(&mut self, gene_symbol: &'a str) -> &mut Self {
        self.gene_symbol = Some(gene_symbol);
        self
    }

    /// specify the `Strand` of the transcript
    pub fn strand(&mut self, strand: Strand) -> &mut Self {
        self.strand = strand;
        self
    }

    /// specify the `cdsStartStat`
    ///
    /// *Attention:* This field does not neccessarily specify the _start_
    /// of the transcript. It specifies the _genomic leftmost_ element.
    /// This can be the transcript termination sequence, if the transcript
    /// is on the minus strand.
    ///
    /// See and use `cds_start_codon_stat` instead if you want to specify the
    /// actual start-codon.
    pub fn cds_start_stat(&mut self, cds_start_stat: CdsStat) -> &mut Self {
        self.cds_start_stat = cds_start_stat;
        self
    }

    /// specify the `cdsStartStat` or `cdsEndStat`, depending on the
    /// transcript's strand.
    ///
    /// Use this method if you want to specify information about the _start codon_.
    pub fn cds_start_codon_stat(
        &mut self,
        stat: CdsStat,
    ) -> Result<&mut Self, BuildTranscriptError> {
        match self.strand {
            Strand::Plus => Ok(self.cds_start_stat(stat)),
            Strand::Minus => Ok(self.cds_end_stat(stat)),
            _ => Err(BuildTranscriptError::new(
                "Cannot set CDS-Startcodon-Stat without defined strand",
            )),
        }
    }

    /// specify the `cdsEndStat`
    ///
    /// *Attention:* This field does not neccessarily specify the _stop_
    /// of the transcript. It specifies the _genomic rightmost_ element.
    /// This can be the transcript start codon, if the transcript
    /// is on the minus strand.
    ///
    /// See and use `cds_stop_codon_stat` instead if you want to specify the
    /// actual stop-codon.
    pub fn cds_end_stat(&mut self, cds_end_stat: CdsStat) -> &mut Self {
        self.cds_end_stat = cds_end_stat;
        self
    }

    /// specify the `cdsStartStat` or `cdsEndStat`, depending on the
    /// transcript's strand.
    ///
    /// Use this method if you want to specify information about the _stop codon_.
    pub fn cds_stop_codon_stat(
        &mut self,
        stat: CdsStat,
    ) -> Result<&mut Self, BuildTranscriptError> {
        match self.strand {
            Strand::Plus => Ok(self.cds_end_stat(stat)),
            Strand::Minus => Ok(self.cds_start_stat(stat)),
            _ => Err(BuildTranscriptError::new(
                "Cannot set CDS-Startcodon-Stat without defined strand",
            )),
        }
    }

    /// specify the score of the transcript
    pub fn score(&mut self, score: Option<f32>) -> &mut Self {
        self.score = score;
        self
    }

    /// Builds and returns a `Transcript`
    pub fn build(&self) -> Result<Transcript, BuildTranscriptError> {
        let t = Transcript {
            bin: self.bin,
            name: match self.name {
                Some(x) => x.to_string(),
                None => return Err(BuildTranscriptError::new("No name specified")),
            },
            chrom: match self.chrom {
                Some(x) => x.to_string(),
                None => return Err(BuildTranscriptError::new("No chromosome specified")),
            },
            strand: self.strand,
            cds_start_stat: self.cds_start_stat,
            cds_end_stat: self.cds_end_stat,
            exons: vec![],
            gene_symbol: match self.gene_symbol {
                Some(x) => x.to_string(),
                None => return Err(BuildTranscriptError::new("No gene symbol specified")),
            },
            score: self.score,
        };
        Ok(t)
    }
}
