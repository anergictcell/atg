use std::fmt;

use crate::models::{Exon, Strand, Transcript};

#[derive(Clone, Debug, PartialEq)]
pub struct BedLine {
    chrom: String,
    start: u32,
    end: u32,
    name: Option<String>,
    score: Option<f32>,
    strand: Option<Strand>,
    thick_start: Option<u32>,
    thick_end: Option<u32>,
    item_rgb: Option<String>,
    block_count: Option<usize>,
    block_sizes: Option<Vec<u32>>,
    block_starts: Option<Vec<u32>>,
}

impl BedLine {
    /// Create one row of a bed file for a single exon
    ///
    /// Please note that this is different from converting from
    /// a `Transcript`, which will contain all Exons
    pub fn from_exon(exon: Exon, chrom: &str) -> Self {
        Self {
            chrom: chrom.to_owned(),
            start: exon.start(),
            end: exon.end(),
            thick_start: *exon.cds_start(),
            thick_end: *exon.cds_end(),
            name: None,
            score: None,
            strand: None,
            item_rgb: None,
            block_count: None,
            block_sizes: None,
            block_starts: None,
        }
    }

    pub fn name_mut(&mut self) -> &mut Option<String> {
        &mut self.name
    }
}

impl From<&Transcript> for BedLine {
    /// Note: Coordinates in bed format are 0-based, so we must
    /// subtract 1 from each coordinate
    fn from(transcript: &Transcript) -> Self {
        Self {
            chrom: transcript.chrom().to_owned(),
            start: transcript.tx_start() - 1,
            end: transcript.tx_end(),
            thick_start: transcript.cds_start().map(|s| s - 1),
            thick_end: transcript.cds_end(),
            name: Some(format!("{}:{}", transcript.gene(), transcript.name())),
            score: transcript.score(),
            strand: Some(transcript.strand()),
            item_rgb: None,
            block_count: Some(transcript.exon_count()),
            block_sizes: Some(transcript.exons().iter().map(|exon| exon.len()).collect()),
            block_starts: Some(
                transcript
                    .exons()
                    .iter()
                    .map(|exon| exon.start() - transcript.tx_start())
                    .collect(),
            ),
        }
    }
}

impl fmt::Display for BedLine {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.chrom,
            self.start,
            self.end,
            match &self.name {
                Some(x) => x,
                _ => "",
            },
            match self.score {
                Some(x) => x.to_string(),
                _ => "".to_string(),
            },
            match self.strand {
                Some(x) => x.to_string(),
                _ => "".to_string(),
            },
            match self.thick_start {
                Some(x) => x.to_string(),
                _ => "".to_string(),
            },
            match self.thick_end {
                Some(x) => x.to_string(),
                _ => "".to_string(),
            },
            match &self.item_rgb {
                Some(x) => x,
                _ => "212,16,48", // Centogene red :)
            },
            match self.block_count {
                Some(x) => x.to_string(),
                _ => "".to_string(),
            },
            match &self.block_sizes {
                Some(x) => x
                    .iter()
                    .map(|i| i.to_string())
                    .collect::<Vec<String>>()
                    .join(","),
                _ => "".to_string(),
            },
            match &self.block_starts {
                Some(x) => x
                    .iter()
                    .map(|i| i.to_string())
                    .collect::<Vec<String>>()
                    .join(","),
                _ => "".to_string(),
            },
        )
    }
}
