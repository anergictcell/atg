use core::cmp::{max, min};
use core::str::FromStr;
use std::fmt;

use crate::gtf::constants::*;
use crate::gtf::utils::Attributes;
use crate::models;
use crate::models::{Exon, Frame, Strand};
use crate::utils::errors::ParseGtfError;

/// Describes the actual feature type (Exon, CDS, etc) of [`GtfRecord`]
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GtfFeature {
    Exon,
    CDS,
    StartCodon,
    StopCodon,
    UTR,
    UTR5,
    UTR3,
    Inter,
    InterCNS,
    IntronCNS,
    Gene,
    Transcript,
    Selenocysteine,
}

impl FromStr for GtfFeature {
    type Err = ParseGtfError;
    fn from_str(s: &str) -> Result<Self, ParseGtfError> {
        match s {
            "CDS" => Ok(Self::CDS),
            "start_codon" => Ok(Self::StartCodon),
            "stop_codon" => Ok(Self::StopCodon),
            "5UTR" => Ok(Self::UTR5),
            "3UTR" => Ok(Self::UTR3),
            "UTR" => Ok(Self::UTR),
            "inter" => Ok(Self::Inter),
            "inter_CNS" => Ok(Self::InterCNS),
            "intron_CNS" => Ok(Self::IntronCNS),
            "exon" => Ok(Self::Exon),
            "gene" => Ok(Self::Gene),
            "transcript" => Ok(Self::Transcript),
            "Selenocysteine" => Ok(Self::Selenocysteine),
            _ => Err(ParseGtfError {
                message: format!("invalid feature type {}", s),
            }),
        }
    }
}

impl fmt::Display for GtfFeature {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::CDS => "CDS",
                Self::StartCodon => "start_codon",
                Self::StopCodon => "stop_codon",
                Self::UTR5 => "5UTR",
                Self::UTR3 => "3UTR",
                Self::UTR => "UTR",
                Self::Inter => "inter",
                Self::InterCNS => "inter_CNS",
                Self::IntronCNS => "intron_CNS",
                Self::Exon => "exon",
                Self::Gene => "gene",
                Self::Transcript => "transcript",
                Self::Selenocysteine => "Selenocysteine",
            }
        )
    }
}

/// Represents a single line of a GTF file (or other input)
///
/// One record *does not* equal a transcript, but only one subset feature
/// e.g.: Exon, Start-Codon etc
///
/// Use [`GtfRecordBuilder`] to create [`GtfRecord`]
#[derive(Debug, PartialEq)]
pub struct GtfRecord {
    chrom: String,
    source: String,
    feature: GtfFeature,
    start: u32,
    end: u32,
    score: Option<f32>,
    strand: Strand,
    frame_offset: Frame,
    attributes: Attributes,
    comments: Option<String>,
}

impl GtfRecord {
    /// The associated gene symbol
    pub fn gene(&self) -> &str {
        self.attributes.gene()
    }

    /// The associated transcript name
    pub fn transcript(&self) -> &str {
        self.attributes.transcript()
    }

    /// The genomic seqname, in most cases the chromosome
    pub fn chrom(&self) -> &str {
        &self.chrom
    }

    /// The [`Strand`](crate::models::Strand) of the transcript
    pub fn strand(&self) -> &Strand {
        &self.strand
    }

    /// The confidence score of the annotation
    pub fn score(&self) -> &Option<f32> {
        &self.score
    }

    /// The type of [feature](GtfFeature) (CDS, Exon, Start-Codon, etc)
    pub fn feature(&self) -> &GtfFeature {
        &self.feature
    }

    /// the start (genomic left-most) position of the Record
    pub fn start(&self) -> u32 {
        self.start
    }

    /// the end (genomic right-most) position of the Record
    pub fn end(&self) -> u32 {
        self.end
    }

    fn attributes_to_string(&self) -> String {
        let mut res: Vec<String> = vec![];
        for attr in &self.attributes.all() {
            res.push(format!("{} \"{}\";", attr.0, attr.1))
        }
        res.join(" ")
    }

    /// Modifies an [`Exon`](crate::models::Exon) to include the [`GtfRecord`]
    ///
    /// use this method if you want to update an [`Exon`](crate::models::Exon)
    /// with the data of the [`GtfRecord`]
    pub fn add_to_exon(self, mut exon: Exon) -> Exon {
        *exon.start_mut() = min(exon.start(), self.start);
        *exon.end_mut() = max(exon.end(), self.end);

        match self.feature {
            GtfFeature::CDS => {
                *exon.cds_start_mut() =
                    Some(min(exon.cds_start().unwrap_or(self.start), self.start));
                *exon.cds_end_mut() = Some(max(exon.cds_end().unwrap_or(self.end), self.end));
                exon.set_frame(self.frame_offset);
            }
            GtfFeature::StopCodon => match self.strand {
                Strand::Plus => {
                    *exon.cds_start_mut() =
                        Some(min(self.start, exon.cds_start().unwrap_or(self.start)));
                    *exon.cds_end_mut() = Some(self.end);
                }
                Strand::Minus => {
                    *exon.cds_start_mut() = Some(self.start);
                    *exon.cds_end_mut() = Some(max(self.end, exon.cds_end().unwrap_or(self.end)));
                }
                _ => {}
            },
            GtfFeature::StartCodon => match self.strand {
                Strand::Plus => {
                    *exon.cds_start_mut() = Some(self.start);
                    *exon.cds_end_mut() = Some(max(self.end, exon.cds_end().unwrap_or(self.end)));
                }
                Strand::Minus => {
                    *exon.cds_start_mut() =
                        Some(min(self.start, exon.cds_start().unwrap_or(self.start)));
                    *exon.cds_end_mut() = Some(self.end);
                }
                _ => {}
            },
            _ => {}
        }

        if exon.is_coding() && !exon.frame_offset().is_known() {
            exon.set_frame(self.frame_offset);
        }

        exon
    }
}

impl fmt::Display for GtfRecord {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut columns: Vec<String> = vec![
            "".to_string();
            match self.comments {
                Some(_) => MAX_GTF_COLUMNS,
                None => MIN_GTF_COLUMNS,
            }
        ];

        columns[CHROMOSOME_COL] = self.chrom.to_string();
        columns[SOURCE_COL] = self.source.to_string();
        columns[FEATURE_COL] = self.feature.to_string();
        columns[START_COL] = self.start.to_string();
        columns[END_COL] = self.end.to_string();
        columns[SCORE_COL] = match self.score {
            Some(x) => x.to_string(),
            _ => ".".to_string(),
        };
        columns[STRAND_COL] = self.strand.to_string();
        columns[FRAME_COL] = self.frame_offset.to_string();
        columns[ATTRIBUTES_COL] = self.attributes_to_string();
        if let Some(x) = &self.comments {
            columns[COMMENTS_COL] = x.to_string();
        }

        write!(f, "{}", columns.join("\t"))
    }
}

impl FromStr for GtfRecord {
    type Err = ParseGtfError;

    fn from_str(s: &str) -> Result<Self, ParseGtfError> {
        let cols: Vec<&str> = s.split('\t').collect();
        if cols.len() < MIN_GTF_COLUMNS || cols.len() > MAX_GTF_COLUMNS {
            return Err(ParseGtfError {
                message: format!(
                    "Wrong number of columns in GTF line. Expected {}-{}, received {}",
                    MIN_GTF_COLUMNS,
                    MAX_GTF_COLUMNS,
                    cols.len()
                ),
            });
        };

        let res = GtfRecord {
            chrom: models::parse_chrom(cols[CHROMOSOME_COL]),
            source: cols[SOURCE_COL].to_string(),
            feature: match GtfFeature::from_str(cols[FEATURE_COL]) {
                Ok(x) => x,
                Err(x) => {
                    return Err(ParseGtfError {
                        message: x.to_string(),
                    })
                }
            },
            start: match cols[START_COL].parse::<u32>() {
                Ok(x) => x,
                Err(x) => {
                    return Err(ParseGtfError {
                        message: x.to_string(),
                    })
                }
            },
            end: match cols[END_COL].parse::<u32>() {
                Ok(x) => x,
                Err(x) => {
                    return Err(ParseGtfError {
                        message: x.to_string(),
                    })
                }
            },
            score: match cols[SCORE_COL].parse::<f32>() {
                Ok(x) => Some(x),
                _ => None,
            },
            strand: match Strand::from_str(cols[STRAND_COL]) {
                Ok(x) => x,
                Err(x) => return Err(ParseGtfError { message: x }),
            },
            frame_offset: match Frame::from_str(cols[FRAME_COL]) {
                Ok(x) => x,
                Err(x) => return Err(ParseGtfError { message: x }),
            },
            attributes: match Attributes::from_str(cols[ATTRIBUTES_COL]) {
                Ok(x) => x,
                Err(err) => {
                    return Err(ParseGtfError::from_chain(
                        err,
                        &format!("Error in line:\n{}\n", s),
                    ))
                }
            },
            comments: match cols.len() {
                MAX_GTF_COLUMNS => Some(cols[COMMENTS_COL].to_string()),
                _ => None,
            },
        };
        Ok(res)
    }
}

impl From<GtfRecord> for models::Exon {
    fn from(feature: GtfRecord) -> Self {
        let mut exon = Exon::new(feature.start, feature.end, None, None, feature.frame_offset);
        match feature.feature {
            GtfFeature::CDS => {
                *exon.cds_start_mut() = Some(feature.start);
                *exon.cds_end_mut() = Some(feature.end);
            }
            GtfFeature::UTR | GtfFeature::UTR3 | GtfFeature::UTR5 => {
                *exon.frame_offset_mut() = Frame::None;
            }
            GtfFeature::StopCodon => {
                *exon.cds_start_mut() = Some(feature.start);
                *exon.cds_end_mut() = Some(feature.end);
            }
            GtfFeature::StartCodon => {
                *exon.cds_start_mut() = Some(feature.start);
                *exon.cds_end_mut() = Some(feature.end);
            }
            _ => {}
        }
        exon
    }
}

/// A builder for [`GtfRecord`]s, providing a cleaner API
///
/// # Examples
///
/// ```rust
/// use atg::gtf::{Attributes, GtfFeature, GtfRecordBuilder};
/// use atg::models::{Frame, Strand};
/// use std::str::FromStr;
///
/// let record = GtfRecordBuilder::new()
///     .chrom("chr1")
///     .source("local source")
///     .feature(GtfFeature::CDS)
///     .start(1)
///     .end(3)
///     .score(1.4)
///     // alternative:
///     // score_option(Some(1.4))
///     .strand(Strand::Plus)
///     .frame_offset(Frame::Zero)
///     .attributes(Attributes::from_str("gene_id Foo; transcript_id Bar").unwrap())
///     .build()
///     .unwrap();
///
/// assert_eq!(record.chrom(), "chr1");
/// ```
pub struct GtfRecordBuilder<'a> {
    chrom: Option<&'a str>,
    source: Option<&'a str>,
    feature: Option<GtfFeature>,
    start: Option<u32>,
    end: Option<u32>,
    score: Option<f32>,
    strand: Strand,
    frame_offset: Frame,
    attributes: Option<Attributes>,
    comments: Option<String>,
}

impl<'a> Default for GtfRecordBuilder<'a> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> GtfRecordBuilder<'a> {
    pub fn new() -> Self {
        Self {
            chrom: None,
            source: None,
            feature: None,
            start: None,
            end: None,
            score: None,
            strand: Strand::Unknown,
            frame_offset: Frame::Dot,
            attributes: None,
            comments: None,
        }
    }

    pub fn chrom(&mut self, chrom: &'a str) -> &mut Self {
        self.chrom = Some(chrom);
        self
    }
    pub fn source(&mut self, source: &'a str) -> &mut Self {
        self.source = Some(source);
        self
    }
    pub fn feature(&mut self, feature: GtfFeature) -> &mut Self {
        self.feature = Some(feature);
        self
    }
    pub fn start(&mut self, start: u32) -> &mut Self {
        self.start = Some(start);
        self
    }
    pub fn end(&mut self, end: u32) -> &mut Self {
        self.end = Some(end);
        self
    }

    /// Adds Some(score) to the GtfRecord
    ///
    /// Use this method if you are certain that you do have an actual
    /// score value. If you are unsure, use [`score_option`](GtfRecordBuilder::score_option) instead.
    pub fn score(&mut self, score: f32) -> &mut Self {
        self.score = Some(score);
        self
    }

    /// Adds an Option to the score of the GtfRecord
    ///
    /// Use this method if you are unsure if you have an actual
    /// score value or `None` in your input data
    pub fn score_option(&mut self, score: Option<f32>) -> &mut Self {
        self.score = score;
        self
    }
    pub fn strand(&mut self, strand: Strand) -> &mut Self {
        self.strand = strand;
        self
    }
    pub fn frame_offset(&mut self, frame_offset: Frame) -> &mut Self {
        self.frame_offset = frame_offset;
        self
    }
    pub fn attributes(&mut self, attributes: Attributes) -> &mut Self {
        self.attributes = Some(attributes);
        self
    }
    pub fn comments(&mut self, comments: &'a str) -> &mut Self {
        self.comments = Some(comments.to_string());
        self
    }

    /// Builds and returns a [`Transcript`](crate::models::Transcript)
    pub fn build(&mut self) -> Result<GtfRecord, String> {
        let r = GtfRecord {
            chrom: match self.chrom {
                Some(x) => x.to_string(),
                None => return Err("Missing chrom".to_string()),
            },
            source: match self.source {
                Some(x) => x.to_string(),
                None => return Err("Missing source".to_string()),
            },
            feature: match self.feature {
                Some(x) => x,
                None => return Err("Missing feature".to_string()),
            },
            start: match self.start {
                Some(x) => x,
                None => return Err("Missing start".to_string()),
            },
            end: match self.end {
                Some(x) => x,
                None => return Err("Missing end".to_string()),
            },
            score: self.score,
            strand: self.strand,
            frame_offset: self.frame_offset,
            attributes: match self.attributes.take() {
                Some(x) => x,
                None => return Err("Missing attributes".to_string()),
            },
            comments: self.comments.take(),
        };
        Ok(r)
    }
}
