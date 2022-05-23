use core::cmp::{max, min};
use core::str::FromStr;
use std::fmt;

use crate::gtf::constants::*;
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
    gene: String,
    transcript: String,
    exon_number: Option<usize>,
}

impl GtfRecord {
    /// The associated gene symbol
    pub fn gene(&self) -> &str {
        &self.gene
    }

    /// The associated transcript name
    pub fn transcript(&self) -> &str {
        &self.transcript
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
        let mut attributes = format!(
            "gene_id \"{}\"; transcript_id \"{}\"; gene_name \"{}\";",
            self.gene(),
            self.transcript(),
            self.gene()
        );

        if let Some(exon_number) = self.exon_number {
            let additional_attrs = format!(
                " exon_number \"{}\"; exon_id \"{}.{}\";",
                exon_number,
                self.transcript(),
                exon_number
            );
            attributes.push_str(&additional_attrs)
        }
        attributes
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
        let mut columns: Vec<String> = vec![String::new(); MIN_GTF_COLUMNS];

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

        write!(f, "{}", columns.join("\t"))
    }
}

impl FromStr for GtfRecord {
    type Err = ParseGtfError;

    fn from_str(s: &str) -> Result<Self, ParseGtfError> {
        let mut rb = GtfRecordBuilder::new();
        let mut last_idx = 0;

        // Going through the GTF lines column by column
        // manually, since this is faster than splitting
        // the string into columns and iterating
        last_idx = rb.chrom_from_str(s, last_idx)?;
        last_idx = rb.source_from_str(s, last_idx)?;
        last_idx = rb.feature_from_str(s, last_idx)?;
        last_idx = rb.start_from_str(s, last_idx)?;
        last_idx = rb.end_from_str(s, last_idx)?;
        last_idx = rb.score_from_str(s, last_idx)?;
        last_idx = rb.strand_from_str(s, last_idx)?;
        last_idx = rb.frame_from_str(s, last_idx)?;

        // Attributes can be the last or second last column
        // so we check for both cases here
        let (gene, transcript) = match s[last_idx..].find('\t') {
            Some(idx) => parse_attributes(&s[last_idx..idx + last_idx])?,
            None => parse_attributes(s[last_idx..].trim_end())?,
        };
        rb.gene(gene).transcript(transcript);

        match rb.build() {
            Ok(x) => Ok(x),
            Err(err) => Err(ParseGtfError::new(err)),
        }
    }
}

fn parse_attributes(mut attrs: &str) -> Result<(&str, &str), ParseGtfError> {
    let mut gene: &str = "";
    let mut transcript: &str = "";

    while let Some(idx) = attrs.find(';') {
        match parse_attribute(attrs[..idx].trim()) {
            Ok(("gene_id", value)) => gene = value,
            Ok(("transcript_id", value)) => transcript = value,
            Ok((_, _)) => {} // ignore all other attributes
            Err(err) => {
                return Err(ParseGtfError::from_chain(
                    err,
                    &format!(
                        "Unable to parse the attribute column\n\n>>>{}<<<\n\n",
                        attrs
                    ),
                ));
            }
        }

        if !gene.is_empty() && !transcript.is_empty() {
            return Ok((gene, transcript));
        }
        attrs = attrs[(idx + 1)..].trim();
    }
    Err(ParseGtfError::new(format!(
        "Missing gene_id or transcript_id\n>>>{}<<<",
        attrs
    )))
}

fn parse_attribute(attr: &str) -> Result<(&str, &str), ParseGtfError> {
    if let Some(idx) = attr.find(' ') {
        Ok((&attr[..idx], attr[idx + 1..].trim_matches('\"')))
    } else {
        Err(ParseGtfError {
            message: format!(
                "Unable to parse the attribute\n\n{}\nPlease check your GTF input.",
                attr
            ),
        })
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
/// use atg::gtf::{GtfFeature, GtfRecordBuilder};
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
///     .gene("Foo")
///     .transcript("Bar")
///     .build()
///     .unwrap();
///
/// assert_eq!(record.chrom(), "chr1");
/// ```
pub struct GtfRecordBuilder<'a> {
    gene: Option<&'a str>,
    transcript: Option<&'a str>,
    chrom: Option<&'a str>,
    source: Option<&'a str>,
    feature: Option<GtfFeature>,
    start: Option<u32>,
    end: Option<u32>,
    score: Option<f32>,
    strand: Strand,
    frame_offset: Frame,
    exon_number: Option<&'a usize>,
}

impl<'a> Default for GtfRecordBuilder<'a> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'a> GtfRecordBuilder<'a> {
    pub fn new() -> Self {
        Self {
            gene: None,
            transcript: None,
            chrom: None,
            source: None,
            feature: None,
            start: None,
            end: None,
            score: None,
            strand: Strand::Unknown,
            frame_offset: Frame::None,
            exon_number: None,
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

    pub fn gene(&mut self, gene: &'a str) -> &mut Self {
        self.gene = Some(gene);
        self
    }

    pub fn transcript(&mut self, transcript: &'a str) -> &mut Self {
        self.transcript = Some(transcript);
        self
    }

    pub fn exon_number(&mut self, exon_number: &'a usize) -> &mut Self {
        self.exon_number = Some(exon_number);
        self
    }

    /// Uses the next tab-separated substring as chrom
    /// and returns the start-index of the next substring
    fn chrom_from_str(&mut self, s: &'a str, start_idx: usize) -> Result<usize, ParseGtfError> {
        match s[start_idx..].find('\t') {
            Some(idx) => {
                let end_idx = idx + start_idx;
                self.chrom(&s[start_idx..end_idx]);
                Ok(end_idx + 1)
            }
            None => Err(ParseGtfError::new(format!(
                "too few columns in line: {}",
                s
            ))),
        }
    }

    /// Uses the next tab-separated substring as source
    /// and returns the start-index of the next substring
    fn source_from_str(&mut self, s: &'a str, start_idx: usize) -> Result<usize, ParseGtfError> {
        match s[start_idx..].find('\t') {
            Some(idx) => {
                let end_idx = idx + start_idx;
                self.source(&s[start_idx..end_idx]);
                Ok(end_idx + 1)
            }
            None => Err(ParseGtfError::new(format!(
                "too few columns in line: {}",
                s
            ))),
        }
    }

    /// Uses the next tab-separated substring as feature
    /// and returns the start-index of the next substring
    fn feature_from_str(&mut self, s: &'a str, start_idx: usize) -> Result<usize, ParseGtfError> {
        match s[start_idx..].find('\t') {
            Some(idx) => {
                let end_idx = idx + start_idx;
                let feature = GtfFeature::from_str(&s[start_idx..end_idx])?;
                self.feature(feature);
                Ok(end_idx + 1)
            }
            None => Err(ParseGtfError::new(format!(
                "too few columns in line: {}",
                s
            ))),
        }
    }

    /// Uses the next tab-separated substring as start
    /// and returns the start-index of the next substring
    fn start_from_str(&mut self, s: &'a str, start_idx: usize) -> Result<usize, ParseGtfError> {
        match s[start_idx..].find('\t') {
            Some(idx) => {
                let end_idx = idx + start_idx;
                let start = s[start_idx..end_idx].parse::<u32>()?;
                self.start(start);
                Ok(end_idx + 1)
            }
            None => Err(ParseGtfError::new(format!(
                "too few columns in line: {}",
                s
            ))),
        }
    }

    /// Uses the next tab-separated substring as end
    /// and returns the start-index of the next substring
    fn end_from_str(&mut self, s: &'a str, start_idx: usize) -> Result<usize, ParseGtfError> {
        match s[start_idx..].find('\t') {
            Some(idx) => {
                let end_idx = idx + start_idx;
                let end = s[start_idx..end_idx].parse::<u32>()?;
                self.end(end);
                Ok(end_idx + 1)
            }
            None => Err(ParseGtfError::new(format!(
                "too few columns in line: {}",
                s
            ))),
        }
    }

    /// Uses the next tab-separated substring as score
    /// and returns the start-index of the next substring
    fn score_from_str(&mut self, s: &'a str, start_idx: usize) -> Result<usize, ParseGtfError> {
        match s[start_idx..].find('\t') {
            Some(idx) => {
                let end_idx = idx + start_idx;
                let score = s[start_idx..end_idx].parse::<f32>().ok();
                self.score_option(score);
                Ok(end_idx + 1)
            }
            None => Err(ParseGtfError::new(format!(
                "too few columns in line: {}",
                s
            ))),
        }
    }

    /// Uses the next tab-separated substring as strand
    /// and returns the start-index of the next substring
    fn strand_from_str(&mut self, s: &'a str, start_idx: usize) -> Result<usize, ParseGtfError> {
        match s[start_idx..].find('\t') {
            Some(idx) => {
                let end_idx = idx + start_idx;
                let strand = Strand::from_str(&s[start_idx..end_idx])?;
                self.strand(strand);
                Ok(end_idx + 1)
            }
            None => Err(ParseGtfError::new(format!(
                "too few columns in line: {}",
                s
            ))),
        }
    }

    /// Uses the next tab-separated substring as frame
    /// and returns the start-index of the next substring
    fn frame_from_str(&mut self, s: &'a str, start_idx: usize) -> Result<usize, ParseGtfError> {
        match s[start_idx..].find('\t') {
            Some(idx) => {
                let end_idx = idx + start_idx;
                let frame = Frame::from_str(&s[start_idx..end_idx])?;
                self.frame_offset(frame);
                Ok(end_idx + 1)
            }
            None => Err(ParseGtfError::new(format!(
                "too few columns in line: {}",
                s
            ))),
        }
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
            gene: match self.gene {
                Some(x) => x.to_string(),
                None => return Err("Missing gene".to_string()),
            },
            transcript: match self.transcript {
                Some(x) => x.to_string(),
                None => return Err("Missing transcript".to_string()),
            },
            exon_number: self.exon_number.copied(),
        };
        Ok(r)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_attributes_parsing() {
        let col = "gene_id \"ZBTB16\"; transcript_id \"NM_001354751.2\"; exon_number \"2\"; exon_id \"NM_001354751.2.2\"; gene_name \"ZBTB16\";";
        let attr = parse_attributes(col).unwrap();
        assert_eq!(attr.0, "ZBTB16");
        assert_eq!(attr.1, "NM_001354751.2");

        let col = "gene_id \"ZBTB16\"; transcript_id \"NM_001354752.1\"; gene_name \"ZBTB16\";";
        let attr = parse_attributes(col).unwrap();
        assert_eq!(attr.0, "ZBTB16");
        assert_eq!(attr.1, "NM_001354752.1");
    }
    // let line =  "chr11\tncbiRefSeq.2021-05-17\texon\t113933933\t113935290\t.\t+\t.\tgene_id \"ZBTB16\"; transcript_id \"NM_001354751.2\"; exon_number \"2\"; exon_id \"NM_001354751.2.2\"; gene_name \"ZBTB16\";"

    #[test]
    fn test_single_attribute_parsing() {
        let res = parse_attribute("gene_id \"ZBTB16\"").unwrap();
        assert_eq!(res.0, "gene_id");
        assert_eq!(res.1, "ZBTB16");

        let res = parse_attribute("transcript_id \"NM_001354751.2\"").unwrap();
        assert_eq!(res.0, "transcript_id");
        assert_eq!(res.1, "NM_001354751.2");
    }
}
