use core::convert::TryFrom;
use std::fmt;

use crate::gtf::{GtfFeature, GtfRecord};
use crate::models::{CdsStat, Exon, Strand, Transcript, TranscriptBuilder};
use crate::utils::errors::ParseGtfError;

/// Groups all [`GtfRecord`] object that belong to one Transcript
///
/// Most transcripts are composed of several [`GtfRecord`]s.
/// [`GtfRecordsGroup`] handles the grouping and allows
/// conversions into [`Transcript`](Transcript).
pub struct GtfRecordsGroup {
    transcript: String,
    exons: Vec<GtfRecord>,
    sorted: bool,
}

impl GtfRecordsGroup {
    pub fn new(transcript_id: &str) -> Self {
        Self {
            transcript: transcript_id.to_string(),
            exons: vec![],
            sorted: false,
        }
    }

    pub fn add_exon(&mut self, exon: GtfRecord) {
        self.exons.push(exon)
    }

    pub fn transcript(&self) -> &str {
        &self.transcript
    }

    pub fn gene(&self) -> &str {
        self.exons[0].gene()
    }

    pub fn chrom(&self) -> &str {
        self.exons[0].chrom()
    }

    pub fn strand(&self) -> &Strand {
        self.exons[0].strand()
    }

    pub fn score(&self) -> &Option<f32> {
        self.exons[0].score()
    }

    fn prepare(&mut self) {
        self.exons.sort_unstable_by_key(|x| x.start());
        self.exons.reverse();
        self.sorted = true;
    }

    fn next_exon(&mut self) -> Option<Exon> {
        if self.exons.is_empty() {
            return None;
        }
        let mut exon = Exon::from(self.exons.pop().unwrap()); // cannot fail, we test for emptyness of the exon vec

        loop {
            let line = self.exons.pop();
            if let Some(x) = line {
                // TODO: Change this to only merge book-ended
                // features if the two features are CDS + StopCodon.
                // All other features should be kept separate, even if they
                // are right next to each other.
                // Genetically speaking, they would be considered as one exon
                // and not two separate exons. But to keep the input data intact,
                // this should be changed.
                if x.start() <= (exon.end() + 1) {
                    exon = x.add_to_exon(exon);
                } else {
                    self.exons.push(x);
                    break;
                }
            } else {
                break;
            }
        }
        Some(exon)
    }

    /// Returns all exons of the transcript as `Vector`
    ///
    /// All rows of the GFT file are grouped by genomic location
    pub fn exons(&mut self) -> Vec<Exon> {
        if !self.sorted {
            self.prepare();
        }
        let mut exons: Vec<Exon> = vec![];
        while let Some(x) = self.next_exon() {
            exons.push(x)
        }

        exons
    }

    fn cds_start_stat(&self) -> CdsStat {
        self.cds_stat(GtfFeature::StartCodon)
    }

    fn cds_end_stat(&self) -> CdsStat {
        self.cds_stat(GtfFeature::StopCodon)
    }

    fn cds_stat(&self, start_stop: GtfFeature) -> CdsStat {
        let mut cds_present = false;
        for exon in &self.exons {
            if exon.feature() == &start_stop {
                return CdsStat::Complete;
            }
            if exon.feature() == &GtfFeature::CDS {
                cds_present = true;
            }
        }
        if cds_present {
            CdsStat::Incomplete
        } else {
            CdsStat::Unknown
        }
    }
}

impl Default for GtfRecordsGroup {
    fn default() -> Self {
        Self::new("Unknown")
    }
}

impl fmt::Display for GtfRecordsGroup {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "GTFT {} with {} exons",
            self.transcript,
            self.exons.len()
        )
    }
}

impl TryFrom<GtfRecordsGroup> for Transcript {
    type Error = ParseGtfError;
    /// Returns a `Transcript` based on features of the GTF file,
    /// belonging to one transcript
    fn try_from(mut gtf_transcript: GtfRecordsGroup) -> Result<Self, ParseGtfError> {
        if gtf_transcript.exons.is_empty() {
            return Err(ParseGtfError {
                message: format!("No exons in {}", gtf_transcript),
            });
        }
        let transcript = TranscriptBuilder::new()
            .name(gtf_transcript.transcript())
            .gene(gtf_transcript.gene())
            .chrom(gtf_transcript.chrom())
            .strand(*gtf_transcript.strand())
            .cds_start_codon_stat(gtf_transcript.cds_start_stat())?
            .cds_stop_codon_stat(gtf_transcript.cds_end_stat())?
            .score(*gtf_transcript.score())
            .build();
        match transcript {
            Ok(mut x) => {
                x.append_exons(&mut gtf_transcript.exons());
                Ok(x)
            }
            _ => Err(ParseGtfError {
                message: "Unable to build Transcript".to_string(),
            }),
        }
    }
}
