use core::convert::TryFrom;
use std::fmt;

use crate::gtf::{GtfFeature, GtfRecord};
use crate::utils::errors::ParseGtfError;
use crate::models::{CdsStat, Exon, Transcript, TranscriptBuilder};

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

    pub fn exons(&self) -> &Vec<GtfRecord> {
        &self.exons
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
        let mut exon = Exon::from(self.exons.pop().unwrap());

        loop {
            let line = self.exons.pop();
            if let Some(x) = line {
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

    /// Returns all exons of the transcript as ```Vector```
    ///
    /// All rows of the GFT file are grouped by genomic location
    pub fn to_exons(&mut self) -> Vec<Exon> {
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
    /// Returns a ```Transcript``` based on features of the GTF file,
    /// belonging to one transcript
    fn try_from(mut gtf_transcript: GtfRecordsGroup) -> Result<Self, ParseGtfError> {
        if gtf_transcript.exons.is_empty() {
            return Err(ParseGtfError {
                message: format!("No exons in {}", gtf_transcript),
            });
        }
        let transcript = TranscriptBuilder::new()
            .name(gtf_transcript.exons()[0].transcript())
            .gene(gtf_transcript.exons()[0].gene())
            .chrom(gtf_transcript.exons()[0].chrom())
            .strand(*gtf_transcript.exons()[0].strand())
            .cds_start_codon_stat(gtf_transcript.cds_start_stat())?
            .cds_stop_codon_stat(gtf_transcript.cds_end_stat())?
            .score(*gtf_transcript.exons()[0].score())
            .build();
        match transcript {
            Ok(mut x) => {
                x.append_exons(&mut gtf_transcript.to_exons());
                Ok(x)
            }
            _ => Err(ParseGtfError {
                message: "Unable to build Transcript".to_string(),
            }),
        }
    }
}
