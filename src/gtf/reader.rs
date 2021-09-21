use std::collections::HashMap;
use std::convert::TryFrom;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::str::FromStr;

use crate::gtf::{GtfFeature, GtfRecord, GtfRecordsGroup, ParseGtfError};
use crate::models::{Transcript, TranscriptRead, Transcripts};
use crate::utils::errors::ReadWriteError;

/// Reads a GTF file
pub struct Reader<R> {
    inner: std::io::BufReader<R>,
}

impl Reader<File> {
    /// Creates a Reader instance that reads from a File
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, ReadWriteError> {
        match File::open(path.as_ref()) {
            Ok(file) => Ok(Self::new(file)),
            Err(err) => Err(ReadWriteError::new(err)),
        }
    }
}

impl<R: std::io::Read> Reader<R> {
    pub fn new(reader: R) -> Self {
        Reader {
            inner: BufReader::new(reader),
        }
    }

    pub fn with_capacity(capacity: usize, reader: R) -> Self {
        Reader {
            inner: BufReader::with_capacity(capacity, reader),
        }
    }

    /// Returns one line of a GTF file as `GtfRecord`
    ///
    /// This method should rarely be used. GTF files can contain unordered
    /// records and handling lines individually is rarely desired.
    fn line(&mut self) -> Option<Result<GtfRecord, ParseGtfError>> {
        let mut line = String::new();
        match self.inner.read_line(&mut line) {
            Ok(_) => {}
            Err(x) => {
                return Some(Err(ParseGtfError {
                    message: x.to_string(),
                }))
            }
        }

        if line.starts_with('#') {
            return self.line();
        }

        if line.is_empty() {
            None
        } else {
            Some(GtfRecord::from_str(line.trim_end()))
        }
    }
}

impl<R: std::io::Read> TranscriptRead for Reader<R> {
    /// Reads in GTF data and returns the final list of `Transcripts`
    fn transcripts(&mut self) -> Result<Transcripts, ReadWriteError> {
        let mut transcript_hashmap: HashMap<String, GtfRecordsGroup> = HashMap::new();
        while let Some(line) = self.line() {
            let gtf_record = match line {
                Err(x) => return Err(ReadWriteError::from(x)),
                Ok(line) => line,
            };
            let key = &gtf_record.transcript().to_string();
            if transcript_hashmap.get_mut(key).is_none() {
                transcript_hashmap.insert(key.to_string(), GtfRecordsGroup::new(key));
            }
            let transcript = transcript_hashmap.get_mut(key).unwrap();
            match gtf_record.feature() {
                GtfFeature::Exon => transcript.add_exon(gtf_record),
                GtfFeature::CDS => transcript.add_exon(gtf_record),
                GtfFeature::UTR | GtfFeature::UTR3 | GtfFeature::UTR5 => {
                    transcript.add_exon(gtf_record)
                }
                GtfFeature::StartCodon => {
                    transcript.add_exon(gtf_record);
                }
                GtfFeature::StopCodon => transcript.add_exon(gtf_record),
                _ => {}
            }
        }

        let mut res = Transcripts::with_capacity(transcript_hashmap.len());
        for (_, gtf_transcript) in transcript_hashmap.drain() {
            match Transcript::try_from(gtf_transcript) {
                Ok(transcript) => res.push(transcript),
                Err(err) => {
                    return Err(ReadWriteError::new(&format!("Error parsing {}", err)));
                }
            }
        }
        Ok(res)
    }
}

#[cfg(test)]
mod test_nm_001385228 {
    use super::*;

    #[test]
    fn test_read() {
        let mut reader = Reader::from_file("tests/data/NM_001385228.1_2.gtf").unwrap();
        let tr = match reader.transcripts() {
            Ok(res) => res,
            _ => panic!("No transcripts could be read"),
        };

        assert_eq!(tr.len(), 1);
        let t = tr.by_name("NM_001385228.1_2").unwrap();
        assert_eq!(t.exons().len(), 9);
        assert_eq!(t.cds_start().unwrap(), 206105119);
        assert_eq!(t.cds_end().unwrap(), 206135359);

        let start = t.start_codon();
        assert_eq!(start.len(), 1);

        let stop = t.stop_codon();
        assert_eq!(stop.len(), 2);

        assert_eq!(t.exons()[0].is_coding(), false);
        assert_eq!(t.exons()[1].is_coding(), false);
        assert_eq!(t.exons()[2].is_coding(), false);
        assert_eq!(t.exons()[3].is_coding(), false);
        assert_eq!(t.exons()[4].is_coding(), false);
        assert_eq!(t.exons()[5].is_coding(), true);
        assert_eq!(t.exons()[6].is_coding(), true);
        assert_eq!(t.exons()[7].is_coding(), true);
        assert_eq!(t.exons()[8].is_coding(), false);

        assert_eq!(t.exons()[5].cds_start().unwrap(), 206105119);
        assert_eq!(t.exons()[6].cds_start().unwrap(), 206105123);
        assert_eq!(t.exons()[7].cds_start().unwrap(), 206135293);

        assert_eq!(t.exons()[5].cds_end().unwrap(), 206105120);
        assert_eq!(t.exons()[6].cds_end().unwrap(), 206105206);
        assert_eq!(t.exons()[7].cds_end().unwrap(), 206135359);
    }
}

#[cfg(test)]
mod test_nm_201550 {
    use super::*;

    #[test]
    fn test_read() {
        let mut reader = Reader::from_file("tests/data/NM_201550.4.gtf")
            .expect("Something failed reading the GTF file NM_201550.4.gtf");
        let tr = match reader.transcripts() {
            Ok(res) => res,
            _ => panic!("No transcripts could be read"),
        };
        assert_eq!(tr.len(), 1);
        let t = tr.by_name("NM_201550.4").unwrap();
        assert_eq!(t.exons().len(), 1);
        assert_eq!(t.cds_start().unwrap(), 70003785);
        assert_eq!(t.cds_end().unwrap(), 70004618);

        let start = t.start_codon();
        assert_eq!(start.len(), 1);

        let stop = t.stop_codon();
        assert_eq!(stop.len(), 1);

        assert_eq!(t.exons()[0].is_coding(), true);
    }
}

/*
This one fails due to book-ended exons in the input file
#[cfg(test)]
mod test_nm_001371720 {
    use super::*;

    #[test]
    fn test_read() {
        let mut reader = Reader::from_file("tests/data/NM_001371720.1.gtf")
            .expect("Something failed reading the GTF file NM_001371720.1.gtf");
        let tr = match reader.transcripts() {
            Ok(res) => res,
            _ => panic!("No transcripts could be read")
        };
        assert_eq!(tr.len(), 1);
        let t = tr.by_name("NM_001371720.1").unwrap();
        assert_eq!(t.exons().len(), 8);
        assert_eq!(t.cds_start().unwrap(), 155158611);
        assert_eq!(t.cds_end().unwrap(), 155162634);

        let start = t.start_codon();
        assert_eq!(start.len(), 1);

        let stop = t.stop_codon();
        assert_eq!(stop.len(), 1);

        assert_eq!(t.exons()[0].is_coding(), true);
        assert_eq!(t.exons()[1].is_coding(), true);
        assert_eq!(t.exons()[2].is_coding(), true);
        assert_eq!(t.exons()[3].is_coding(), true);
        assert_eq!(t.exons()[4].is_coding(), true);
        assert_eq!(t.exons()[5].is_coding(), true);
        assert_eq!(t.exons()[6].is_coding(), true);
        assert_eq!(t.exons()[7].is_coding(), true);
    }
}
*/

#[cfg(test)]
mod test_ighm {
    use super::*;
    use crate::models::CdsStat;

    #[test]
    fn test_read() {
        let mut reader = Reader::from_file("tests/data/id-IGHM.gtf")
            .expect("Something failed reading the GTF file id-IGHM.gtf");
        let tr = match reader.transcripts() {
            Ok(res) => res,
            _ => panic!("No transcripts could be read"),
        };
        assert_eq!(tr.len(), 1);
        let t = tr.by_name("id-IGHM").unwrap();
        assert_eq!(t.exons().len(), 6);
        assert_eq!(t.cds_start().unwrap(), 106318298);
        assert_eq!(t.cds_end().unwrap(), 106322322);

        let start = t.start_codon();
        assert_eq!(start.len(), 1);

        let stop = t.stop_codon();
        assert_eq!(stop.len(), 1);

        assert_eq!(t.exons()[0].is_coding(), true);
        assert_eq!(t.exons()[1].is_coding(), true);
        assert_eq!(t.exons()[2].is_coding(), true);
        assert_eq!(t.exons()[3].is_coding(), true);
        assert_eq!(t.exons()[4].is_coding(), true);

        assert_eq!(t.cds_start_codon_stat(), CdsStat::Complete);
        assert_eq!(t.cds_stop_codon_stat(), CdsStat::Incomplete);
    }
}