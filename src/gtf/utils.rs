use crate::models::Transcript;
use crate::utils::errors::ParseGtfError;

use std::str::FromStr;

#[derive(Debug, PartialEq)]
pub struct Attributes {
    transcript: String,
    gene: String,
    others: Vec<(String, String)>,
}

impl Attributes {
    pub fn gene(&self) -> &str {
        &self.gene
    }

    pub fn transcript(&self) -> &str {
        &self.transcript
    }

    pub fn all(&self) -> Vec<(&str, &str)> {
        let mut res = Vec::with_capacity(self.others.len() + 2);
        res.push(("gene_id", self.gene()));
        res.push(("transcript_id", self.transcript()));
        for attr in &self.others {
            res.push((&attr.0, &attr.1));
        }
        res
    }

    pub fn from_transcript(transcript: &Transcript) -> Self {
        let name = transcript.name().to_string();
        let gene = transcript.gene().to_string();
        Self {
            transcript: name,
            gene,
            others: vec![],
        }
    }
}

impl FromStr for Attributes {
    type Err = ParseGtfError;

    fn from_str(s: &str) -> Result<Self, ParseGtfError> {
        let mut gene: Option<String> = None;
        let mut transcript: Option<String> = None;
        let mut others: Vec<(String, String)> = vec![];

        for attr in s.trim_end_matches(';').split(';') {
            match parse_attribute(attr) {
                Ok(("gene_id", value)) => gene = Some(value.to_string()),
                Ok(("transcript_id", value)) => transcript = Some(value.to_string()),
                Ok((key, value)) => others.push((key.to_string(), value.to_string())),
                Err(err) => {
                    return Err(ParseGtfError::from_chain(
                        err,
                        &format!("Unable to parse the attribute column\n\n>>>{}<<<\n\n", s),
                    ));
                }
            }
        }
        match (gene, transcript) {
            (Some(gene), Some(transcript)) => Ok(Attributes {
                gene,
                transcript,
                others,
            }),
            (None, None) => {
                return Err(ParseGtfError {
                    message: format!("missing gene_id and transcript_id in {}", s),
                })
            }
            (None, Some(_)) => {
                return Err(ParseGtfError {
                    message: format!("missing gene_id in {}", s),
                })
            }
            (Some(_), None) => {
                return Err(ParseGtfError {
                    message: format!("missing transcript_id in {}", s),
                })
            }
        }
    }
}

fn parse_attribute(attr: &str) -> Result<(&str, &str), ParseGtfError> {
    let items: Vec<&str> = attr.trim().splitn(2, ' ').collect();
    match items.len() {
        2 => Ok((
            items[0],
            items[1].trim_end_matches('\"').trim_start_matches('\"'),
        )),
        _ => Err(ParseGtfError {
            message: format!(
                "Unable to parse the attribute\n\n{}\nPlease check your GTF input.",
                attr
            ),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_attributes_parsing() {
        let col = "gene_id \"ZBTB16\"; transcript_id \"NM_001354751.2\"; exon_number \"2\"; exon_id \"NM_001354751.2.2\"; gene_name \"ZBTB16\";";
        let attr = Attributes::from_str(col).unwrap();
        assert_eq!(attr.transcript(), "NM_001354751.2");
        assert_eq!(attr.gene(), "ZBTB16");
        assert_eq!(attr.all().len(), 5);

        let col = "gene_id \"ZBTB16\"; transcript_id \"NM_001354752.1\"; gene_name \"ZBTB16\";";
        let attr = Attributes::from_str(col).unwrap();
        assert_eq!(attr.transcript(), "NM_001354752.1");
        assert_eq!(attr.gene(), "ZBTB16");
        assert_eq!(attr.all().len(), 3);
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
