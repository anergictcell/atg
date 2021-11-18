//! Convert from/to GTF
//!
//! The GTF section is written for GTF2.2 as defined by
//! [the Brent lab](https://mblab.wustl.edu/GTF22.html) and
//! [UCSC](http://genome.ucsc.edu/FAQ/FAQformat.html#format4).
//!
//! # GTF specs
//! GTF2.2 files must contain the following columns:
//!
//! | Column | Mandatory | Type | Explanation |
//! | --- | --- | --- | --- |
//! | seqname | Yes | str | The name of the sequence. Commonly, this is the chromosome ID or contig ID. Note that the coordinates used must be unique within each sequence name in all GTFs for an annotation set. |
//! | source  | Yes | str | The source column should be a unique label indicating where the annotations came from --- typically the name of either a prediction program or a public database. |
//! | feature  | Yes | enum("CDS", "start_codon", "stop_codon", "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS", "exon") | The following feature types are required: "CDS", "start_codon", "stop_codon". The features "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS" and "exon" are optional. All other features will be ignored. The types must have the correct capitalization shown here. |
//! | start | Yes | u32 | The leftmost position of the feature. Inclusive |
//! | end | Yes | u32 | The rightmost position of the feature. Inclusive |
//! | score | Yes | float | The score field indicates a degree of confidence in the feature's existence and coordinates |
//! | strand | Yes | enum("+", "-") | The strand of the feature |
//! | frame | Yes | enum(0, 1, 2) | The frame-offset of the feature |
//! | attributes | Yes | str | key=value list (separated by ; ) **Must contain gene_id and transcript_id** (all other fields will be ignored)|
//! | comments  | Optional | str | Additional comments about the feature (this column is ignored)|
//!

mod constants;
mod reader;
mod record;
mod transcript;
mod writer;

pub use crate::gtf::reader::Reader;
pub use crate::gtf::record::{GtfFeature, GtfRecord, GtfRecordBuilder};
use crate::gtf::transcript::GtfRecordsGroup;
pub use crate::gtf::writer::Writer;
