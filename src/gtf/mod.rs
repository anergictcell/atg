// The main source and spec for GTF file is here:
// https://mblab.wustl.edu/GTF22.html
// and here: http://genome.ucsc.edu/FAQ/FAQformat.html#format4
// It defines the following columns
// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
// seqname == chrom     varchar     The name of the sequence. Commonly, this is the chromosome ID or contig ID. Note that the coordinates used must be unique within each sequence name in all GTFs for an annotation set.
// source               varchar     The source column should be a unique label indicating where the annotations came from --- typically the name of either a prediction program or a public database.
// feature              enum("CDS", "start_codon", "stop_codon", "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS" and "exon") The following feature types are required: "CDS", "start_codon", "stop_codon". The features "5UTR", "3UTR", "inter", "inter_CNS", "intron_CNS" and "exon" are optional. All other features will be ignored. The types must have the correct capitalization shown here.
// <start>              u32
// <end>                u32
// <score>              u32
// <strand>             char        + or -
// <frame>              enum(0, 1, 2)
// [attributes]         varchar     key=value list (separated by ; ); Must contain gene_id and transcript_id
// comments             varchar     (optional)

pub mod constants;
pub mod reader;
pub mod record;
pub mod transcript;
pub mod utils;
pub mod writer;

pub use crate::gtf::reader::Reader;
pub use crate::gtf::record::{GtfFeature, GtfRecord, GtfRecordBuilder};
use crate::gtf::transcript::GtfRecordsGroup;
pub use crate::gtf::utils::{Attributes, ParseGtfError};
pub use crate::gtf::writer::{Composer, Writer};
pub use crate::models;
pub use crate::models::{Exon, Frame, Strand, Transcript};
