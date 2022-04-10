//! Convert from/to GenePredExt
//!
//!
//! The GenePred(ext) format is described by [UCSC](http://genome.ucsc.edu/FAQ/FAQformat#format9)
//!
//!  # Schema for NCBI RefSeq - RefSeq genes from NCBI
//! | Column | Type | Example | Description |
//! | --- | --- | --- | --- |
//! | name | str |  NR_046018.2 | Name of gene (usually transcript_id from GTF) |
//! | chrom | str | chr1 | Reference sequence chromosome or scaffold |
//! | strand | enum("+", "-") | + | + or - for strand |
//! | txStart | int | 11873 | Transcription start position (or end position for minus strand item) |
//! | txEnd | int | 14409 | Transcription end position (or start position for minus strand item) |
//! | cdsStart | int | 14409 | Coding region start (or end position for minus strand item) |
//! | cdsEnd | int | 4409 | Coding region end (or start position for minus strand item) |
//! | exonCount | int | 3 | Number of exons |
//! | exonStarts | List of int | 1873,12612,13220, | Exon start positions (or end positions for minus strand item) (with trailing comma) |
//! | exonEnds | List of int | 12227,12721,14409, | Exon end positions (or start positions for minus strand item) (with trailing comma) |
//! | score | int | 0 | The score field indicates a degree of confidence in the feature's existence and coordinates |
//! | name2 | str | DDX11L1 | Alternate name (e.g. gene_id from GTF) |
//! | cdsStartStat | enum("none", "unk", "incmpl", "cmpl") | compl | Status of CDS start annotation (none, unknown, incomplete, or complete) |
//! | cdsEndStat | enum("none", "unk", "incmpl", "cmpl") | compl | Status of CDS end annotation (none, unknown, incomplete, or complete) |
//! | exonFrames | List of enum(-1, 0, 1, 2) | -1,0,2,1 | Exon frame {0,1,2}, or -1 if no frame for exon |
//!
//! The format is almost identical to RefGene, it's only missing the `bin` column. So, instead of reinventing the wheel,
//! we copied most of refgene Reader and Writer code and just added/removed the `bin` column.

mod reader;
mod writer;

pub use crate::genepredext::reader::Reader;
pub use crate::genepredext::writer::Writer;
