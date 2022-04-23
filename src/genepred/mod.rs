//! Convert from/to GenePredExt
//!
//!
//! The GenePred format is described by [UCSC](http://genome.ucsc.edu/FAQ/FAQformat#format9)
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
//!
//! The format is almost identical to RefGene, it's only missing some column. So, instead of reinventing the wheel,
//! we copied most of refgene Writer code and just removed the extra columns.
//!
//! At the moment, there is only a GenePred `Writer`. `Reader` is not yet implemented.
//! Parsing GeneProd is not yet possible due to the missing exonFrames columns. This
//! could be calculated during parsing, but this is not yet done.

mod writer;

pub use crate::genepred::writer::Writer;
