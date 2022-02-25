//! Read and write Fasta files
//!
//! Please note that this module is special in that Reader and Writer are
//! not complimentary methods:
//!
//! The [`Writer`] writes the genomic, cDNA or coding sequence of transcripts
//! into fasta format. This writer is a standard `ATG` writer, similar to
//! [`GtfWriter`](`crate::gtf::Writer`) and [`RefGeneWriter`](`crate::refgene::Writer`).
//! It is slightly different in that it also requires access to a reference genome -
//! through the [`FastaReader`].
//!
//! The [`FastaReader`] is only meant to read Fasta sequences of reference genomes
//! and it cannot be used to build [`Transcripts`](`crate::models::Transcript`).

mod reader;
mod writer;
pub use reader::FastaReader;
pub use writer::Writer;
