//! Contains some functions and Errors that are used all across `ATG`
pub mod errors;
mod genomic_relations;

pub use crate::utils::genomic_relations::{
    exon_cds_overlap, intersect, relation, subtract, union, GenomicRelation,
};
