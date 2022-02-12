//! Convert to bed format
//!
//! bed does not provide a Reader instance, can only be used as output
mod record;
mod writer;

pub use crate::bed::record::BedLine;
pub use crate::bed::writer::Writer;
