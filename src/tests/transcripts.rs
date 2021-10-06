//! This module contains some standard transcripts to use for testing
//!
//! This simplifies writing unit- and integration tests.
//! Some transcripts in here have been choosen deliberately, because they
//! contain some edge cases that I tripped across during development.

use crate::models;
use crate::models::{Exon, Transcript, TranscriptBuilder};

/// NM_001365057.2
///
/// The last exon contains the last nucleotide of the stop codon.
/// Thus the exon is considered `non-coding` in GTF specs, but RefGene requires
/// a frame-offset value for the exon regardless.
pub fn nm_001365057() -> Transcript {
    let mut transcript = TranscriptBuilder::new()
        .name("NM_001365057.2")
        .chrom("chr9")
        .gene("C9orf85")
        .strand(models::Strand::Plus)
        .cds_start_codon_stat(models::CdsStat::Complete)
        .unwrap()
        .cds_stop_codon_stat(models::CdsStat::Complete)
        .unwrap()
        .build()
        .unwrap();

    transcript.push_exon(Exon::new(
        74526555,
        74526752,
        Some(74526651),
        Some(74526752),
        models::Frame::Zero,
    ));

    transcript.push_exon(Exon::new(
        74561922,
        74562028,
        Some(74561922),
        Some(74562028),
        models::Frame::Zero,
    ));

    transcript.push_exon(Exon::new(
        74597573,
        74600974,
        Some(74597573),
        Some(74597573),
        models::Frame::One,
    ));

    transcript
}

/// NM_001365408.1
///
/// The start codon is split across two exons.
pub fn nm_001365408() -> Transcript {
    let mut transcript = TranscriptBuilder::new()
        .name("NM_001365408.1")
        .chrom("chr16")
        .gene("CES2")
        .strand(models::Strand::Plus)
        .cds_start_codon_stat(models::CdsStat::Complete)
        .unwrap()
        .cds_stop_codon_stat(models::CdsStat::Complete)
        .unwrap()
        .build()
        .unwrap();

    transcript.push_exon(Exon::new(
        66969419,
        66969778,
        None,
        None,
        models::Frame::None,
    ));
    transcript.push_exon(Exon::new(
        66971940,
        66972144,
        Some(66972143),
        Some(66972144),
        models::Frame::Zero,
    ));
    transcript.push_exon(Exon::new(
        66973120,
        66973261,
        Some(66973120),
        Some(66973261),
        models::Frame::One,
    ));
    transcript.push_exon(Exon::new(
        66974125,
        66974258,
        Some(66974125),
        Some(66974258),
        models::Frame::Zero,
    ));
    transcript.push_exon(Exon::new(
        66974340,
        66974598,
        Some(66974340),
        Some(66974598),
        models::Frame::One,
    ));
    transcript.push_exon(Exon::new(
        66975027,
        66975125,
        Some(66975027),
        Some(66975125),
        models::Frame::Zero,
    ));
    transcript.push_exon(Exon::new(
        66975409,
        66975549,
        Some(66975409),
        Some(66975549),
        models::Frame::Zero,
    ));
    transcript.push_exon(Exon::new(
        66975671,
        66975751,
        Some(66975671),
        Some(66975751),
        models::Frame::Zero,
    ));
    transcript.push_exon(Exon::new(
        66976008,
        66976152,
        Some(66976008),
        Some(66976152),
        models::Frame::Zero,
    ));
    transcript.push_exon(Exon::new(
        66976551,
        66976640,
        Some(66976551),
        Some(66976640),
        models::Frame::Two,
    ));
    transcript.push_exon(Exon::new(
        66977202,
        66977274,
        Some(66977202),
        Some(66977274),
        models::Frame::Two,
    ));
    transcript.push_exon(Exon::new(
        66977742,
        66978999,
        Some(66977742),
        Some(66977928),
        models::Frame::One,
    ));

    transcript
}

/// NM_001371720.1
///
/// This transcript contains two book-ended exons.
/// This is causing issues in GTF parsing.
/// To help fixing the issue, this function requires an extra
/// boolean parameter to indicate if it should be parsed as GTF.
/// This should be removed once GTF parsing is fixed.
pub fn nm_001371720(gtf: bool) -> Transcript {
    let mut transcript = TranscriptBuilder::new()
        .name("NM_001371720.1")
        .chrom("chr1")
        .gene("MUC1")
        .strand(models::Strand::Minus)
        .cds_start_codon_stat(models::CdsStat::Complete)
        .unwrap()
        .cds_stop_codon_stat(models::CdsStat::Complete)
        .unwrap()
        .build()
        .unwrap();

    transcript.push_exon(Exon::new(
        155158300,
        155158685,
        Some(155158611),
        Some(155158685),
        models::Frame::Zero,
    ));
    transcript.push_exon(Exon::new(
        155159701,
        155159850,
        Some(155159701),
        Some(155159850),
        models::Frame::Zero,
    ));
    transcript.push_exon(Exon::new(
        155159931,
        155160052,
        Some(155159931),
        Some(155160052),
        models::Frame::Two,
    ));
    transcript.push_exon(Exon::new(
        155160198,
        155160334,
        Some(155160198),
        Some(155160334),
        models::Frame::One,
    ));
    transcript.push_exon(Exon::new(
        155160484,
        155160539,
        Some(155160484),
        Some(155160539),
        models::Frame::Zero,
    ));

    if gtf {
        transcript.push_exon(Exon::new(
            155160639,
            155162101,
            Some(155160639),
            Some(155162101),
            models::Frame::Two,
        ));
    } else {
        transcript.push_exon(Exon::new(
            155160639,
            155161619,
            Some(155160639),
            Some(155161619),
            models::Frame::Zero,
        ));
        transcript.push_exon(Exon::new(
            155161620,
            155162101,
            Some(155161620),
            Some(155162101),
            models::Frame::Two,
        ));
    }

    transcript.push_exon(Exon::new(
        155162577,
        155162700,
        Some(155162577),
        Some(155162634),
        models::Frame::Zero,
    ));
    transcript
}

/// NM_201550.4
///
/// Single exon transcript
pub fn nm_201550() -> Transcript {
    let mut transcript = TranscriptBuilder::new()
        .name("NM_201550.4")
        .chrom("chr12")
        .gene("LRRC10")
        .strand(models::Strand::Minus)
        .cds_start_codon_stat(models::CdsStat::Complete)
        .unwrap()
        .cds_stop_codon_stat(models::CdsStat::Complete)
        .unwrap()
        .build()
        .unwrap();

    transcript.push_exon(Exon::new(
        70002344,
        70004687,
        Some(70003785),
        Some(70004618),
        models::Frame::Zero,
    ));

    transcript
}

/// Generates a transcript to be used for tests.
/// It contains 5 exons, 3 of which are coding
///
/// ```text
/// Nucleotide positions:
///
///    1....   2....   3....   4....   5....
///    12345   12345   12345   12345   12345
/// ---=====---===XX---XXXXX---XXXX=---=====---
///
/// ---  Intron
/// ===  Exon (non-coding)
/// XXX  CDS
/// ```
pub fn standard_transcript() -> models::Transcript {
    let mut transcript = models::TranscriptBuilder::new()
        .name("Test-Transcript")
        .chrom("chr1")
        .strand(models::Strand::Plus)
        .gene("Test-Gene")
        .cds_start_stat(models::CdsStat::None)
        .cds_end_stat(models::CdsStat::None)
        .build()
        .unwrap();

    transcript.append_exons(&mut exons());
    transcript
}

fn exons() -> Vec<models::Exon> {
    vec![
        models::Exon::new(11, 15, None, None, models::Frame::None),
        models::Exon::new(21, 25, Some(24), Some(25), models::Frame::Zero),
        models::Exon::new(31, 35, Some(31), Some(35), models::Frame::One),
        models::Exon::new(41, 45, Some(41), Some(44), models::Frame::Two),
        models::Exon::new(51, 55, None, None, models::Frame::None),
    ]
}
