//! Holds the main data models and structs that are used in ATG
//!

mod codon;
mod exon;
mod frame;
mod sequence;
mod transcript;
mod transcripts;
mod utils;

pub use crate::models::codon::Codon;
pub use crate::models::exon::Exon;
pub use crate::models::frame::Frame;
pub use crate::models::sequence::{Nucleotide, Sequence};
pub use crate::models::transcript::{Coordinate, CoordinateVector, Transcript, TranscriptBuilder};
pub use crate::models::transcripts::Transcripts;
pub use crate::models::utils::{CdsStat, Strand, TranscriptRead, TranscriptWrite};

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_transcript() {
        let a = TranscriptBuilder::new()
            .name("Test-Transcript")
            .chrom("chr8")
            .strand(Strand::Plus)
            .gene("Test-Gene")
            .cds_start_stat(CdsStat::None)
            .cds_end_stat(CdsStat::None)
            .build()
            .unwrap();
        assert_eq!(a.name(), "Test-Transcript");
        assert_eq!(a.chrom(), "chr8");
        // assert_eq!(a.strand, Strand::Plus);
        assert_eq!(a.gene(), "Test-Gene");
    }

    #[test]
    fn test_special_transcripts() {
        let _incompl = "33 NM_001042758.1  chr1    +   206516199   206629505   206516199   206628374   18  206516199,206557368,206566045,206566901,206574779,206578607,206579735,206589248,206592724,206594600,206603511,206610314,206611313,206613325,206619420,206623730,206626545,206628223,    206516362,206557431,206566261,206567030,206575004,206578709,206579933,206589333,206592752,206594625,206603572,206610392,206611448,206613431,206619645,206623810,206626726,206629505,    0   SRGAP2  incmpl  cmpl    2,0,0,0,0,0,0,0,1,2,0,1,1,1,2,2,1,2,";
    }

    #[test]
    fn test_frame_addition() {
        assert_eq!((Frame::None + Frame::Zero).unwrap(), Frame::Zero);
        assert_eq!((Frame::Zero + Frame::None).unwrap(), Frame::Zero);
        assert_eq!((Frame::None + Frame::Zero).unwrap(), Frame::Zero);
        assert_eq!((Frame::Zero + Frame::None).unwrap(), Frame::Zero);
        assert_eq!((Frame::Zero + Frame::Zero).unwrap(), Frame::Zero);

        assert_eq!((Frame::None + Frame::One).unwrap(), Frame::One);
        assert_eq!((Frame::One + Frame::None).unwrap(), Frame::One);
        assert_eq!((Frame::None + Frame::One).unwrap(), Frame::One);
        assert_eq!((Frame::One + Frame::None).unwrap(), Frame::One);
        assert_eq!((Frame::Zero + Frame::One).unwrap(), Frame::One);
        assert_eq!((Frame::One + Frame::Zero).unwrap(), Frame::One);

        assert_eq!((Frame::One + Frame::One).unwrap(), Frame::Two);
        assert_eq!((Frame::One + Frame::Two).unwrap(), Frame::Zero);
        assert_eq!((Frame::Two + Frame::One).unwrap(), Frame::Zero);

        assert_eq!((Frame::Two + Frame::Two).unwrap(), Frame::One);
    }

    #[test]
    fn test_exon_downstream_frame() {
        assert_eq!(
            Exon::new(1, 1, Some(1), Some(3), Frame::Zero)
                .downstream_frame()
                .unwrap(),
            Frame::Zero
        );

        assert_eq!(
            Exon::new(1, 1, Some(1), Some(4), Frame::Zero)
                .downstream_frame()
                .unwrap(),
            Frame::Two
        );

        assert_eq!(
            Exon::new(1, 1, Some(1), Some(5), Frame::Zero)
                .downstream_frame()
                .unwrap(),
            Frame::One
        );

        assert_eq!(
            Exon::new(1, 1, Some(1), Some(6), Frame::Zero)
                .downstream_frame()
                .unwrap(),
            Frame::Zero
        );

        assert_eq!(
            Exon::new(1, 1, Some(1), Some(3), Frame::One)
                .downstream_frame()
                .unwrap(),
            Frame::One
        );

        // XATG
        assert_eq!(
            Exon::new(1, 1, Some(1), Some(4), Frame::One)
                .downstream_frame()
                .unwrap(),
            Frame::Zero
        );

        // XATGX
        assert_eq!(
            Exon::new(1, 1, Some(1), Some(5), Frame::One)
                .downstream_frame()
                .unwrap(),
            Frame::Two
        );

        // XATGXX
        assert_eq!(
            Exon::new(1, 1, Some(1), Some(6), Frame::One)
                .downstream_frame()
                .unwrap(),
            Frame::One
        );

        // XXA
        assert_eq!(
            Exon::new(1, 1, Some(1), Some(3), Frame::Two)
                .downstream_frame()
                .unwrap(),
            Frame::Two
        );

        // XXAT
        assert_eq!(
            Exon::new(1, 1, Some(1), Some(4), Frame::Two)
                .downstream_frame()
                .unwrap(),
            Frame::One
        );

        // XXATG
        assert_eq!(
            Exon::new(1, 1, Some(1), Some(5), Frame::Two)
                .downstream_frame()
                .unwrap(),
            Frame::Zero
        );

        // XXATGX
        assert_eq!(
            Exon::new(1, 1, Some(1), Some(6), Frame::Two)
                .downstream_frame()
                .unwrap(),
            Frame::Two
        );
    }
}

#[cfg(test)]
// TODO: Add test for negative strand
mod test_start_codon {
    //       1....   2....   3....   4....   5....
    //       12345   12345   12345   12345   12345
    //    ---=====---===XX---XXXXX---XXXX=---=====---
    // 1. ---=====---ATGXX---XXXXX---XXXX=---=====--- >> all in one exon
    // 2. ---=====---=ATGX---XXXXX---XXXX=---=====--- >> all in one exon
    // 3. ---=====---==ATG---XXXXX---XXXX=---=====--- >> all in one exon
    // 4. ---=====---===AT---GXXXX---XXXX=---=====--- >> split
    // 5. ---=====---====A---TGXXX---XXXX=---=====--- >> split
    // 6. ---=====---====A-----T-----GXXX=---=====--- >> not really possible, but for the sake of it, let's consider it as well

    use super::*;
    use crate::tests::transcripts::standard_transcript;

    #[test]
    fn case_1() {
        let mut transcript = standard_transcript();
        *transcript.exons_mut()[1].cds_start_mut() = Some(21);
        let a = transcript.start_codon();

        assert_eq!(a.len(), 1);
        assert_eq!(a[0].0, 21);
        assert_eq!(a[0].1, 23);
    }

    #[test]
    fn case_2() {
        let mut transcript = standard_transcript();
        *transcript.exons_mut()[1].cds_start_mut() = Some(22);
        let a = transcript.start_codon();

        assert_eq!(a.len(), 1);
        assert_eq!(a[0].0, 22);
        assert_eq!(a[0].1, 24);
        assert_eq!(a[0].2, Frame::Zero);
    }

    #[test]
    fn case_3() {
        let mut transcript = standard_transcript();
        *transcript.exons_mut()[1].cds_start_mut() = Some(23);
        let a = transcript.start_codon();

        assert_eq!(a.len(), 1);
        assert_eq!(a[0].0, 23);
        assert_eq!(a[0].1, 25);
    }

    #[test]
    fn case_4() {
        let mut transcript = standard_transcript();
        *transcript.exons_mut()[1].cds_start_mut() = Some(24);
        let a = transcript.start_codon();

        assert_eq!(a.len(), 2);
        assert_eq!(a[0].0, 24);
        assert_eq!(a[0].1, 25);

        assert_eq!(a[1].0, 31);
        assert_eq!(a[1].1, 31);
        assert_eq!(a[1].2, Frame::One);
    }

    #[test]
    fn case_5() {
        let mut transcript = standard_transcript();
        *transcript.exons_mut()[1].cds_start_mut() = Some(25);
        let a = transcript.start_codon();

        assert_eq!(a.len(), 2);
        assert_eq!(a[0].0, 25);
        assert_eq!(a[0].1, 25);

        assert_eq!(a[1].0, 31);
        assert_eq!(a[1].1, 32);
        assert_eq!(a[1].2, Frame::Two);
    }

    #[test]
    fn case_6() {
        let mut transcript = standard_transcript();
        *transcript.exons_mut()[1].cds_start_mut() = Some(25);
        *transcript.exons_mut()[2].cds_start_mut() = Some(33);
        *transcript.exons_mut()[2].cds_end_mut() = Some(33);
        let a = transcript.start_codon();

        assert_eq!(a.len(), 3);
        assert_eq!(a[0].0, 25);
        assert_eq!(a[0].1, 25);

        assert_eq!(a[1].0, 33);
        assert_eq!(a[1].1, 33);

        assert_eq!(a[2].0, 41);
        assert_eq!(a[2].1, 41);
    }
}

#[cfg(test)]
// TODO: Add test for negative strand
mod test_stop_codon {
    //       1....   2....   3....   4....   5....
    //       12345   12345   12345   12345   12345
    //    ---=====---===XX---XXXXX---XXXX=---=====---
    // 1. ---=====---===XX---XXXXX---XXUAG---=====--- >> all in one exon
    // 2. ---=====---=XXXX---XXXXX---XUAG=---=====--- >> all in one exon
    // 3. ---=====---==XXX---XXXXX---UAGX=---=====--- >> all in one exon
    // 4. ---=====---===XX---XXXXU---AGXX=---=====--- >> split
    // 5. ---=====---====X---XXXUA---GXXX=---=====--- >> split
    // 6. ---=====---===XU-----A-----G====---=====--- >> not really possible, but for the sake of it, let's consider it as well

    use super::*;
    use crate::tests::transcripts::standard_transcript;

    #[test]
    fn case_1() {
        let mut transcript = standard_transcript();
        *transcript.exons_mut()[3].cds_end_mut() = Some(45);
        let a = transcript.stop_codon();

        assert_eq!(a.len(), 1);
        assert_eq!(a[0].0, 43);
        assert_eq!(a[0].1, 45);
        assert_eq!(a[0].2, Frame::Zero);
    }

    #[test]
    fn case_2() {
        let mut transcript = standard_transcript();
        *transcript.exons_mut()[3].cds_end_mut() = Some(44);
        let a = transcript.stop_codon();

        assert_eq!(a.len(), 1);
        assert_eq!(a[0].0, 42);
        assert_eq!(a[0].1, 44);
    }

    #[test]
    fn case_3() {
        let mut transcript = standard_transcript();
        *transcript.exons_mut()[3].cds_end_mut() = Some(43);
        let a = transcript.stop_codon();

        assert_eq!(a.len(), 1);
        assert_eq!(a[0].0, 41);
        assert_eq!(a[0].1, 43);
    }

    #[test]
    fn case_4() {
        let mut transcript = standard_transcript();
        *transcript.exons_mut()[3].cds_end_mut() = Some(42);
        let a = transcript.stop_codon();

        assert_eq!(a.len(), 2);
        assert_eq!(a[0].0, 35);
        assert_eq!(a[0].1, 35);
        assert_eq!(a[0].2, Frame::Zero);

        assert_eq!(a[1].0, 41);
        assert_eq!(a[1].1, 42);
        assert_eq!(a[1].2, Frame::Two);
    }

    #[test]
    fn case_5() {
        let mut transcript = standard_transcript();
        *transcript.exons_mut()[3].cds_end_mut() = Some(41);
        let a = transcript.stop_codon();

        assert_eq!(a.len(), 2);
        assert_eq!(a[0].0, 34);
        assert_eq!(a[0].1, 35);
        assert_eq!(a[0].2, Frame::Zero);

        assert_eq!(a[1].0, 41);
        assert_eq!(a[1].1, 41);
        assert_eq!(a[1].2, Frame::One);
    }

    #[test]
    fn case_6() {
        let mut transcript = standard_transcript();
        *transcript.exons_mut()[2].cds_start_mut() = Some(33);
        *transcript.exons_mut()[2].cds_end_mut() = Some(33);
        *transcript.exons_mut()[3].cds_end_mut() = Some(41);
        let a = transcript.stop_codon();

        assert_eq!(a.len(), 3);
        assert_eq!(a[0].0, 25);
        assert_eq!(a[0].1, 25);

        assert_eq!(a[1].0, 33);
        assert_eq!(a[1].1, 33);

        assert_eq!(a[2].0, 41);
        assert_eq!(a[2].1, 41);
    }
}

#[cfg(test)]
mod test_codon {
    use super::*;
    use crate::models::codon::Codon;
    use crate::tests::transcripts::standard_transcript;

    #[test]
    fn test_codon_builder() {
        let transcript = standard_transcript();

        let codon = Codon::from_transcript(&transcript, &25).unwrap();
        assert_eq!(codon.fragments().len(), 2);
        assert_eq!(codon.fragments()[0].frame_offset(), Frame::Zero);
        assert_eq!(codon.fragments()[1].frame_offset(), Frame::Two);

        let codon = Codon::from_transcript(&transcript, &34).unwrap();
        assert_eq!(codon.fragments().len(), 2);
        assert_eq!(codon.fragments()[0].frame_offset(), Frame::Zero);
        assert_eq!(codon.fragments()[1].frame_offset(), Frame::One);
    }
}
