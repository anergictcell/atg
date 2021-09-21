// I could not find a proper documentation for the refGene format anywhere
// There are various trustworthy and untrustworthy defintions
// I mainly relied on
// https://www.biostars.org/p/18480/
// and the tip to "Describe table schema" on http://genome.ucsc.edu/cgi-bin/hgTables
//
//  Schema for NCBI RefSeq - RefSeq genes from NCBI
//     Database: hg38    Primary Table: ncbiRefSeq    Row Count: 173,266   Data last updated: 2021-05-28
// Format description: A gene prediction with some additional info.
// field   example SQL type    info    description
// bin     585 smallint(5) unsigned    range   Indexing field to speed chromosome range queries.
// name    NR_046018.2 varchar(255)    values  Name of gene (usually transcript_id from GTF)
// chrom   chr1    varchar(255)    values  Reference sequence chromosome or scaffold
// strand  +   char(1)     values  + or - for strand
// txStart     11873   int(10) unsigned    range   Transcription start position (or end position for minus strand item)
// txEnd   14409   int(10) unsigned    range   Transcription end position (or start position for minus strand item)
// cdsStart    14409   int(10) unsigned    range   Coding region start (or end position for minus strand item)
// cdsEnd  14409   int(10) unsigned    range   Coding region end (or start position for minus strand item)
// exonCount   3   int(10) unsigned    range   Number of exons
// exonStarts  11873,12612,13220,  longblob        Exon start positions (or end positions for minus strand item)
// exonEnds    12227,12721,14409,  longblob        Exon end positions (or start positions for minus strand item)
// score   0   int(11)     range   score
// name2   DDX11L1 varchar(255)    values  Alternate name (e.g. gene_id from GTF)
// cdsStartStat    none    enum('none', 'unk', 'incmpl', 'cmpl')   values  Status of CDS start annotation (none, unknown, incomplete, or complete)
// cdsEndStat  none    enum('none', 'unk', 'incmpl', 'cmpl')   values  Status of CDS end annotation (none, unknown, incomplete, or complete)
// exonFrames  -1,-1,-1,   longblob        Exon frame {0,1,2}, or -1 if no frame for exon

pub mod constants;
pub mod reader;
pub mod utils;
pub mod writer;

pub use crate::refgene::constants::*;
pub use crate::refgene::reader::Reader;
pub use crate::refgene::utils::ParseRefGeneError;
pub use crate::refgene::writer::Writer;

#[cfg(test)]
mod tests {
    use crate::models::Transcript;
    use crate::models::TranscriptWrite;
    use crate::refgene::Writer;
    use std::convert::TryFrom;
    use std::str;

    fn compose_line(transcript: &Transcript) -> String {
        let buf = Vec::new();
        let mut writer = Writer::new(buf);
        match writer.write_single_transcript(transcript) {
            Ok(_) => String::from(str::from_utf8(&writer.into_inner().unwrap()).unwrap()),
            Err(_x) => "Error".to_string(),
        }
    }

    #[test]
    fn test_transcript_building() {
        let cols = vec![
            "585",
            "NR_046018.2",
            "chr1",
            "+",
            "11873",
            "14409",
            "14409",
            "14409",
            "3",
            "11873,12612,13220,",
            "12227,12721,14409,",
            "0",
            "DDX11L1",
            "unk",
            "unk",
            "-1,-1,-1,",
        ];
        let transcript = Transcript::try_from(cols).unwrap();
        assert_eq!(transcript.name(), "NR_046018.2");
        assert_eq!(transcript.chrom(), "chr1");
        assert_eq!(transcript.exon_count(), 3);
    }

    #[test]
    fn test_parse_exons_with_cds() {
        let line = "1233\tNM_005274.2\tchr1\t-\t84964005\t84972262\t84967527\t84971774\t4\t84964005,84967508,84971693,84972118,\t84964231,84967653,84971984,84972262,\t0\tGNG5\tcmpl\tcmpl\t-1,0,0,-1,";

        let cols: Vec<&str> = line.trim().split('\t').collect();
        let transcript = Transcript::try_from(cols).unwrap();

        assert_eq!(transcript.exons().len(), 4);

        assert_eq!(transcript.exons()[0].start, 84964006);
        assert_eq!(transcript.exons()[0].end, 84964231);
        assert_eq!(transcript.exons()[0].cds_start, None);
        assert_eq!(transcript.exons()[0].cds_end, None);

        assert_eq!(transcript.exons()[1].start, 84967509);
        assert_eq!(transcript.exons()[1].end, 84967653);
        assert_eq!(transcript.exons()[1].cds_start, Some(84967528));
        assert_eq!(transcript.exons()[1].cds_end, Some(84967653));

        assert_eq!(transcript.exons()[2].start, 84971694);
        assert_eq!(transcript.exons()[2].end, 84971984);
        assert_eq!(transcript.exons()[2].cds_start, Some(84971694));
        assert_eq!(transcript.exons()[2].cds_end, Some(84971774));

        assert_eq!(transcript.exons()[3].start, 84972119);
        assert_eq!(transcript.exons()[3].end, 84972262);
        assert_eq!(transcript.exons()[3].cds_start, None);
        assert_eq!(transcript.exons()[3].cds_end, None);
    }

    #[test]
    fn test_parsing_and_composing() {
        let lines = vec![
            "1233\tNM_005274.2\tchr1\t-\t84964005\t84972262\t84967527\t84971774\t4\t84964005,84967508,84971693,84972118,\t84964231,84967653,84971984,84972262,\t0\tGNG5\tcmpl\tcmpl\t-1,0,0,-1,",
            "1091\tNM_002896.3\tchr11\t+\t66406087\t66413944\t66407182\t66411603\t4\t66406087,66407170,66410920,66413497,\t66406223,66407594,66411611,66413944,\t0\tRBM4\tcmpl\tcmpl\t-1,0,1,-1,",
            "185\tNM_001253849.1\tchr1\t-\t117686208\t117746458\t117690279\t117699355\t6\t117686208,117690234,117695712,117699195,117712728,117746239,\t117687847,117690404,117695991,117699543,117712793,117746458,\t0\tVTCN1\tcmpl\tcmpl\t-1,1,1,0,-1,-1,",
            "1368\tNM_001163814.1\tchr10\t+\t102747292\t102754158\t102749519\t102750782\t5\t102747292,102749400,102750192,102750625,102752946,\t102747370,102749641,102750300,102750811,102754158,\t0\tC10orf2\tcmpl\tcmpl\t-1,0,2,2,-1,",
            "1374\tNM_001206389.1\tchr10\t-\t103529886\t103540126\t103530085\t103534513\t5\t103529886,103531219,103534488,103535496,103540052,\t103530376,103531326,103534636,103535533,103540126,\t0\tFGF8\tcmpl\tcmpl\t0,1,0,-1,-1,",
            "1374\tNM_006119.4\tchr10\t-\t103529886\t103535759\t103530085\t103535657\t5\t103529886,103531219,103534488,103535496,103535625,\t103530376,103531326,103534669,103535533,103535759,\t0\tFGF8\tcmpl\tcmpl\t0,1,0,2,0,",
            "1392\tNM_145247.4\tchr10\t+\t105881815\t105886143\t105882748\t105885462\t4\t105881815,105882722,105883471,105885270,\t105881912,105882844,105883882,105886143,\t0\tSFR1\tcmpl\tcmpl\t-1,0,0,0,",
            "1431\tNR_046928.1\tchr10\t+\t110930414\t110999490\t110999490\t110999490\t3\t110930414,110930439,110999457,\t110930427,110930442,110999490,\t0\tRNU6-53\tunk\tunk\t-1,-1,-1,",
            "1437\tNR_038943.1\tchr10\t-\t111705316\t111768139\t111768139\t111768139\t5\t111705316,111706953,111713513,111765714,111767938,\t111706053,111707034,111713649,111765839,111768139,\t0\tLOC100505933\tunk\tunk\t-1,-1,-1,-1,-1,",
            "1444\tNR_026932.1\tchr10\t-\t112628647\t112630662\t112630662\t112630662\t1\t112628647,\t112630662,\t0\tLOC282997\tunk\tunk\t-1,"
        ];
        for line in lines {
            let cols: Vec<&str> = line.trim().split('\t').collect();
            let transcript = Transcript::try_from(cols).unwrap();

            assert_eq!(compose_line(&transcript), line);
        }
    }
}
