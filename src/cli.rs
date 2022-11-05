use atglib::qc::QcCheck;
use atglib::qc::QcResult;
use clap::{Parser, ValueEnum};

/// Convert transcript data from and to different file formats
///
/// More detailed usage instructions on Github: <https://github.com/anergictcell/atg>
#[derive(Parser, Debug)]
#[command(author, version, about)]
pub struct Args {
    /// Format of input file
    #[arg(short, long, value_name = "FORMAT")]
    pub from: InputFormat,

    /// Output format
    #[arg(short, long, value_name = "FORMAT")]
    pub to: OutputFormat,

    /// Path to input file
    #[arg(short, long, default_value = "/dev/stdin", value_name = "FILE")]
    pub input: String,

    /// Path to output file
    #[arg(short, long, default_value = "/dev/stdout", value_name = "FILE")]
    pub output: String,

    /// The feature source to indicate in GTF files (optional with `--output gtf`)
    #[arg(short, long, default_value = env!("CARGO_PKG_NAME"), value_name = "FILE")]
    pub gtf_source: String,

    /// Path to reference genome fasta file. (required with `--output [fasta | fasta-split | feature-sequence | qc]`)
    ///
    /// You can also specify an S3 Uri (s3://mybucket/myfile.fasta), but reading from S3 is currently quite slow
    #[arg(short, long, value_name = "FASTA_FILE", required_if_eq_any([("to", "fasta"),("to", "fasta-split"),("to", "feature-sequence"),("to", "qc")]))]
    pub reference: Option<String>,

    /// Which part of the transcript to transcribe
    ///
    /// This option is only needed when generating fasta output.
    #[arg(long, default_value = "cds")]
    pub fasta_format: FastaFormat,

    /// Sets the level of verbosity
    #[arg(short, action = clap::ArgAction::Count)]
    pub verbose: u8,

    /// Specify which genetic code to use for translating the transcripts
    ///
    /// Genetic codes can be specified globally (e.g. `-c standard`)
    ///
    /// or chromosome specific (e..g `-c "chrM:vertebrate mitochondrial"`).
    ///
    /// Specify by name or amino acid lookup table (e.g. `FFLLSSSSYY**CC*....`)
    ///
    /// Defaults to the standard genetic code for all transcripts. Suggested use for vertebrates:
    ///
    /// `-c standard -c "chrM:vertebrate mitochondrial"`)
    ///
    /// (optional with `--output qc`)
    #[arg(short = 'c', long, action = clap::ArgAction::Append, value_name = "GENETIC CODE")]
    pub genetic_code: Vec<String>,

    /// Remove all variants from the output that fail QC-checks
    ///
    /// You can specify one or multiple QC-checks. Only `NOK` results will be removed. `OK` and `NA` will remain.
    #[arg(short = 'q', long = "qc-check", action = clap::ArgAction::Append, value_name = "QC CHECKS", requires = "reference")]
    pub qc_check: Vec<QcFilter>,
}

#[derive(Clone, Debug, ValueEnum)]
pub enum FastaFormat {
    /// The full genomic sequence of the transcript, including introns. (similar to pre-processed mRNA)
    Transcript,
    /// All exons of the transcript. (similar to processed mRNA)
    Exons,
    /// The coding sequence of the transcript
    Cds,
}

impl FastaFormat {
    pub fn as_str(&self) -> &str {
        match self {
            FastaFormat::Transcript => "transcript",
            FastaFormat::Exons => "exons",
            FastaFormat::Cds => "cds",
        }
    }
}

#[derive(Clone, Debug, ValueEnum)]
pub enum InputFormat {
    /// GTF2.2 format
    Gtf,
    /// RefGene format (one transcript per line)
    Refgene,
    /// GenePredExt format (one transcript per line)
    Genepredext,
    /// ATG-specific binary format
    Bin,
}

impl std::fmt::Display for InputFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[derive(Clone, Debug, ValueEnum)]
pub enum OutputFormat {
    /// GTF2.2 format
    Gtf,
    /// RefGene format (one transcript per line)
    Refgene,
    /// GenePred format (one transcript per line)
    Genepred,
    /// GenePredExt format (one transcript per line)
    Genepredext,
    /// Bedfile (one transcript per line)
    Bed,
    /// Nucleotide sequence. There are multiple formatting options available, see --fasta-format
    Fasta,
    /// Like 'fasta', but every transcript is written to its own file. (--output must be the path to a folder)
    FastaSplit,
    /// Nucleotide sequence for every 'feature' (UTR, CDS or non-coding exons)
    FeatureSequence,
    /// Custom format, as needed for SpliceAI
    Spliceai,
    /// ATG-specific binary format
    Bin,
    /// Performs QC checks on all Transcripts
    Qc,
    /// No output
    None,
    /// This only makes sense for debugging purposes
    Raw,
}

impl std::fmt::Display for OutputFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[derive(Clone, Debug, Eq, PartialEq, ValueEnum)]
pub enum QcFilter {
    /// Transcript contains at least one exon
    Exon,
    /// The length of the CDS is divisible by 3
    CdsLength,
    /// The CDS starts with ATG
    Start,
    /// The CDS ends with a Stop codon
    Stop,
    /// TThe 5’UTR does not contain another start codon ATG
    UpstreamStart,
    /// The CDS does not contain another in-frame stop-codon
    UpstreamStop,
    /// The transcript is within the coordinates of the reference genome
    Coordinates,
}

impl QcFilter {
    pub fn remove(&self, qc: &QcCheck) -> bool {
        match self {
            QcFilter::Exon => qc.contains_exon() == QcResult::NOK,
            QcFilter::CdsLength => qc.correct_cds_length() == QcResult::NOK,
            QcFilter::Start => qc.correct_start_codon() == QcResult::NOK,
            QcFilter::Stop => qc.correct_stop_codon() == QcResult::NOK,
            QcFilter::UpstreamStart => qc.no_upstream_start_codon() == QcResult::NOK,
            QcFilter::UpstreamStop => qc.no_upstream_stop_codon() == QcResult::NOK,
            QcFilter::Coordinates => qc.correct_coordinates() == QcResult::NOK,
        }
    }
}

impl std::fmt::Display for QcFilter {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}
