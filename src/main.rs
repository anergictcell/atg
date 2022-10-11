#[macro_use]
extern crate log;
use std::fs::File;
use std::process;

use bincode::{deserialize_from, serialize_into};
use clap::{Parser, ValueEnum};

use s3reader::{S3ObjectUri, S3Reader};

use atglib::bed;
use atglib::fasta;
use atglib::fasta::FastaReader;
use atglib::genepred;
use atglib::genepredext;
use atglib::gtf;
use atglib::models::{GeneticCode, TranscriptWrite, Transcripts};
use atglib::qc;
use atglib::read_transcripts;
use atglib::refgene;
use atglib::spliceai;
use atglib::utils::errors::AtgError;

/// Convert transcript data from and to different file formats
///
/// More detailed usage instructions on Github: <https://github.com/anergictcell/atg>
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// Format of input file
    #[arg(short, long, value_name = "FORMAT")]
    from: InputFormat,

    /// Output format
    #[arg(short, long, value_name = "FORMAT")]
    to: OutputFormat,

    /// Path to input file
    #[arg(short, long, default_value = "/dev/stdin", value_name = "FILE")]
    input: String,

    /// Path to output file
    #[arg(short, long, default_value = "/dev/stdout", value_name = "FILE")]
    output: String,

    /// The feature source to indicate in GTF files (optional with `--output gtf`)
    #[arg(short, long, default_value = env!("CARGO_PKG_NAME"), value_name = "FILE")]
    gtf_source: String,

    /// Path to reference genome fasta file. (required with `--output [fasta | fasta-split | feature-sequence | qc]`)
    ///
    /// You can also specify an S3 Uri (s3://mybucket/myfile.fasta), but reading from S3 is currently quite slow
    #[arg(short, long, value_name = "FASTA_FILE", required_if_eq_any([("to", "fasta"),("to", "fasta-split"),("to", "feature-sequence"),("to", "qc")]))]
    reference: Option<String>,

    /// Which part of the transcript to transcribe
    ///
    /// This option is only needed when generating fasta output.
    #[arg(long, default_value = "cds")]
    fasta_format: FastaFormat,

    /// Sets the level of verbosity
    #[arg(short, action = clap::ArgAction::Count)]
    verbose: u8,

    /// Specify which genetic code to use for translating the transcripts
    ///
    /// Genetic codes can be specified globally (e.g. `-c standard`)
    ///
    /// or chromosome specific (e..g `-c chrM:vertebrate_mitochondria`).
    ///
    /// Specify by name or amino acid lookup table (e.g. `FFLLSSSSYY**CC*....`)
    ///
    /// Defaults to the standard genetic code for all transcripts. Suggested use for vertebrates:
    ///
    /// `-c standard -c chrM:vertebrate_mitochondria`)
    ///
    /// (optional with `--output qc`)
    #[arg(short = 'c', long, action = clap::ArgAction::Append, value_name = "GENETIC CODE")]
    genetic_code: Vec<String>,
}

#[derive(Clone, Debug, ValueEnum)]
enum FastaFormat {
    /// The full genomic sequence of the transcript, including introns. (similar to pre-processed mRNA)
    Transcript,
    /// All exons of the transcript. (similar to processed mRNA)
    Exons,
    /// The coding sequence of the transcript
    Cds,
}

impl FastaFormat {
    fn as_str(&self) -> &str {
        match self {
            FastaFormat::Transcript => "transcript",
            FastaFormat::Exons => "exons",
            FastaFormat::Cds => "cds",
        }
    }
}

#[derive(Clone, Debug, ValueEnum)]
enum InputFormat {
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
enum OutputFormat {
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


// There will be only a single instance of this enum
// so we can allow a large variant
#[allow(clippy::large_enum_variant)]
/// ReadSeekWrapper is an enum to allow dynamic assignment of either File or S3 Readers
/// to be used in the Reader objects of Atglib.
enum ReadSeekWrapper {
    File(File, String),
    S3(S3Reader, String),
}

impl ReadSeekWrapper {
    pub fn from_filename(filename: &str) -> Self {
        if filename.starts_with("s3://") {
            let uri = S3ObjectUri::new(filename).unwrap();
            let s3obj = S3Reader::open(uri).unwrap();
            Self::S3(s3obj, filename.to_string())
        } else {
            Self::File(File::open(filename).unwrap(), filename.to_string())
        }
    }

    pub fn from_cli_arg(filename: &Option<&str>) -> Result<ReadSeekWrapper, AtgError> {
        if let Some(filename) = filename {
            if filename.starts_with("s3://") {
                let uri = S3ObjectUri::new(filename).unwrap();
                let s3obj = S3Reader::open(uri).unwrap();
                Ok(Self::S3(s3obj, filename.to_string()))
            } else {
                Ok(Self::File(
                    File::open(filename).unwrap(),
                    filename.to_string(),
                ))
            }
        } else {
            Err(AtgError::new("No file specified"))
        }
    }

    pub fn filename(&self) -> &str {
        match self {
            ReadSeekWrapper::File(_, fname) => fname,
            ReadSeekWrapper::S3(_, fname) => fname,
        }
    }
}

// forward all custom implementations straight to the actual reader
impl std::io::Read for ReadSeekWrapper {
    fn read(&mut self, buf: &mut [u8]) -> Result<usize, std::io::Error> {
        match self {
            ReadSeekWrapper::S3(r, _) => r.read(buf),
            ReadSeekWrapper::File(r, _) => r.read(buf),
        }
    }

    fn read_to_end(&mut self, buf: &mut Vec<u8>) -> Result<usize, std::io::Error> {
        match self {
            ReadSeekWrapper::S3(r, _) => r.read_to_end(buf),
            ReadSeekWrapper::File(r, _) => r.read_to_end(buf),
        }
    }

    fn read_to_string(&mut self, buf: &mut String) -> Result<usize, std::io::Error> {
        match self {
            ReadSeekWrapper::S3(r, _) => r.read_to_string(buf),
            ReadSeekWrapper::File(r, _) => r.read_to_string(buf),
        }
    }
}

impl std::io::Seek for ReadSeekWrapper {
    fn seek(&mut self, pos: std::io::SeekFrom) -> Result<u64, std::io::Error> {
        match self {
            ReadSeekWrapper::S3(r, _) => r.seek(pos),
            ReadSeekWrapper::File(r, _) => r.seek(pos),
        }
    }
}


fn read_input_file(args: &Args) -> Result<Transcripts, AtgError> {
    let input_format = &args.from;
    let input_fd = &args.input;
    debug!("Reading {} transcripts from {}", input_format, input_fd);

    let transcripts = match input_format {
        InputFormat::Refgene => read_transcripts(refgene::Reader::from_file(input_fd))?,
        InputFormat::Genepredext => read_transcripts(genepredext::Reader::from_file(input_fd))?,
        InputFormat::Gtf => read_transcripts(gtf::Reader::from_file(input_fd))?,
        InputFormat::Bin => {
            let reader = File::open(input_fd)?;
            match deserialize_from(reader) {
                Ok(res) => res,
                Err(err) => return Err(AtgError::new(err)),
            }
        }
    };

    debug!(
        "Finished parsing input data. Found {} transcripts",
        transcripts.len()
    );
    Ok(transcripts)
}

fn write_output(args: &Args, transcripts: Transcripts) -> Result<(), AtgError> {
    let output_fd = &args.output;
    let output_format = &args.to;

    let fasta_format = &args.fasta_format;
    let fasta_reference = &args.reference;

    debug!("Writing transcripts as {} to {}", output_format, output_fd);

    let fastareader = get_fasta_reader(&fasta_reference.as_deref());

    match output_format {
        OutputFormat::Refgene => {
            let mut writer = refgene::Writer::from_file(output_fd)?;
            writer.write_transcripts(&transcripts)?
        }
        OutputFormat::Genepred => {
            let mut writer = genepred::Writer::from_file(output_fd)?;
            writer.write_transcripts(&transcripts)?
        }
        OutputFormat::Genepredext => {
            let mut writer = genepredext::Writer::from_file(output_fd)?;
            writer.write_transcripts(&transcripts)?
        }
        OutputFormat::Gtf => {
            let mut writer = gtf::Writer::from_file(output_fd)?;
            writer.set_source(&args.gtf_source);
            writer.write_transcripts(&transcripts)?
        }
        OutputFormat::Bed => {
            let mut writer = bed::Writer::from_file(output_fd)?;
            writer.write_transcripts(&transcripts)?
        }
        OutputFormat::Fasta => {
            let mut writer = fasta::Writer::from_file(output_fd)?;
            writer.fasta_reader(
                fastareader.ok_or_else(|| AtgError::new("unable to open fasta file"))?,
            );
            writer.fasta_format(fasta_format.as_str());
            writer.write_transcripts(&transcripts)?
        }
        OutputFormat::FastaSplit => {
            let outdir = std::path::Path::new(&output_fd);
            if !outdir.is_dir() {
                return Err(AtgError::new(
                    "fasta-split requires a directory as --output option",
                ));
            }
            let mut writer = fasta::Writer::from_file("/dev/null")?;
            writer.fasta_reader(
                fastareader.ok_or_else(|| AtgError::new("unable to open fasta file"))?,
            );
            writer.fasta_format(fasta_format.as_str());

            for tx in transcripts {
                let outfile = outdir.join(format!("{}.fasta", tx.name()));
                *writer.inner_mut() = std::io::BufWriter::new(File::create(outfile)?);
                writer.writeln_single_transcript(&tx)?;
            }
        }
        OutputFormat::FeatureSequence => {
            let mut writer = fasta::Writer::from_file(output_fd)?;
            writer.fasta_reader(
                fastareader.ok_or_else(|| AtgError::new("unable to open fasta file"))?,
            );
            for tx in transcripts {
                writer.write_features(&tx)?
            }
        }
        OutputFormat::Spliceai => {
            let mut writer = spliceai::Writer::from_file(output_fd)?;
            writer.write_transcripts(&transcripts)?
        }
        OutputFormat::Qc => {
            let mut writer = qc::Writer::from_file(output_fd)?;
            let genetic_code_arg = &args.genetic_code;
            add_genetic_code(genetic_code_arg, &mut writer)?;
            writer.fasta_reader(
                fastareader.ok_or_else(|| AtgError::new("unable to open fasta file"))?,
            );
            writer.write_header()?;
            writer.write_transcripts(&transcripts)?
        }
        OutputFormat::Bin => {
            let writer = File::create(output_fd)?;
            match serialize_into(&writer, &transcripts) {
                Ok(res) => res,
                Err(err) => return Err(AtgError::new(err)),
            }
        }
        OutputFormat::Raw => {
            for t in transcripts {
                println!("{}", t);
                for exon in t.exons() {
                    println!("{}", exon)
                }
            }
        }
        OutputFormat::None => {}
    };

    Ok(())
}

fn get_fasta_reader(filename: &Option<&str>) -> Option<FastaReader<ReadSeekWrapper>> {
    // Both fasta_reader and fai_reader are Result<ReadSeekWrapper> instances
    // so they can be accessed via `?`
    let fasta_reader = ReadSeekWrapper::from_cli_arg(filename);
    let fai_reader = match &fasta_reader {
        Ok(r) => Ok(ReadSeekWrapper::from_filename(&format!(
            "{}.fai",
            r.filename()
        ))),
        Err(err) => Err(AtgError::new(format!("invalid fasta reader: {}", err))),
    };

    if let (Ok(fasta_wrapper), Ok(fai_wrapper)) = (fasta_reader, fai_reader) {
        Some(FastaReader::from_reader(fasta_wrapper, fai_wrapper).unwrap())
    } else {
        None
    }
}

fn add_genetic_code<W: std::io::Write, R: std::io::Read + std::io::Seek>(
    genetic_code_arg: &Vec<String>,
    writer: &mut qc::Writer<W, R>,
) -> Result<(), AtgError> {
    for genetic_code_value in genetic_code_arg {
        match genetic_code_value.split_once(':') {
            // if the value contains a `:`, it is a key:value pair
            // for chromosome:genetic_code.
            Some((chrom, seq)) => {
                let gen_code = GeneticCode::guess(seq)?;
                debug!("Adding genetic code {} for {}", gen_code, chrom);
                writer.add_genetic_code(chrom.to_string(), gen_code);
            }
            // Without `:` the genetic code is used as default
            None => {
                let gen_code = GeneticCode::guess(genetic_code_value)?;
                debug!("Setting default genetic code to {}", gen_code);
                writer.default_genetic_code(gen_code);
            }
        }
    }
    Ok(())
}

fn main() {
    let cli_commands = Args::parse();

    loggerv::init_with_verbosity(cli_commands.verbose.into()).unwrap();

    let transcripts = match read_input_file(&cli_commands) {
        Ok(x) => x,
        Err(err) => {
            println!("\x1b[1;31mError:\x1b[0m {}", err);
            println!("\nPlease check `atg --help` for more options\n");
            process::exit(1);
        }
    };

    match write_output(&cli_commands, transcripts) {
        Ok(_) => debug!("All done here."),
        Err(err) => {
            println!("\x1b[1;31mError:\x1b[0m {}", err);
            println!("\nPlease check `atg --help` for more options\n");
            process::exit(1);
        }
    }
}
