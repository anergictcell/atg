#[macro_use]
extern crate log;
use std::fs::File;
use std::process;

use bincode::{deserialize_from, serialize_into};
use clap::{App, Arg, ArgMatches};

use atglib::fasta::FastaReader;
use atglib::models::Transcripts;
use atglib::utils::errors::AtgError;

use atglib::bed;
use atglib::fasta;
use atglib::genepred;
use atglib::genepredext;
use atglib::gtf;
use atglib::models::{GeneticCode, TranscriptWrite};
use atglib::qc;
use atglib::read_transcripts;
use atglib::refgene;
use atglib::spliceai;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

fn parse_cli_args() -> ArgMatches {
    App::new(env!("CARGO_PKG_NAME"))
        .version(VERSION)
        .author(env!("CARGO_PKG_AUTHORS"))
        .about(
            "Convert transcript data from and to different file formats\n\n\
            ATG supports the following formats:\n\
              * genepred         - GenePred format (one transcript per line\n\
              * genepredext      - GenePredExt format (one transcript per line\n\
              * refgene          - RefGene format (one transcript per line)\n\
              * gtf              - GTF2.2 format\n\
              * bed              - Bedfile (one transcript per line)\n\
              * fasta            - Nucleotide sequence. There are multiple formatting options available\n\
              * fasta-split      - Like 'fasta', but every transcript is written to its own file. \
              (--output must be the path to a folder)\n\
              * feature-sequence - Nucleotide sequence for every 'feature' (UTR, CDS or non-coding exons)\n\
              * spliceai         - Custom format, as needed for SpliceAI\n\
              * qc               - Performs QC checks on all Transcripts\n\
              * bin              - Binary format"
        )
        .after_help(format!(
            "For more detailed usage information visit {}\n",
            env!("CARGO_PKG_HOMEPAGE")
        ).as_str())
        .arg(
            Arg::with_name("from")
                .short('f')
                .long("from")
                .possible_values(&["genepredext", "refgene", "gtf", "bin"])
                .case_insensitive(true)
                .value_name("file-format")
                .help("Data format of input file")
                .display_order(1)
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("to")
                .short('t')
                .long("to")
                .possible_values(&["genepred", "genepredext", "refgene", "gtf", "bed", "fasta", "fasta-split", "feature-sequence", "spliceai", "raw", "bin", "qc", "none"])
                .case_insensitive(true)
                .value_name("file-format")
                .help("data format of the output")
                .display_order(3)
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("input")
                .short('i')
                .long("input")
                .value_name("/path/to/input/file")
                .help("Path to input file")
                .display_order(2)
                .takes_value(true)
                .default_value("/dev/stdin"),
        )
        .arg(
            Arg::with_name("output")
                .short('o')
                .long("output")
                .value_name("/path/to/output/file")
                .help("Path to output file")
                .display_order(4)
                .takes_value(true)
                .default_value("/dev/stdout"),
        )
        .arg(
            Arg::with_name("gtf_source")
                .long("gtf-source")
                .short('g')
                .value_name("source")
                .help("The feature source to indicate in GTF files (optional with `--output gtf`)")
                .display_order(1000)
                .takes_value(true)
                .required_if("to", "gtf")
                .default_value(env!("CARGO_PKG_NAME")),
        )
        .arg(
            Arg::with_name("fasta_reference")
                .long("reference")
                .short('r')
                .value_name("/path/to/reference.fasta")
                .help("Path to reference genome fasta file. (required with `--output [fasta | fasta-split | feature-sequence]`)")
                .display_order(1001)
                .takes_value(true)
                .required_ifs(&[
                    ("to", "fasta"),
                    ("to", "fasta-split"),
                    ("to", "feature-sequence"),
                    ("to", "qc")
                ])
        )
        .arg(
            Arg::with_name("fasta_format")
                .long("fasta-format")
                .possible_values(&["transcript", "exons", "cds"])
                .value_name("format")
                .help(
                    "Which part of the transcript to transcribe. \
                    (This option is only needed when generating \
                    fasta output)\n\
                    * transcript - The full genomic sequence of the transcript \
                    including introns. (similar to pre-processed mRNA)\n\
                    * cds        - The coding sequence of the transcript.\n\
                    * exons      - All exons of the transcript. (similar to processed mRNA)",
                )
                .display_order(1002)
                .takes_value(true)
                .default_value("cds")
                .required_ifs(&[
                    ("to", "fasta"),
                    ("to", "fasta-split")
                ]),
        )
        .arg(
            Arg::with_name("v")
                .short('v')
                .action(clap::ArgAction::Count)
                .help("Sets the level of verbosity"),
        )
        .arg(
            Arg::with_name("genetic_code")
            .long("genetic-code")
            .short('c')
            .value_name("vertebrate_mitochondria")
            .display_order(1003)
            .takes_value(true)
            .action(clap::ArgAction::Append)
            .help(
                "Specify which genetic code \
                to use for translating the transcripts. \
                Genetic codes can be specified globally (e.g. `-c standard`) \
                or chromosome specific (e..g `-c chrM:vertebrate_mitochondria`). \
                Genetic codes can be specified via their name or by specifying the \
                amino acid lookup table directly. \
                You can specify multiple genetic codes. \
                Default is the standard genetic code for all transcripts. Most \
                likely you want to use the default and specify the mitochondria code \
                separately: `--genetic-code standard --genetic-code chrM:vertebrate_mitochondria`)
                (optional with `--output qc`)"
            )
        )
        .get_matches()
}

fn read_input_file(args: &ArgMatches) -> Result<Transcripts, AtgError> {
    let input_format = args.value_of("from").unwrap();
    let input_fd = args.value_of("input").unwrap();
    debug!("Reading {} transcripts from {}", input_format, input_fd);

    let transcripts = match input_format {
        "refgene" => read_transcripts(refgene::Reader::from_file(input_fd))?,
        "genepredext" => read_transcripts(genepredext::Reader::from_file(input_fd))?,
        "gtf" => read_transcripts(gtf::Reader::from_file(input_fd))?,
        "bin" => {
            let reader = File::open(input_fd)?;
            match deserialize_from(reader) {
                Ok(res) => res,
                Err(err) => return Err(AtgError::new(err)),
            }
        }
        _ => {
            return Err(AtgError::from(format!(
                "Invalid file-format: {}",
                input_format
            )))
        }
    };

    debug!(
        "Finished parsing input data. Found {} transcripts",
        transcripts.len()
    );
    Ok(transcripts)
}

fn write_output(args: &ArgMatches, transcripts: Transcripts) -> Result<(), AtgError> {
    let output_fd = args.value_of("output").unwrap();
    let output_format = args.value_of("to").unwrap();

    let fasta_format = args.value_of("fasta_format");
    let fasta_reference = args.value_of("fasta_reference");

    debug!("Writing transcripts as {} to {}", output_format, output_fd);

    match output_format {
        "refgene" => {
            let mut writer = refgene::Writer::from_file(output_fd)?;
            writer.write_transcripts(&transcripts)?
        }
        "genepred" => {
            let mut writer = genepred::Writer::from_file(output_fd)?;
            writer.write_transcripts(&transcripts)?
        }
        "genepredext" => {
            let mut writer = genepredext::Writer::from_file(output_fd)?;
            writer.write_transcripts(&transcripts)?
        }
        "gtf" => {
            let mut writer = gtf::Writer::from_file(output_fd)?;
            writer.set_source(args.value_of("gtf_source").unwrap());
            writer.write_transcripts(&transcripts)?
        }
        "bed" => {
            let mut writer = bed::Writer::from_file(output_fd)?;
            writer.write_transcripts(&transcripts)?
        }
        "fasta" => {
            let mut writer = fasta::Writer::from_file(output_fd)?;
            writer.fasta_reader(FastaReader::from_file(fasta_reference.unwrap())?);
            writer.fasta_format(fasta_format.unwrap());
            writer.write_transcripts(&transcripts)?
        }
        "fasta-split" => {
            let outdir = std::path::Path::new(output_fd);
            if !outdir.is_dir() {
                return Err(AtgError::new(
                    "fasta-split requires a directory as --output option",
                ));
            }
            for tx in transcripts {
                let outfile = outdir.join(format!("{}.fasta", tx.name()));
                let mut writer = fasta::Writer::from_file(
                    outfile
                        .to_str()
                        .expect("ATG does not support non-UTF8 filename"),
                )?;
                writer.fasta_reader(FastaReader::from_file(fasta_reference.unwrap())?);
                writer.fasta_format(fasta_format.unwrap());
                writer.writeln_single_transcript(&tx)?;
            }
        }
        "feature-sequence" => {
            let mut writer = fasta::Writer::from_file(output_fd)?;
            writer.fasta_reader(FastaReader::from_file(fasta_reference.unwrap())?);
            for tx in transcripts {
                writer.write_features(&tx)?
            }
        }
        "spliceai" => {
            let mut writer = spliceai::Writer::from_file(output_fd)?;
            writer.write_transcripts(&transcripts)?
        }
        "qc" => {
            let mut writer = qc::Writer::from_file(output_fd)?;
            let genetic_code_arg = args.get_many("genetic_code");
            add_genetic_code(genetic_code_arg, &mut writer)?;
            writer.fasta_reader(FastaReader::from_file(fasta_reference.unwrap())?);
            writer.write_header()?;
            writer.write_transcripts(&transcripts)?
        }
        "bin" => {
            let writer = File::create(output_fd)?;
            match serialize_into(&writer, &transcripts) {
                Ok(res) => res,
                Err(err) => return Err(AtgError::new(err)),
            }
        }
        "raw" => {
            for t in transcripts {
                println!("{}", t);
                for exon in t.exons() {
                    println!("{}", exon)
                }
            }
        }
        "none" => {}
        _ => return Err(AtgError::new("Invalid >>to<< parameter")),
    };

    Ok(())
}

fn add_genetic_code<W: std::io::Write>(
    genetic_code_arg: Option<clap::parser::ValuesRef<'_, String>>,
    writer: &mut qc::Writer<W>,
) -> Result<(), AtgError> {
    for genetic_code_value in genetic_code_arg.unwrap_or_default() {
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
    let cli_commands = parse_cli_args();

    loggerv::init_with_verbosity(cli_commands.get_count("v").into()).unwrap();

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
