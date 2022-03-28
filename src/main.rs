#[macro_use]
extern crate log;
use std::fs::File;
use std::process;

use bincode::{deserialize_from, serialize_into};
use clap::{App, Arg, ArgMatches};

use atg::fasta::FastaReader;
use atg::models::Transcripts;
use atg::utils::errors::AtgError;

use atg::bed;
use atg::fasta;
use atg::gtf;
use atg::models::TranscriptWrite;
use atg::read_transcripts;
use atg::refgene;

fn parse_cli_args() -> ArgMatches<'static> {
    App::new(env!("CARGO_PKG_NAME"))
        .version(atg::VERSION)
        .author(env!("CARGO_PKG_AUTHORS"))
        .about(
            "Convert transcript data from and to different file formats\n\n\
            ATG supports the following formats:\n\
              * refgene          - RefGene format (one transcript per line)\n\
              * gtf              - GTF2.2 format\n\
              * bed              - Bedfile (one transcript per line)\n\
              * fasta            - Nucleotide sequence. There are multiple formatting options available\n\
              * feature-sequence - Nucleotide sequence for every feature\n\
              * bin              - Binary format"
        )
        .after_help(format!(
            "For more detailed usage information visit {}\n",
            env!("CARGO_PKG_HOMEPAGE")
        ).as_str())
        .arg(
            Arg::with_name("from")
                .short("f")
                .long("from")
                .possible_values(&["refgene", "gtf", "bin"])
                .case_insensitive(true)
                .value_name("file-format")
                .help("Data format of input file")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("to")
                .short("t")
                .long("to")
                .possible_values(&["refgene", "gtf", "bed", "fasta", "fasta-split", "feature-sequence", "raw", "bin", "none"])
                .case_insensitive(true)
                .value_name("file-format")
                .help("data format of the output")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("input")
                .short("i")
                .long("input")
                .value_name("/path/to/input/file")
                .help("Path to input file")
                .takes_value(true)
                .default_value("/dev/stdin"),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .value_name("/path/to/output/file")
                .help("Path to output file")
                .takes_value(true)
                .default_value("/dev/stdout"),
        )
        .arg(
            Arg::with_name("gtf_source")
                .long("gtf-source")
                .short("-g")
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
                .short("-r")
                .value_name("/path/to/reference.fasta")
                .help("Path to reference genome fasta file. (required with `--output [fasta | fasta-split | feature-sequence]`)")
                .display_order(1000)
                .takes_value(true)
                .required_ifs(&[
                    ("to", "fasta"),
                    ("to", "fasta-split"),
                    ("to", "feature-sequence")
                ])
        )
        .arg(
            Arg::with_name("fasta_format")
                .long("fasta-format")
                .possible_values(&["transcript", "exons", "cds"])
                .value_name("format")
                .help(
                    "Which part of the transcript to translate. \
                       (This option is only needed when generating \
                       fasta output)",
                )
                .display_order(1001)
                .takes_value(true)
                .default_value("cds")
                .required_ifs(&[
                    ("to", "fasta"),
                    ("to", "fasta-split")
                ]),
        )
        .arg(
            Arg::with_name("v")
                .short("v")
                .multiple(true)
                .help("Sets the level of verbosity"),
        )
        .get_matches()
}

fn read_input_file(input_format: &str, input_fd: &str) -> Result<Transcripts, AtgError> {
    debug!("Reading {} transcripts from {}", input_format, input_fd);

    let transcripts = match input_format {
        "refgene" => read_transcripts(refgene::Reader::from_file(input_fd))?,
        "gtf" => read_transcripts(gtf::Reader::from_file(input_fd))?,
        "bin" => {
            let reader = File::open(input_fd)?;
            deserialize_from(reader)?
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

fn write_output(
    output_format: &str,
    output_fd: &str,
    gtf_source: Option<&str>,
    fasta_reference: Option<&str>,
    fasta_format: Option<&str>,
    transcripts: Transcripts,
) -> Result<(), AtgError> {
    debug!("Writing transcripts as {} to {}", output_format, output_fd);

    let _ = match output_format {
        "refgene" => {
            let mut writer = refgene::Writer::from_file(output_fd)?;
            writer.write_transcripts(&transcripts)?
        }
        "gtf" => {
            let mut writer = gtf::Writer::from_file(output_fd)?;
            writer.set_source(gtf_source.unwrap());
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
        "bin" => {
            let writer = File::create(output_fd)?;
            serialize_into(&writer, &transcripts)?;
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

fn main() {
    let cli_commands = parse_cli_args();

    loggerv::init_with_verbosity(cli_commands.occurrences_of("v")).unwrap();

    let input_fd = cli_commands.value_of("input").unwrap();
    let output_fd = cli_commands.value_of("output").unwrap();

    let input_format = cli_commands.value_of("from").unwrap();
    let output_format = cli_commands.value_of("to").unwrap();

    let gtf_source = cli_commands.value_of("gtf_source");
    let fasta_format = cli_commands.value_of("fasta_format");
    let fasta_reference = cli_commands.value_of("fasta_reference");

    let transcripts = match read_input_file(input_format, input_fd) {
        Ok(x) => x,
        Err(err) => {
            println!("\x1b[1;31mError:\x1b[0m {}", err);
            println!("\nPlease check `atg --help` for more options\n");
            process::exit(1);
        }
    };

    match write_output(
        output_format,
        output_fd,
        gtf_source,
        fasta_reference,
        fasta_format,
        transcripts,
    ) {
        Ok(_) => debug!("All done here."),
        Err(err) => {
            println!("\x1b[1;31mError:\x1b[0m {}", err);
            println!("\nPlease check `atg --help` for more options\n");
            process::exit(1);
        }
    }
}
