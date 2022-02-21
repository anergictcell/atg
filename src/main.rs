#[macro_use]
extern crate log;

use clap::{App, Arg, ArgMatches};

use std::process;

use atg::bed;
use atg::gtf;
use atg::models::TranscriptWrite;
use atg::read_transcripts;
use atg::refgene;

fn parse_cli_args() -> ArgMatches<'static> {
    App::new(env!("CARGO_PKG_NAME"))
        .version(atg::VERSION)
        .author(env!("CARGO_PKG_AUTHORS"))
        .about("Convert transcript data from and to different file formats")
        .after_help(env!("CARGO_PKG_HOMEPAGE"))
        .arg(
            Arg::with_name("from")
                .short("f")
                .long("from")
                .possible_values(&["refgene", "gtf"])
                .case_insensitive(true)
                .value_name("file-format")
                .help("Defines the input data format")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("to")
                .short("t")
                .long("to")
                .possible_values(&["refgene", "gtf", "bed", "fasta"])
                .case_insensitive(true)
                .value_name("file-format")
                .help("Defines the output data format")
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
            Arg::with_name("fasta_format")
                .long("fasta-format")
                .possible_values(&["transcript", "exons", "cds"])
                .value_name("format")
                .help(
                    "Which part of the transcript to translate. \
                       (This open is only needed when generating \
                       fasta output)",
                )
                .next_line_help(true)
                .display_order(1000)
                .takes_value(true)
                .default_value("cds")
                .required_if("to", "fasta"),
        )
        .arg(
            Arg::with_name("v")
                .short("v")
                .multiple(true)
                .help("Sets the level of verbosity"),
        )
        .get_matches()
}

fn main() {
    let cli_commands = parse_cli_args();

    loggerv::init_with_verbosity(cli_commands.occurrences_of("v")).unwrap();

    let input_fd = cli_commands.value_of("input").unwrap();
    let output_fd = cli_commands.value_of("output").unwrap();

    let input_format = cli_commands.value_of("from").unwrap();
    let output_format = cli_commands.value_of("to").unwrap();
    let fasta_format = cli_commands.value_of("fasta_format").unwrap();

    debug!("pid is {}", process::id());

    debug!("Parsed CLI arguments. Starting processing");

    trace!("Parameters:");
    trace!(
        "Input Format: {} | Input Source: {}",
        input_format,
        input_fd
    );
    trace!(
        "Output Format: {} | Output Target: {}",
        output_format,
        output_fd
    );

    debug!("Reading input data");

    let transcripts = match input_format {
        "refgene" => read_transcripts(refgene::Reader::from_file(input_fd)),
        "gtf" => read_transcripts(gtf::Reader::from_file(input_fd)),
        _ => panic!("Invalid command"),
    }
    .unwrap_or_else(|e| panic!("Error parsing the input data: {}", e));

    debug!(
        "Finished parsing input data. Found {} transcripts",
        transcripts.len()
    );

    debug!("Writing ouput data");
    let _res = match output_format {
        "refgene" => match refgene::Writer::from_file(output_fd) {
            Ok(mut writer) => writer.write_transcripts(&transcripts),
            Err(err) => panic!("Error writing GTF: {}", err),
        },
        "gtf" => match gtf::Writer::from_file(output_fd) {
            Ok(mut writer) => {
                writer.set_source("atg");
                writer.write_transcripts(&transcripts)
            }
            Err(err) => panic!("Error writing GTF: {}", err),
        },
        "bed" => match bed::Writer::from_file(output_fd) {
            Ok(mut writer) => writer.write_transcripts(&transcripts),
            Err(err) => panic!("Error writing GTF: {}", err),
        },
        "fasta" => {
            println!(
                "Fasta ourput is not yet implemented. Target {}",
                fasta_format
            );
            Ok(())
        }
        "raw" => {
            for t in transcripts {
                println!("{}", t);
                for exon in t.exons() {
                    println!("{}", exon)
                }
            }
            Ok(())
        }
        _ => panic!("Invalid >>to<< parameter"),
    };
    debug!("Finshed writing output data");
}
