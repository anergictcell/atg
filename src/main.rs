#[macro_use]
extern crate log;

use clap::{App, Arg};

use atg::gtf;
use atg::models::TranscriptWrite;
use atg::read_transcripts;
use atg::refgene;

fn main() {
    let cli_commands = App::new("Transcripts")
        .version("0.1.0")
        .author("Jonas Marcello <jonas.marcello@centogene.com>")
        .about("Convert transcript data from and to different file formats")
        .arg(
            Arg::with_name("from")
                .short("f")
                .long("from")
                .value_name("refgene|gtf")
                .help("Defines the input data format")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("to")
                .short("t")
                .long("to")
                .value_name("refgene|gtf")
                .help("Defines the output data format")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("input")
                .short("i")
                .long("input")
                .value_name("/path/to/input/fike")
                .help("Path to input file")
                .takes_value(true)
                .default_value("/dev/stdin"),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .value_name("/path/to/output/fike")
                .help("Path to output file")
                .takes_value(true)
                .default_value("/dev/stdout"),
        )
        .arg(
            Arg::with_name("v")
                .short("v")
                .multiple(true)
                .help("Sets the level of verbosity"),
        )
        .get_matches();

    loggerv::init_with_verbosity(cli_commands.occurrences_of("v")).unwrap();

    let input_fd = cli_commands.value_of("input").unwrap();
    let output_fd = cli_commands.value_of("output").unwrap();

    let input_format = cli_commands.value_of("from").unwrap();
    let output_format = cli_commands.value_of("to").unwrap();

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
                writer.set_source("ncbiRefSeq.2021-05-17");
                writer.write_transcripts(&transcripts)
            }
            Err(err) => panic!("Error writing GTF: {}", err),
        },
        _ => panic!("Invalid >>to<< parameter"),
    };
    debug!("Finshed writing output data");
}
