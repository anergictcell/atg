#[macro_use]
extern crate log;
use std::fs::File;
use std::process;

use bincode::{deserialize_from, serialize_into};
use clap::Parser;

use atglib::bed;
use atglib::fasta;
use atglib::fasta::FastaReader;
use atglib::genepred;
use atglib::genepredext;
use atglib::gtf;
use atglib::models::{GeneticCode, TranscriptWrite, Transcripts};
use atglib::qc;
use atglib::qc::QcCheck;
use atglib::read_transcripts;
use atglib::refgene;
use atglib::spliceai;
use atglib::utils::errors::AtgError;

mod cli;
use cli::{Args, InputFormat, OutputFormat};

mod reader_wrapper;
use reader_wrapper::ReadSeekWrapper;

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
    let fastareader = get_fasta_reader(&fasta_reference.as_deref());

    debug!("Writing transcripts as {} to {}", output_format, output_fd);

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
            writer.fasta_reader(fastareader?);
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
            writer.fasta_reader(fastareader?);
            writer.fasta_format(fasta_format.as_str());

            for tx in transcripts {
                let outfile = outdir.join(format!("{}.fasta", tx.name()));
                *writer.inner_mut() = std::io::BufWriter::new(File::create(outfile)?);
                writer.writeln_single_transcript(&tx)?;
            }
        }
        OutputFormat::FeatureSequence => {
            let mut writer = fasta::Writer::from_file(output_fd)?;
            writer.fasta_reader(fastareader?);
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
            add_genetic_code(&args.genetic_code, &mut writer)?;
            writer.fasta_reader(fastareader?);
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

/// Helper function to get a FastaReader that can read both local files and S3 objects
fn get_fasta_reader(filename: &Option<&str>) -> Result<FastaReader<ReadSeekWrapper>, AtgError> {
    if filename.is_none() {
        return Err(AtgError::new("no Fasta filename specified"));
    }
    // Both fasta_reader and fai_reader are Result<ReadSeekWrapper> instances
    let fasta_reader = ReadSeekWrapper::from_cli_arg(filename)?;
    let fai_reader = ReadSeekWrapper::from_filename(&format!("{}.fai", fasta_reader.filename()))?;

    Ok(FastaReader::from_reader(fasta_reader, fai_reader)?)
}

/// Attaches the chromosome-specific and default genetic code to the QC-Writer
fn add_genetic_code<W: std::io::Write, R: std::io::Read + std::io::Seek>(
    genetic_code_arg: &Vec<String>,
    writer: &mut qc::Writer<W, R>,
) -> Result<(), AtgError> {
    let codes = GeneticCodeSelecter::from_cli(genetic_code_arg)?;

    debug!("Setting default genetic code to {}", codes.default);
    writer.default_genetic_code(codes.default);

    for (chrom, code) in codes.custom {
        debug!("Adding genetic code {} for {}", &code, &chrom);
        writer.add_genetic_code(chrom, code);
    }
    Ok(())
}

#[derive(Default)]
/// Helper struct for parsing the genetic-code CLI arguments
///
/// The CLI argument can specify both one generic/default genetic code
/// and several chromosomse-specific genetic codes
struct GeneticCodeSelecter {
    default: GeneticCode,
    custom: Vec<(String, GeneticCode)>,
}

impl GeneticCodeSelecter {
    fn from_cli(genetic_code_arg: &Vec<String>) -> Result<Self, AtgError> {
        let mut code = GeneticCodeSelecter::default();
        for genetic_code_value in genetic_code_arg {
            match genetic_code_value.split_once(':') {
                // if the value contains a `:`, it is a key:value pair
                // for chromosome:genetic_code.
                Some((chrom, seq)) => {
                    let gen_code = GeneticCode::guess(seq)?;
                    debug!("Specified custom genetic code {} for {}", gen_code, chrom);
                    code.custom.push((chrom.to_string(), gen_code));
                }
                // Without `:` the genetic code is used as default
                None => {
                    let gen_code = GeneticCode::guess(genetic_code_value)?;
                    debug!("Specified default genetic code {}", gen_code);
                    code.default = gen_code;
                }
            }
        }
        Ok(code)
    }
}

/// Returns a filtered `Transcript`s object based on CLI-provided filter criteria
///
/// If a transcript fails one of the QC checks, it is removed from the output
///
/// Some QC checks might need the fasta file. To keep the logic simple,
/// the filter function will always run all QC checks (using `QcCheck`)
/// and then filter based on only the requested criteria.
/// This might not be the best performance approach, but other approaches
/// would add a lot more logic complexity.
/// The performance hit does not impact the most frequent use cases, where Fasta
/// data is needed anyway
fn filter_transcripts(transcripts: Transcripts, args: &Args) -> Result<Transcripts, AtgError> {
    let len_start = transcripts.len();

    let fasta_reference = &args.reference;
    let mut fastareader = get_fasta_reader(&fasta_reference.as_deref())?;

    // To collect all transcripts that pass the filter
    let mut filtered_transcripts = Transcripts::new();

    let codes = GeneticCodeSelecter::from_cli(&args.genetic_code)?;
    let mut custom_code: Option<&GeneticCode>;

    'tx_loop: for tx in transcripts.to_vec() {
        let qc = match codes.custom.is_empty() {
            true => QcCheck::new(&tx, &mut fastareader, &codes.default),
            false => {
                custom_code = None;
                for cc in &codes.custom {
                    if cc.0 == tx.chrom() {
                        custom_code = Some(&cc.1);
                        break;
                    }
                }
                QcCheck::new(&tx, &mut fastareader, custom_code.unwrap_or(&codes.default))
            }
        };

        for check in &args.qc_check {
            if check.remove(&qc) {
                debug!("Removing {} for failing QC filter {}", tx.name(), check);
                // Transcript fails the QC check, move on to the next transcript
                continue 'tx_loop;
            }
        }

        // only keep transcripts that did not fail any QC test
        filtered_transcripts.push(tx)
    }
    info!(
        "Filtered out {} transcripts.",
        len_start - filtered_transcripts.len()
    );
    Ok(filtered_transcripts)
}

fn main() {
    let cli_commands = Args::parse();

    loggerv::init_with_verbosity(cli_commands.verbose.into()).unwrap();

    let mut transcripts = match read_input_file(&cli_commands) {
        Ok(x) => x,
        Err(err) => {
            println!("\x1b[1;31mError:\x1b[0m {}", err);
            println!("\nPlease check `atg --help` for more options\n");
            process::exit(1);
        }
    };

    if !cli_commands.qc_check.is_empty() {
        debug!("Filtering transcripts");
        transcripts = match filter_transcripts(transcripts, &cli_commands) {
            Ok(t) => t,
            Err(err) => {
                println!("\x1b[1;31mError:\x1b[0m {}", err);
                println!("\nPlease check `atg --help` for more options\n");
                process::exit(1);
            }
        };
    }

    match write_output(&cli_commands, transcripts) {
        Ok(_) => debug!("All done here."),
        Err(err) => {
            println!("\x1b[1;31mError:\x1b[0m {}", err);
            println!("\nPlease check `atg --help` for more options\n");
            process::exit(1);
        }
    }
}
