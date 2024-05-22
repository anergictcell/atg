# Changelog

## 0.8.5
- Update dependencies
- Include [update from `atglib`](https://github.com/anergictcell/atglib/pull/19) to fix bug in GTF parsing of `gene` records

## 0.8.4.
- Unpin `serde`, the issue from 0.8.3 has been fixed
- Update other dependencies, fix security issue in dependencies

## 0.8.3.
- Pin serde version (See https://github.com/serde-rs/serde/issues/2538 for context)

## 0.8.0
- Add QC-filter option
- Add option to specify custom genetic codes
- Allow using S3 objects for Fasta reference files (using [`S3Reader`](https://crates.io/crates/s3reader))
- Switch to [`ATGLib` 0.2.0](https://crates.io/crates/atglib)
- Switch to `clap` v4 and use `Derive`

## 0.7.0
- Switch to ATGlib 0.1.4
- Add `spliceai` ouput format

## 0.6.0
- Switch to ATGlib 0.1.3
    - ncExons in `feature-sequence` for non-coding transcripts
    - Add QC module

## 0.5.1
- Uncouple the lib from the CLI, make an extra atglib crate
- Fix exon sorting to take the transcirpt orientation into account

## 0.5
- Add exon numbers in GTF output
- Refactor GTF parsing to improve performance
- Add GenePredExt Reader and Writer
- Add GenePred Writer
- Remove chr prefixing during parsing of chromosomes/contigs
- Improve docs
- Add more integration tests

## 0.4
- Support writing Fasta files
- Support sequence output per feature
- Improved error handling

## 0.3
- Support writing bed files
- Update CICD and proper release guidelines

## 0.2
- Improve performance and memory usage
- Remove unused GTF attributes fields
- Improve logging

## 0.1
- Inital release
