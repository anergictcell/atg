# Changelog

# 0.6.0
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
