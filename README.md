# ATG
_ATG_ is a library to handle and convert different data formats used in Genomics and Transcriptomics. The library provides convenient APIs to parse GTF and RefGene data and work with the resulting transcripts for all kind of downstream analyses. 

The binary can be used to convert GTF into RefGene data and vice versa.

The main purpose is actually just that - convert a GTF file into a RefGene file. Surprsingly, there are not many tools to do this properly. Even _atg_ does not handle all edge cases of GTF - but I tried to handle as much as possible. In addition, transcripts can also be written in bed format as well.

The project started only because I wanted to learn Rust. You will see that some sections have really bad code, others will have some better and more improved code. Overall, I'm still very new to Rust and I'm sure I fell for many traps and use lots of unidiomatic code. I'm happy for any feedback and improvement suggestions.

The library is still in its infancy but works so far and can handle what it's supposed to do. The current API is probably going to change a lot in future updates, so be careful of using _atg_ in production or other critical workflows.

## Usage
### ATG command line tool

#### Install
The easiest way to install _ATG_ is to use `cargo` (if you have `cargo` and `rust` installed)
```bash
cargo install atg
```

### Download
You can download pre-build binaries for Linux and Mac (M1) of all recent releases from [Github](https://github.com/anergictcell/atg/releases)

#### Run
Convert a GTF file to a RefGene file
```bash
atg --from gtf --to refgene --input /path/to/input.gtf --output /path/to/output.refgene
```

Convert RefGene to GTF
```bash
atg --from refgene --to gtf --input /path/to/input.refgene --output /path/to/output.gtf
```

Convert RefGene to bed
```bash
atg --from refgene --to bed --input /path/to/input.refgene --output /path/to/output.bed
```

##### Other `--output` formats
###### gtf
Output in GTF format.

You can specify a `source` using the `--gtf-source`/`-g` option. Defaults to `atg`

###### fasta
cDNA sequence of all transcripts into one file. Please note that the sequence is stranded.

This target format requires a reference genome fasta file that must be specified using `-r <PATH>` or `--reference <PATH>`.

*This output allows different `--fasta-format` options:*
- `transcript` --> The full transcript sequence (from the genomic start to end position)
- `exons` --> The cDNA sequence of the processed transcript, i.e. the sequence of all exons, including non-coding exons.
- `cds` --> The CDS of the transcript

```text
>NM_007298.3 BRCA1
ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGC
TATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAAC
CTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAA
CTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGA
TATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTAGTCAACTTGTTG
...
>NM_001365057.2 C9orf85
ATGAGCTCCCAGAAAGGCAACGTGGCTCGTTCCAGACCTCAGAAGCACCA
GAATACGTTTAGCTTCAAAAATGACAAGTTCGATAAAAGTGTGCAGACCA
AGAAAATTAATGCAAAACTTCATGATGGAGTATGTCAGCGCTGTAAAGAA
GTTCTTGAGTGGCGTGTAAAATACAGCAAATACAAACCATTATCAAAACC
TAAAAAGTGA
...
```

###### fasta-split
Like `fasta` above, but one file for each transcript.

*This output also allows different `--fasta-format` options:*
- `transcript` --> The full transcript sequence (from the genomic start to end position)
- `exons` --> The cDNA sequence of the processed transcript, i.e. the sequence of all exons, including non-coding exons.
- `cds` --> The CDS of the transcript

###### feature-sequence
cDNA sequence of each feature (5' UTR, CDS, 3'UTR), each in a separate row.

**Example**:
```text
BRCA1   NM_007298.3 chr17   41196311    41197694    -   3UTR    CTGCAGCCAGCCAC...
BRCA1   NM_007298.3 chr17   41197694    41197819    -   CDS CAATTGGGCAGATGTGTG...
BRCA1   NM_007298.3 chr17   41199659    41199720    -   CDS GGTGTCCACCCAATTGTG...
BRCA1   NM_007298.3 chr17   41201137    41201211    -   CDS ATCAACTGGAATGGATGG...
BRCA1   NM_007298.3 chr17   41203079    41203134    -   CDS ATCTTCAGGGGGCTAGAA...
BRCA1   NM_007298.3 chr17   41209068    41209152    -   CDS CATGATTTTGAAGTCAGA...
BRCA1   NM_007298.3 chr17   41215349    41215390    -   CDS GGGTGACCCAGTCTATTA...
BRCA1   NM_007298.3 chr17   41215890    41215968    -   CDS ATGCTGAGTTTGTGTGTG...
BRCA1   NM_007298.3 chr17   41219624    41219712    -   CDS ATGCTCGTGTACAAGTTT...
BRCA1   NM_007298.3 chr17   41222944    41223255    -   CDS AGGGAACCCCTTACCTGG...
C9orf85 NM_001365057.2  chr9    74526555    74526650    +   5UTR    ATTGACAGAA...
C9orf85 NM_001365057.2  chr9    74526651    74526752    +   CDS ATGAGCTCCCAGAA...
C9orf85 NM_001365057.2  chr9    74561922    74562028    +   CDS AAAATTAATGCAAA...
C9orf85 NM_001365057.2  chr9    74597573    74597573    +   CDS A
C9orf85 NM_001365057.2  chr9    74597574    74600974    +   3UTR    TGGAGTCTCC...
```

###### raw
This is mainly useful for debugging, as it gives a quick glimpse into the Exons and CDS coordinates of the transcripts.


Yes, that's it. Nothing more at the moment.

Ok, you can also change the verbosity, by adding `-v` (show info messages), `-vv` (debug), `-vvv` (trace)

On most Linux systems, you can use `--input /dev/stdin` and/or `--output /dev/stdout` to pipe into and out of atg.

Of course, all commands also have shorthand parameters:
- `-f`, `--from`
- `-t`, `--to`
- `-i`, `--input`
- `-o`, `--output`

### ATG as library
[The library API is mostly documented inline and available on docs.rs](https://docs.rs/atg)

#### Examples

##### Convert GTF to RefGene
```no_run
use atg::gtf::Reader;
use atg::refgene::Writer;
use atg::models::{TranscriptRead, TranscriptWrite};

let mut reader = Reader::from_file("path/to/input.gtf")
    .unwrap_or_else(|_| panic!("Error opening input file."));

let mut writer = Writer::from_file("path/to/output.refgene")
    .unwrap_or_else(|_| panic!("Unable to open output file"));

let transcripts = reader.transcripts()
    .unwrap_or_else(|err| panic!("Error parsing GTF: {}", err));

match writer.write_transcripts(&transcripts) {
    Ok(_) => println!("Success"),
    Err(err) => panic!("Error writing RefGene file: {}", err)
};
```


## ToDo / Next tasks
- [x] Add to crates.io
- [x] Bed module to generate bed files with exons and introns
- [ ] Compare transcripts from two different inputs
- [x] Add fasta reading for nt and aa sequence outputs

## Known issues
### GTF parsing
- [ ] NM_001371720.1 has two book-ended exons (155160639-155161619 || 155161620-155162101). During input parsing, book-ended features are merged into one exon

