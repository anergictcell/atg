# ATG
_ATG_ is a library and standalone CLI tool to handle and convert different data formats used in Genomics and Transcriptomics. The library provides convenient APIs to parse GTF and RefGene data and work with the resulting transcripts for all kind of downstream analyses.

The binary can be used to convert between GTF and RefGene data, generate bed files for transcripts or generate the nucleotide sequence as Fasta.

The main purpose is actually just that - convert between GTF and RefGene file formats. Surprsingly, there are not many tools to do this properly. Even _atg_ does not handle all edge cases of GTF - but I tried to handle as many as possible. In addition, transcripts can also be written in bed format or as fasta.

The project started only because I wanted to learn Rust. You will see that some sections have really bad code, others will have some better and more improved code. Overall, I'm still very new to Rust and I'm sure I fell for many traps and use lots of unidiomatic code. I'm happy for any feedback and improvement suggestions.

The library is still in its infancy but works so far and can handle what it's supposed to do. The current API is probably going to change a lot in future updates, so be careful of using _atg_ in production or other critical workflows.

## Usage
### ATG command line tool

#### Install
##### cargo
The easiest way to install _ATG_ is to use `cargo` (if you have `cargo` and `rust` installed)
```bash
cargo install atg
```

##### Pre-built binaries
You can download pre-built binaries for Linux and Mac (M1) from [Github](https://github.com/anergictcell/atg/releases). Save the downloaded binaries in a folder in your `PATH` or create a symlink as `atg`.

##### From source
You can also build _ATG_ from source (if you have the rust toolchains installed):

```bash
git clone https://github.com/anergictcell/atg.git
cd atg
cargo build --release
````

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

##### Supported `--output` formats

###### gtf
Output in GTF format.

```text
chr9    ncbiRefSeq.2021-05-17   transcript  74526555    74600974    .   +   .   gene_id "C9orf85"; transcript_id "NM_001365057.2";  gene_name "C9orf85";
chr9    ncbiRefSeq.2021-05-17   exon    74526555    74526752    .   +   .   gene_id "C9orf85"; transcript_id "NM_001365057.2"; exon_number "1"; exon_id "NM_001365057.2.1"; gene_name "C9orf85";
chr9    ncbiRefSeq.2021-05-17   5UTR    74526555    74526650    .   +   .   gene_id "C9orf85"; transcript_id "NM_001365057.2"; exon_number "1"; exon_id "NM_001365057.2.1"; gene_name "C9orf85";
chr9    ncbiRefSeq.2021-05-17   CDS 74526651    74526752    .   +   0   gene_id "C9orf85"; transcript_id "NM_001365057.2"; exon_number "1"; exon_id "NM_001365057.2.1"; gene_name "C9orf85";
chr9    ncbiRefSeq.2021-05-17   exon    74561922    74562028    .   +   .   gene_id "C9orf85"; transcript_id "NM_001365057.2"; exon_number "2"; exon_id "NM_001365057.2.2"; gene_name "C9orf85";
chr9    ncbiRefSeq.2021-05-17   CDS 74561922    74562026    .   +   0   gene_id "C9orf85"; transcript_id "NM_001365057.2"; exon_number "2"; exon_id "NM_001365057.2.2"; gene_name "C9orf85";
...
```

You can specify the value of the `source` column manually using the `--gtf-source`/`-g` option. Defaults to `atg`

###### refgene
Output in the refGene format, as used by some UCSC and NCBI RefSeq services 

```text
0   NM_001101.5 chr7    -   5566778 5570232 5567378 5569288 6   5566778,5567634,5567911,5568791,5569165,5570154,    5567522,5567816,5568350,5569031,5569294,5570232,    0   ACTB    cmpl    cmpl    0,1,0,0,0,-1,
0   NM_001203247.2  chr7    -   148504474   148581383   148504737   148544390   20  148504474,148506162,148506401,148507424,148508716,148511050,148512005,148512597,148513775,148514313,148514968,148516687,148523560,148524255,148525831,148526819,148529725,148543561,148544273,148581255,    148504798,148506247,148506482,148507506,148508812,148511229,148512131,148512638,148513870,148514483,148515209,148516779,148523724,148524358,148525972,148526940,148529842,148543690,148544397,148581383,    0   EZH2    cmpl    cmpl    2,1,1,0,0,1,1,2,0,1,0,1,2,1,1,0,0,0,0,-1,
0   NM_001203248.2  chr7    -   148504474   148581383   148504737   148544390   20  148504474,148506162,148506401,148507424,148508716,148511050,148512005,148512597,148513775,148514313,148514968,148516687,148523560,148524255,148525831,148526819,148529725,148543588,148544273,148581255,    148504798,148506247,148506482,148507506,148508812,148511229,148512131,148512638,148513870,148514483,148515209,148516779,148523724,148524358,148525972,148526940,148529842,148543690,148544397,148581383,    0   EZH2    cmpl    cmpl    2,1,1,0,0,1,1,2,0,1,0,1,2,1,1,0,0,0,0,-1,
0   NM_001354750.2  chr11   +   113930432   114127487   113934022   114121277   7   113930432,113933932,114027058,114057673,114112888,114117919,114121047,  113930864,113935290,114027156,114057760,114113059,114118087,114127487,  0   ZBTB16  cmpl    cmpl    -1,0,2,1,1,1,1,
```

###### fasta
Writes the cDNA sequence of all transcripts into one file. Please note that the sequence is stranded.

This target format requires a reference genome fasta file that must be specified using `--reference`/`-r`.

*This output allows different `--fasta-format` options:*
- `transcript`: The full transcript sequence (from the genomic start to end position, including introns)
- `exons`: The cDNA sequence of the processed transcript, i.e. the sequence of all exons, including non-coding exons.
- `cds` (default): The CDS of the transcript

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
Like `fasta` above, but one file for each transcript. Instead of an output file, you must specify an output directory, _ATG_ will save each transcript as `<Transcript_name>.fasta`, e.g.: `NM_001365057.2.fasta`.

This target format requires a reference genome fasta file that must be specified using `--reference`/`-r`.

*This output allows different `--fasta-format` options:*
- `transcript`: The full transcript sequence (from the genomic start to end position, including introns)
- `exons`: The cDNA sequence of the processed transcript, i.e. the sequence of all exons, including non-coding exons.
- `cds` (default): The CDS of the transcript

###### feature-sequence
cDNA sequence of each feature (5' UTR, CDS, 3'UTR), each in a separate row.

This target format requires a reference genome fasta file that must be specified using `--reference`/`-r`.

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

###### bin
Save Transcripts in _ATG_ binary format for faster re-reading.


#### Tips
Reading in GTF files is rather slow, due to the complexity of the format. If you need to repeatedly read in data from GTF, I recommend to generate a RefGene or binary file once and use this as input for subsequent steps.

You can change the verbosity, by adding `-v` (show info messages), `-vv` (debug), `-vvv` (trace)

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
- [x] Binary data format
- [ ] use [Smartstring](https://crates.io/crates/smartstring) or [Smallstr](https://crates.io/crates/smallstr) for gene-symbol, transcript name and chromosome
- [ ] Remove function to sanitize chromosome name with `chr` prefix
- [ ] Parallelize input parsing
- [ ] Check if exons can be stored in smaller vec
- [ ] Use std::mem::replace to move out of attributes, e.g. in TranscriptBuilder and remove Copy/Clone traits <https://stackoverflow.com/questions/31307680/how-to-move-one-field-out-of-a-struct-that-implements-drop-trait>

## Known issues
### GTF parsing
- [ ] NM_001371720.1 has two book-ended exons (155160639-155161619 || 155161620-155162101). During input parsing, book-ended features are merged into one exon

