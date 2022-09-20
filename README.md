# ATG

Convert your genomic reference data between formats with a single tool. _ATG_ handles the conversion from and to GTF, GenePred(ext) and Refgene. You can generate bed files, fasta sequences or custom feature sequences. A single tool for all your conversion.

| File format | Can be used as source | Can be created |
| ----------- | ------------- | -------- |
| GTF | Yes | Yes |
| GenePred (extended) | Yes | Yes |
| RefGene | Yes | Yes |
| GenePred (simple) | No | Yes |
| Bed | No | Yes |
| Fasta | No | Yes (multiple options) |
| SpliceAI gene annotation | No | Yes |
| Quality Checks | No | Yes |


**Reasons to use _ATG_**
* No need to maintain multiple tools for one-way conversions (`gtfToGenePred`, `genePredToGtf`, etc). _ATG_ handles many formats and can convert in both directions.
* Speed: _ATG_ is really fast - almost twice as fast as `gtfToGenePred`.
* Robust parser: It handles GTF, GenePred with all extras according to spec.
* Low memory footprint: It also runs on machines with little RAM.
* Extra features, such as quality control and correctness checks.
* Open for contributions: Every help is welcome improve ATG or to add more functionality.
* You can also use _ATG_ as a library for your own Rust projects.


## ATG command line tool

### Install
There are currently 3 different options how to install _ATG_:

##### cargo
The easiest way to install _ATG_ is to use `cargo` (if you have `cargo` and `rust` installed)
```bash
cargo install atg
```

##### Pre-built binaries
You can download pre-built binaries for Linux and Mac from [Github](https://github.com/anergictcell/atg/releases).

##### From source
You can also build _ATG_ from source (if you have the rust toolchains installed):

```bash
git clone https://github.com/anergictcell/atg.git
cd atg
cargo build --release
```


### Usage
The main CLI arguments are 
- `-f`, `--from`: Specify the file format of the source (e.g. `gtf`, `genepredext`, `refgene`)
- `-t`, `--to`: Specify the target file format (e.g. `gtf`, `genepred`, `bed`, `fasta` etc)
- `-i`, `--input`: Path to source file. (Use `/dev/stdin` if you are using _atg_ in a pipe)
- `-o`, `--output`: Path to target file. Existing files will be overwritten. (Use `/dev/stdout` if you are using _atg_ in a pipe)
- `-v`, `-vv`, `-vvv`: Verbosity (info, debug, trace)
- `-h`, `--help`: Print the help dialog with detailed usage instructions.

Additional, optional arguments:
- `-g`, `--gtf-source`: Specify the source for GTF output files. Defaults to `atg`
- `-r`, `--reference`: Path of a reference genome fasta file. Required for fasta output
- `-c`, `--genetic-code`: Specify which genetic code to use for translating the transcripts. Genetic codes can be specified per chromosome by specifying the chromsome and the code, separated by `:` (e.g. `-c chrM:vertebrate mitochondrial`). They can also be specified for all chromsomes by omitting the chromosome (e.g. `-c vertebrate mitochondrial`). The argument can be specified multiple times (e.g: `-c "standard" -c "chrM:vertebrate mitochondrial" -c "chrAYN:alternative yeast nuclear"`). The code names are based on the `name` field from the [NCBI specs](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/gc.prt) but all lowercase characters. Alternatively, you can also specify the amino acid lookup table directly: `-c "chrM:FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"`. Defaults to `standard`.

#### Examples:
```bash
## Convert a GTF file to a RefGene file
atg --from gtf --to refgene --input /path/to/input.gtf --output /path/to/output.refgene

## Convert a GTF file to a GenePred file
atg --from gtf --to genepred --input /path/to/input.gtf --output /path/to/output.genepred

## Convert a GTF file to a GenePredExt file
atg --from gtf --to genepredext --input /path/to/input.gtf --output /path/to/output.genepredext

## Convert RefGene to GTF
atg --from refgene --to gtf --input /path/to/input.refgene --output /path/to/output.gtf

## Convert RefGene to bed
atg --from refgene --to bed --input /path/to/input.refgene --output /path/to/output.bed
```

### Supported `--output` formats

#### gtf
Output in [GTF](http://genome.ucsc.edu/FAQ/FAQformat.html#format4) format.

```text
chr9    ncbiRefSeq.2021-05-17   transcript  74526555    74600974    .   +   .   gene_id "C9orf85"; transcript_id "NM_001365057.2";
chr9    ncbiRefSeq.2021-05-17   exon        74526555    74526752    .   +   .   gene_id "C9orf85"; transcript_id "NM_001365057.2";
chr9    ncbiRefSeq.2021-05-17   5UTR        74526555    74526650    .   +   .   gene_id "C9orf85"; transcript_id "NM_001365057.2";
chr9    ncbiRefSeq.2021-05-17   CDS         74526651    74526752    .   +   0   gene_id "C9orf85"; transcript_id "NM_001365057.2";
chr9    ncbiRefSeq.2021-05-17   exon        74561922    74562028    .   +   .   gene_id "C9orf85"; transcript_id "NM_001365057.2";
chr9    ncbiRefSeq.2021-05-17   CDS         74561922    74562026    .   +   0   gene_id "C9orf85"; transcript_id "NM_001365057.2";
...
```

You can specify the value of the `source` column manually using the `--gtf-source`/`-g` option. Defaults to `atg`

#### refgene
Output in the [refGene](http://rohsdb.cmb.usc.edu/GBshape/cgi-bin/hgTables?hgsid=583_AkEae6dMkhjf5kd9BxNksFo9ySiK&hgta_doSchemaDb=mm10&hgta_doSchemaTable=refGene) format, as used by some UCSC and NCBI RefSeq services 

```text
0   NM_001101.5     chr7    -   5566778     5570232     5567378     5569288    6   5566778,5567634,5567911,5568791,5569165,5570154,    5567522,5567816,5568350,5569031,5569294,5570232,    0   ACTB    cmpl    cmpl    0,1,0,0,0,-1,
0   NM_001203247.2  chr7    -   148504474   148581383   148504737   148544390  20  148504474,148506162,148506401,148507424,148508716,148511050,148512005,148512597,148513775,148514313,148514968,148516687,148523560,148524255,148525831,148526819,148529725,148543561,148544273,148581255,    148504798,148506247,148506482,148507506,148508812,148511229,148512131,148512638,148513870,148514483,148515209,148516779,148523724,148524358,148525972,148526940,148529842,148543690,148544397,148581383,    0   EZH2    cmpl    cmpl    2,1,1,0,0,1,1,2,0,1,0,1,2,1,1,0,0,0,0,-1,
0   NM_001203248.2  chr7    -   148504474   148581383   148504737   148544390  20  148504474,148506162,148506401,148507424,148508716,148511050,148512005,148512597,148513775,148514313,148514968,148516687,148523560,148524255,148525831,148526819,148529725,148543588,148544273,148581255,    148504798,148506247,148506482,148507506,148508812,148511229,148512131,148512638,148513870,148514483,148515209,148516779,148523724,148524358,148525972,148526940,148529842,148543690,148544397,148581383,    0   EZH2    cmpl    cmpl    2,1,1,0,0,1,1,2,0,1,0,1,2,1,1,0,0,0,0,-1,
0   NM_001354750.2  chr11   +   113930432   114127487   113934022   114121277  7   113930432,113933932,114027058,114057673,114112888,114117919,114121047,  113930864,113935290,114027156,114057760,114113059,114118087,114127487,  0   ZBTB16  cmpl    cmpl    -1,0,2,1,1,1,1,
```

#### genepred(ext)
Output in the [GenePred(Ext)](http://genome.ucsc.edu/FAQ/FAQformat#format9) format, as used by some UCSC and NCBI RefSeq services 

**GenePred:**
```text
NM_001101.5     chr7    -   5566778     5570232     5567378     5569288     6   5566778,5567634,5567911,5568791,5569165,5570154,    5567522,5567816,5568350,5569031,5569294,5570232,
NM_001203247.2  chr7    -   148504474   148581383   148504737   148544390   20  148504474,148506162,148506401,148507424,148508716,148511050,148512005,148512597,148513775,148514313,148514968,148516687,148523560,148524255,148525831,148526819,148529725,148543561,148544273,148581255,    148504798,148506247,148506482,148507506,148508812,148511229,148512131,148512638,148513870,148514483,148515209,148516779,148523724,148524358,148525972,148526940,148529842,148543690,148544397,148581383,
NM_001203248.2  chr7    -   148504474   148581383   148504737   148544390   20  148504474,148506162,148506401,148507424,148508716,148511050,148512005,148512597,148513775,148514313,148514968,148516687,148523560,148524255,148525831,148526819,148529725,148543588,148544273,148581255,    148504798,148506247,148506482,148507506,148508812,148511229,148512131,148512638,148513870,148514483,148515209,148516779,148523724,148524358,148525972,148526940,148529842,148543690,148544397,148581383,
NM_001354750.2  chr11   +   113930432   114127487   113934022   114121277   7   113930432,113933932,114027058,114057673,114112888,114117919,114121047,  113930864,113935290,114027156,114057760,114113059,114118087,114127487,
```

**GenePredExt**
```text
NM_001101.5     chr7    -   5566778     5570232     5567378     5569288     6   5566778,5567634,5567911,5568791,5569165,5570154,    5567522,5567816,5568350,5569031,5569294,5570232,    0   ACTB    cmpl    cmpl    0,1,0,0,0,-1,
NM_001203247.2  chr7    -   148504474   148581383   148504737   148544390   20  148504474,148506162,148506401,148507424,148508716,148511050,148512005,148512597,148513775,148514313,148514968,148516687,148523560,148524255,148525831,148526819,148529725,148543561,148544273,148581255,    148504798,148506247,148506482,148507506,148508812,148511229,148512131,148512638,148513870,148514483,148515209,148516779,148523724,148524358,148525972,148526940,148529842,148543690,148544397,148581383,    0   EZH2    cmpl    cmpl    2,1,1,0,0,1,1,2,0,1,0,1,2,1,1,0,0,0,0,-1,
NM_001203248.2  chr7    -   148504474   148581383   148504737   148544390   20  148504474,148506162,148506401,148507424,148508716,148511050,148512005,148512597,148513775,148514313,148514968,148516687,148523560,148524255,148525831,148526819,148529725,148543588,148544273,148581255,    148504798,148506247,148506482,148507506,148508812,148511229,148512131,148512638,148513870,148514483,148515209,148516779,148523724,148524358,148525972,148526940,148529842,148543690,148544397,148581383,    0   EZH2    cmpl    cmpl    2,1,1,0,0,1,1,2,0,1,0,1,2,1,1,0,0,0,0,-1,
NM_001354750.2  chr11   +   113930432   114127487   113934022   114121277   7   113930432,113933932,114027058,114057673,114112888,114117919,114121047,  113930864,113935290,114027156,114057760,114113059,114118087,114127487,  0   ZBTB16  cmpl    cmpl    -1,0,2,1,1,1,1,
```

#### bed
Output in [bed](http://genome.ucsc.edu/FAQ/FAQformat#format1) format.

```text
chr7    5566778     5570232     ACTB:NM_001101.5       -   5567378    5569288    212,16,48   6   744,182,439,240,129,78  0,856,1133,2013,2387,3376
chr11   113930432   114127487   ZBTB16:NM_001354750.2  +   113934022  114121277  212,16,48   7   432,1358,98,87,171,168,6440 0,3500,96626,127241,182456,187487,190615
chr17   40852292    40897058    EZH1:NM_001321082.2    -   40854549   40880959   212,16,48   20  2318,85,81,82,96,179,126,41,92,197,181,92,164,103,177,121,129,128,91,30 0,2602,3465,4327,4813,5732,7683,8601,9571,12014,12934,17701,18179,18830,19998,22520,27360,28550,30553,44736
```

#### fasta
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

#### fasta-split
Like `fasta` above, but one file for each transcript. Instead of an output file, you must specify an output directory, _ATG_ will save each transcript as `<Transcript_name>.fasta`, e.g.: `NM_001365057.2.fasta`.

This target format requires a reference genome fasta file that must be specified using `--reference`/`-r`.

*This output allows different `--fasta-format` options:*
- `transcript`: The full transcript sequence (from the genomic start to end position, including introns)
- `exons`: The cDNA sequence of the processed transcript, i.e. the sequence of all exons, including non-coding exons.
- `cds` (default): The CDS of the transcript

#### feature-sequence
cDNA sequence of each feature (5' UTR, CDS, 3'UTR), each in a separate row.

This target format requires a reference genome fasta file that must be specified using `--reference`/`-r`.

```text
BRCA1   NM_007298.3     chr17   41196311    41197694    -   3UTR    CTGCAGCCAGCCAC...
BRCA1   NM_007298.3     chr17   41197694    41197819    -   CDS     CAATTGGGCAGATGTGTG...
BRCA1   NM_007298.3     chr17   41199659    41199720    -   CDS     GGTGTCCACCCAATTGTG...
BRCA1   NM_007298.3     chr17   41201137    41201211    -   CDS     ATCAACTGGAATGGATGG...
BRCA1   NM_007298.3     chr17   41203079    41203134    -   CDS     ATCTTCAGGGGGCTAGAA...
BRCA1   NM_007298.3     chr17   41209068    41209152    -   CDS     CATGATTTTGAAGTCAGA...
BRCA1   NM_007298.3     chr17   41215349    41215390    -   CDS     GGGTGACCCAGTCTATTA...
BRCA1   NM_007298.3     chr17   41215890    41215968    -   CDS     ATGCTGAGTTTGTGTGTG...
BRCA1   NM_007298.3     chr17   41219624    41219712    -   CDS     ATGCTCGTGTACAAGTTT...
BRCA1   NM_007298.3     chr17   41222944    41223255    -   CDS     AGGGAACCCCTTACCTGG...
C9orf85 NM_001365057.2  chr9    74526555    74526650    +   5UTR    ATTGACAGAA...
C9orf85 NM_001365057.2  chr9    74526651    74526752    +   CDS     ATGAGCTCCCAGAA...
C9orf85 NM_001365057.2  chr9    74561922    74562028    +   CDS     AAAATTAATGCAAA...
C9orf85 NM_001365057.2  chr9    74597573    74597573    +   CDS     A
C9orf85 NM_001365057.2  chr9    74597574    74600974    +   3UTR    TGGAGTCTCC...
```

#### spliceai
This is a custom format useful for [SpliceAI](https://github.com/Illumina/SpliceAI)
splice predictions. The repo lists [example files](https://github.com/Illumina/SpliceAI/tree/master/spliceai/annotations).
The output has one gene per row, each gene record contains a consensus transcript, created by merging overlapping exons.

```text
#NAME       CHROM   STRAND  TX_START    TX_END  EXON_START      EXON_END
OR4F5       1       +       69090       70008   69090,          70008,
AL627309.1  1       -       134900      139379  134900,137620,  135802,139379,
```

#### qc
Runs some basic consistency checks on the transcripts:

| QC check | Explanation | Non-Coding vs Coding | requires Fasta File |
| --- | --- | --- | --- |
| Exon | Contains at least one exon | all | no |
| Correct CDS Length | The length of the CDS is divisible by 3 | Coding | no |
| Correct Start Codon | The CDS starts with `ATG` | Coding | yes |
| Correct Stop Codon | The CDS ends with a Stop codon `TAG`, `TAA`, or `TGA` | Coding | yes |
| No upstream Start Codon | The 5'UTR does not contain another start codon `ATG` (This test do not make sense biologically. It is totally fine for a transcript to have upstream `ATG` start cordons that are not utilized but the ribosome.) | Coding | yes |
| No upstream Stop Codon| The CDS does not contain another in-frame stop-codon | Coding | yes |
| No Start codon | The full exon sequence does not contain a start codon `ATG` (Biologically speaking, a non-coding transcript could have `ATG` start codons that are not utilized) | Non-Coding | yes |
| Correct Coordinates | The transcript is within the coordinates of the reference genome | all | yes |

**Test results:**
- `NA` Test could not be performed (e.g. CDS-length for non-coding transcripts), so no conclusion could be drawn
- `OK` The test succeeded with an OK results
- `NOK` The test failed and gave a NOT OK result

```text
Gene     transcript     Exon  CDS Length  Correct Start Codon  Correct Stop Codon  No upstream Start Codon  No upstream Stop Codon  Correct Coordinates
FAM239A  NR_146581.1    OK    N/A         N/A                  N/A                 OK                       N/A                     OK
OR5H2    NM_001005482.1 OK    OK          OK                   OK                  OK                       OK                      OK
SNX20    NM_001144972.2 OK    OK          OK                   OK                  NOK                      OK                      OK
```

#### raw
This is mainly useful for debugging, as it gives a quick glimpse into the Exons and CDS coordinates of the transcripts.

#### bin
Save Transcripts in _ATG_ binary format for faster re-reading.


## ATG as library
_ATG_ uses the _atglib_ library, which is documented inline and available on [docs.rs](https://docs.rs/atglib)


## Known issues
### GTF parsing
- [ ] NM_001371720.1 has two book-ended exons (155160639-155161619 || 155161620-155162101). During input parsing, book-ended features are merged into one exon

