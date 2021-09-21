# ATG
_atg_ is a library to handle and convert different data formats used in Genomics and Transcriptomics. The library provides convenient APIs to parse GTF and RefGene data and work with the resulting transcripts for all kind of downstream analyses. The binary can be used to convert GTF into RefGene data and vice versa.

The main purpose is actually just that - convert a GTF file into a RefGene file. Surprsingly, there are not many tools to do this properly. Even _atg_ does not handle all edge cases of GTF - but I tried to handle as much as possible.

The project started only because I wanted to learn Rust. You will see that some sections have really bad code, others will have some better and more improved code. Overall, I'm still very new to Rust and I'm sure I fell for many traps and use lots of unidiomatic code. I'm happy for any feedback and improvement suggestions.

The library is still in its infancy but works so far and can handle what it's supposed to do. The current API is probably going to change a lot in future updates, so be careful of using _atg_ in production or other critical workflows.

## Usage
### Binary atg tool

Convert a GTF file to a RefGene file
```bash
atg --from gtf --to refgene --input /path/to/input.gtf --output /path/to/output.refgene
```

Convert RefGene to GTF
```bash
atg --from refgene --to gtf --input /path/to/input.refgene --output /path/to/output.gtf
```

Yes, that's it. Nothing more at the moment.

Ok, you can also change the verbosity, by adding `-v` (show info messages), `-vv` (debug), `-vvv` (trace)

On most Linux systems, you can use `--input /dev/stdin` and/or `--output /dev/stdout` to pipe into and out of atg.

Of course, all commands also have shorthand parameters:
- `-f`, `--from`
- `-t`, `--to`
- `-i`, `--input`
- `-o`, `--output`

### ATG as library
The library API is mostly documented inline

#### models
The models mod contains the main data models and structs - `Transcripts`, `Transcript` and `Exon`. In addition, there are some other important structs and functions.

#### utils
Contains utility functions and some shared Error types.

#### gtf
Handles GTF parsing (Reader) and composing (Writer).

#### refgene
Handled RefGene parsing (Reader) and composing (Writer())

## Known issues
### GTF parsing
- [ ] NM_00137172 has two book-ended exons (155160639-155161619 || 155161620-155162101). During input parsing, book-ended features are merged into one exon


## ToDo
- [ ] Generate bed files with exons and introns
- [ ] Compare transcripts from two different inputs
- [ ] Add fasta reading from BioRust for nt and aa sequence outputs