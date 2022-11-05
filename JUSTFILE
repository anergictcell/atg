set positional-arguments

new version:
    #!/usr/bin/env sh
    if git rev-parse --quiet --verify release/{{version}} > /dev/null; then
        echo "Release branch exists already"
    else
      echo "Creating release branch release/{{version}}"
      git checkout main && \
      git pull && \
      git checkout -b release/{{version}} && \
      sed -i .bck "s/^version =.*$/version = \"{{version}}\"/" ./Cargo.toml && \
      cargo build && \
      git commit -am "Prepare release branch {{version}}" && \
      git push -u origin release/{{version}}
    fi

@check version:
    git checkout release/{{version}} && git pull
    echo "Running linter and unittests"
    cargo clippy && cargo fmt && cargo test -q && cargo doc

@build version: (check version)
    cargo build --release --target=aarch64-apple-darwin
    cp target/aarch64-apple-darwin/release/atg target/atg-{{version}}-aarch64-apple-darwin
    cargo build --release --target=x86_64-apple-darwin
    cp target/x86_64-apple-darwin/release/atg target/atg-{{version}}-x86_64-apple-darwin
    cargo build --release --target=x86_64-unknown-linux-gnu
    cp target/x86_64-unknown-linux-gnu/release/atg target/atg-{{version}}-x86_64-unknown-linux-gnu
    cargo build --release --target=x86_64-unknown-linux-musl
    cp target/x86_64-unknown-linux-musl/release/atg target/atg-{{version}}-x86_64-unknown-linux-musl

@release version: (build version)
    git tag {{version}}
    git push --tags
    cargo publish

@testbuild version:
    cargo build --release --target=aarch64-apple-darwin
    cp target/aarch64-apple-darwin/release/atg target/atg-{{version}}-aarch64-apple-darwin
    cargo build --release --target=x86_64-unknown-linux-gnu
    cp target/x86_64-unknown-linux-gnu/release/atg target/atg-{{version}}-x86_64-unknown-linux-gnu
    cargo build --release --target=x86_64-unknown-linux-musl
    cp target/x86_64-unknown-linux-musl/release/atg target/atg-{{version}}-x86_64-unknown-linux-musl

test:
    #!/usr/bin/env zsh
    echo -ne "Checking formatting and doc generation"
    (cargo clippy -q && cargo fmt -q --check && cargo doc -q && \
    echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"
    echo -ne "Checking GTF to RefGene"
    (diff <( cargo run -q -- -f gtf -i tests/data/example.gtf -t refgene -o /dev/stdout 2> /dev/null | sort ) tests/data/example.refgene && \
    echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"
    echo -ne "Checking RefGene to GTF"
    (diff <( cargo run -q -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | cut -f1-5,7-8 | sort ) <( cut -f1-5,7-8 tests/data/example.gtf) && \
    echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"
    echo -ne "Checking GTF to SpliceAI"
    (diff <( cargo run -q -- -f gtf -i tests/data/example.gtf -t spliceai | wc -l | sed "s/ //g" ) <(echo "7")  && \
    echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"

    echo -ne "Checking uniqness of GTF attribute column across exons"
    (diff <( cargo run -q -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\texon\t" | cut -f9 | uniq -c | awk '{print $1}' | sort | uniq -c) <(echo " 695 1") && \
    echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"
    echo -ne "Checking exon numbering in exon, CDS, and UTR"
    (diff <( cargo run -q -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\texon\t" | grep "NM_004015.3.18" | wc -l | sed "s/ //g") <(echo "1") && \
    diff <( cargo run -q -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\tCDS\t" | grep "NM_004015.3.18" | wc -l | sed "s/ //g") <(echo "1") && \
    diff <( cargo run -q -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\t5UTR\t" | grep "NM_004015.3.18" | wc -l | sed "s/ //g") <(echo "0") && \
    diff <( cargo run -q -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\t5UTR\t" | grep "NM_004015.3.0" | wc -l | sed "s/ //g") <(echo "0") && \
    diff <( cargo run -q -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\t5UTR\t" | grep "NM_004015.3.1" | wc -l | sed "s/ //g") <(echo "1") && \
    diff <( cargo run -q -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\t5UTR\t" | grep "NM_004015.3.2" | wc -l | sed "s/ //g") <(echo "0") && \
    diff <( cargo run -q -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\t3UTR\t" | grep "NM_004015.3.18" | wc -l | sed "s/ //g") <(echo "1") && \
    diff <( cargo run -q -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "NM_004015.3.18" | wc -l | sed "s/ //g") <(echo "3") && \
    echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"
    echo -ne "Checking QC stats"
    # This test only checks that QC passes without errors. The transcripts cannot be checked with the small.fasta file
    (diff <( cargo run -q -- -f gtf -i tests/data/example.gtf -r tests/data/small.fasta -t qc | cut -f 3-9 | grep "N/A\tN/A\tN/A\tN/A\tN/A\tN/A\tNOK" | wc -l | sed "s/ //g") <(echo "27") && \
    echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"




    if [ ! -f tests/data/hg19.fasta ]; then
        echo "\e[35m\e[1mSkipping tests that require a reference genome Fasta file\e[0m" >&2
    else

        echo -ne "Fasta CDS"
        (diff <( cargo run -q -- -f gtf -t fasta -r tests/data/hg19.fasta -i tests/data/example.gtf | wc -l | sed "s/ //g") <(echo "1880") && \
        diff <( cargo run -q -- -f gtf -t fasta -r tests/data/hg19.fasta -i tests/data/example.gtf --fasta-format cds | wc -l | sed "s/ //g") <(echo "1880") && \
        echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"

        echo -ne "Fasta Exons"
        (diff <( cargo run -q -- -f gtf -t fasta -r tests/data/hg19.fasta -i tests/data/example.gtf --fasta-format exons | wc -l | sed "s/ //g") <(echo "3322") && \
        echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"

        echo -ne "Fasta Full transcript"
        (diff <( cargo run -q -- -f gtf -t fasta -r tests/data/hg19.fasta -i tests/data/example.gtf --fasta-format transcript | wc -l | sed "s/ //g") <(echo "258605") && \
        echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"

        echo -ne "Feature Sequence"
        (diff <( cargo run -q -- -f gtf -t feature-sequence -r tests/data/hg19.fasta -i tests/data/example.gtf | wc -l | sed "s/ //g") <(echo "749") && \
        features=$(cargo run -q -- -f gtf -t feature-sequence -r tests/data/hg19.fasta -i tests/data/example.gtf | cut -f 7 | sort | uniq -c) && \
        diff <(echo $features | grep "3UTR" | sed "s/ //g") <(echo "273UTR") && \
        diff <(echo $features | grep "5UTR" | sed "s/ //g") <(echo "615UTR") && \
        diff <(echo $features | grep "CDS" | sed "s/ //g") <(echo "661CDS") && \
        echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"

        echo -ne "Test QC"
        (diff <( cargo run -q -- -f gtf -i tests/data/example.gtf -r tests/data/hg19.fasta -t qc | tail -n +2  | cut -f3-9 | sort |  uniq -c | sed "s/ //g") <(echo -ne "14OK\tOK\tOK\tOK\tNOK\tOK\tOK\n13OK\tOK\tOK\tOK\tOK\tOK\tOK\n") && \
        echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"

        echo -ne "Test Filter"
        (diff <( cargo run -q -- -f gtf -i tests/data/example.gtf -t refgene  -r tests/data/hg19.fasta -q exon -q start -q stop -q exon -q cds-length -q upstream-start -q upstream-stop -q coordinates | wc -l | sed "s/ //g") <(echo -ne "13\n") && \
        diff <( cargo run -q -- -f gtf -i tests/data/example.gtf -t refgene  -r tests/data/hg19.fasta -q exon -q start -q stop -q exon -q cds-length -q upstream-stop -q coordinates | wc -l | sed "s/ //g") <(echo -ne "27\n") && \
        diff <( cargo run -q -- -f gtf -i tests/data/example.gtf -t refgene  -r tests/data/hg19.fasta -q upstream-stop -c "chrY:vertebrate mitochondrial" | wc -l | sed "s/ //g") <(echo -ne "26\n") && \
        echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"

    fi

benchmark name:
    cargo build --release --target=aarch64-apple-darwin
    rm target/atg-{{name}}-aarch64-apple-darwin
    cp target/aarch64-apple-darwin/release/atg target/atg-{{name}}-aarch64-apple-darwin
    time target/atg-{{name}}-aarch64-apple-darwin -f gtf -i tests/data/hg19.ncbiRefSeq.gtf -o /dev/null -t none
    flamegraph --root -- target/atg-{{name}}-aarch64-apple-darwin -f gtf -i tests/data/hg19.ncbiRefSeq.gtf -o /dev/null -t none