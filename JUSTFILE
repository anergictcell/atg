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
    (cargo clippy && cargo fmt --check && cargo test -q && cargo doc && \
    echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"
    echo -ne "Checking GTF to RefGene"
    (diff <( cargo run -- -f gtf -i tests/data/example.gtf -t refgene -o /dev/stdout 2> /dev/null | sort ) tests/data/example.refgene && \
    echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"
    echo -ne "Checking RefGene to GTF"
    (diff <( cargo run -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | cut -f1-5,7-8 | sort ) <( cut -f1-5,7-8 tests/data/example.gtf) && \
    echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"
    echo -ne "Checking uniqness of GTF attribute column across exons"
    (diff <( cargo run -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\texon\t" | cut -f9 | uniq -c | awk '{print $1}' | sort | uniq -c) <(echo " 695 1") && \
    echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"
    echo -ne "Checking exon numbering in exon, CDS, and UTR"
    (diff <( cargo run -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\texon\t" | grep "NM_004015.3.18" | wc -l | sed "s/ //g") <(echo "1") && \
    diff <( cargo run -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\tCDS\t" | grep "NM_004015.3.18" | wc -l | sed "s/ //g") <(echo "1") && \
    diff <( cargo run -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\t5UTR\t" | grep "NM_004015.3.18" | wc -l | sed "s/ //g") <(echo "0") && \
    diff <( cargo run -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\t5UTR\t" | grep "NM_004015.3.0" | wc -l | sed "s/ //g") <(echo "0") && \
    diff <( cargo run -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\t5UTR\t" | grep "NM_004015.3.1" | wc -l | sed "s/ //g") <(echo "1") && \
    diff <( cargo run -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\t5UTR\t" | grep "NM_004015.3.2" | wc -l | sed "s/ //g") <(echo "0") && \
    diff <( cargo run -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "\t3UTR\t" | grep "NM_004015.3.18" | wc -l | sed "s/ //g") <(echo "1") && \
    diff <( cargo run -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | grep "NM_004015.3.18" | wc -l | sed "s/ //g") <(echo "3") && \
    echo " \e[32m\e[1mOK\e[0m") || echo "\e[31m\e[1mERROR\e[0m"

benchmark name:
    cargo build --release --target=aarch64-apple-darwin
    rm target/atg-{{name}}-aarch64-apple-darwin
    cp target/aarch64-apple-darwin/release/atg target/atg-{{name}}-aarch64-apple-darwin
    time target/atg-{{name}}-aarch64-apple-darwin -f gtf -i tests/data/hg19.ncbiRefSeq.gtf -o /dev/null -t none
    flamegraph --root -- target/atg-{{name}}-aarch64-apple-darwin -f gtf -i tests/data/hg19.ncbiRefSeq.gtf -o /dev/null -t none