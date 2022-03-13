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
    cargo build --release --target=x86_64-unknown-linux-gnu
    cp target/x86_64-unknown-linux-gnu/release/atg target/atg-{{version}}-x86_64-unknown-linux-gnu
    cargo build --release --target=x86_64-unknown-linux-musl
    cp target/x86_64-unknown-linux-musl/release/atg target/atg-{{version}}-x86_64-unknown-linux-musl

@release version: (build version)
    git tag {{version}}
    git push --tags
    cargo publish

@betabuild version:
    cargo build --release --target=aarch64-apple-darwin
    cp target/aarch64-apple-darwin/release/atg target/atg-{{version}}-beta-aarch64-apple-darwin
    cargo build --release --target=x86_64-unknown-linux-gnu
    cp target/x86_64-unknown-linux-gnu/release/atg target/atg-{{version}}-beta-x86_64-unknown-linux-gnu
    cargo build --release --target=x86_64-unknown-linux-musl
    cp target/x86_64-unknown-linux-musl/release/atg target/atg-{{version}}-beta-x86_64-unknown-linux-musl

test:
    #!/usr/bin/env zsh
    cargo clippy && cargo fmt && cargo test -q && cargo doc
    echo -ne "Checking GTF to RefGene"
    diff <( cargo run -- -f gtf -i tests/data/example.gtf -t refgene -o /dev/stdout 2> /dev/null | sort ) tests/data/example.refgene
    echo " \e[32m\e[1mOK\e[0m"
    echo -ne "Checking RefGene to GTF"
    diff <( cargo run -- -f refgene -i tests/data/example.refgene -t gtf -g ncbiRefSeq.2021-05-17 -o /dev/stdout 2> /dev/null | cut -f1-5,7-8 | sort ) <( cut -f1-5,7-8 tests/data/example.gtf)
    echo " \e[32m\e[1mOK\e[0m"
