# Release Checklist

- `VERSION=<new version>`
- Run unittests, clippy and linting
  - `cargo test && cargo clippy && cargo fmt`
- Update version in [`Cargo.toml`](./Cargo.toml).
- Run `cargo build` to update [`Cargo.lock`](./Cargo.lock).
- Add changes since last release to [`CHANGELOG.md`](./CHANGELOG.md). (You
  should do this with every commit!)
- Commit all changes with commit message: `vX.Y.Z Release`
- Tag commit and push it to GitHub: `git tag $VERSION && git push origin $VERSION`
- Publish new version to crates.io: `cargo publish`
- Generate new binaries:
  - `cargo build --release`  # MacOS M1 binary
  - `cargo build --release --target=x86_64-unknown-linux-gnu`
  - `cargo build --release --target=x86_64-unknown-linux-musl`
- Create GitHub release
  - Click "Create new release"
  - Select tag
  - Copy stuff from `CHANGELOG.md` to description
  - Attach binaries generated above
