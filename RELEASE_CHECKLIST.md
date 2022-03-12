# Release Checklist

- Crate a new release branch `just new <VERSION>`, e.g. `just new 1.0.0`
- Merge all feature branches into release branch
- Add changes since last release to [`CHANGELOG.md`](./CHANGELOG.md). (You
  should do this with every commit!)
- (Optional): Run unittests and checks: `just check <VERSION>`
- (Optional): Build binaries: `just build <VERSION>`
- Release to crates.io (This step run tests and build steps from above): `just release <VERSION>`
- Create GitHub release
  - Click "Create new release"
  - Select tag
  - Copy stuff from `CHANGELOG.md` to description
  - Attach binaries generated above
