# Release instructions

## Checklist

- [ ] pipeline is green
- [ ] everything has been merged into main
- [ ] update `RELEASES.md`
- [ ] update `package.version` in `Cargo.toml`
- [ ] `cargo make pre_commit`
- [ ] commit (include Cargo.toml)
- [ ] you're absolute sure everything has been merged into main
- [ ] `git tag -a $tag -m $tag`
- [ ] `git push`
- [ ] `git push --tags`

## Template

```md
# Version X.Y.Z (YYYY-MM-DD)
- ✨ new features:
- 🙏 quality of life:
- 🏗 api changes:
- 🐛 bug fixes:
- ➕ dependencies:
```