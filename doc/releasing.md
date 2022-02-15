# Release instructions

## Checklist

- [ ] pipeline is green
- [ ] update `RELEASES.md`
- [ ] update `package.version` in `Cargo.toml`
- [ ] `cargo make pre_commit`
- [ ] commit (include Cargo.toml)
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