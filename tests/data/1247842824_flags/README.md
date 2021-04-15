# Generated with

```bash
cotter \
  -m /mnt/data/1247842824_vis/1247842824.metafits \
  -o /mnt/data/1247842824_flags/Flagfile%%.mwaf \
  -nostats \
  -noantennapruning \
  -noflagautos \
  -noflagdcchannels \
  -norfi \
  -nogeom \
  -nosbgains \
  -edgewidth 0 \
  -allowmissing \
  -nocablelength \
  -flag-strategy /usr/local/share/aoflagger/strategies/generic-minimal.lua \
  /mnt/data/1247842824_vis/1247842824_*gpubox*.fits
```
