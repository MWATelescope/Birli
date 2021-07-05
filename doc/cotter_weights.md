

## Generate uvfits

```bash
cotter \
  -m /mnt/data/1247842824_vis/1247842824.metafits \
  -o "/mnt/data/1247842824_vis/1247842824.uvfits" \
  -allowmissing \
  -edgewidth 0 \
  -endflag 0 \
  -initflag 0 \
  -noantennapruning \
  -nocablelength \
  -noflagautos \
  -noflagdcchannels \
  -nogeom \
  -sbpassband tests/data/subband-passband-128ch-unitary.txt \
  -nostats \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox01_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox02_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox03_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox04_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox05_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox06_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox07_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox08_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox09_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox10_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox11_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox12_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox13_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox14_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox15_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox16_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox17_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox18_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox19_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox20_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox21_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox22_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox23_00.fits \
  /mnt/data/1247842824_vis/1247842824_20190722150008_gpubox24_00.fits \
  | tee cotter-1247842824_uvfits_mwa.log
```