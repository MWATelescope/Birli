non-corrected cotter uvfits

```bash
cotter \
  -m /mnt/data/ _vis/1254670392.metafits \
  -o /mnt/data/1254670392_vis/1254670392.cotter.uvfits \
  -allowmissing \
  -edgewidth 0 \
  -endflag 0 \
  -initflag 0 \
  -noantennapruning \
  -nocablelength \
  -noflagautos \
  -noflagdcchannels \
  -nogeom \
  -nosbgains \
  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
  -nostats \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  /mnt/data/1254670392_vis/1254670392_*gpubox*.fits \
  | tee cotter-1254670392-uvfits.log
```

cable corrected cotter uvfits

```bash
cotter \
  -m /mnt/data/1254670392_vis/1254670392.metafits \
  -o /mnt/data/1254670392_vis/1254670392.cotter.cable.uvfits \
  -allowmissing \
  -edgewidth 0 \
  -endflag 0 \
  -initflag 0 \
  -noantennapruning \
  -noflagautos \
  -noflagdcchannels \
  -nogeom \
  -nosbgains \
  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
  -nostats \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  /mnt/data/1254670392_vis/1254670392_*gpubox*.fits \
  | tee cotter-1254670392-uvfits-cable.log
```


geom corrected cotter uvfits

```bash
cotter \
  -m /mnt/data/1254670392_vis/1254670392.metafits \
  -o /mnt/data/1254670392_vis/1254670392.cotter.geom.uvfits \
  -allowmissing \
  -edgewidth 0 \
  -endflag 0 \
  -initflag 0 \
  -noantennapruning \
  -nocablelength \
  -noflagautos \
  -noflagdcchannels \
  -nosbgains \
  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
  -nostats \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  /mnt/data/1254670392_vis/1254670392_*gpubox*.fits \
  | tee cotter-1254670392-uvfits-geom.log
```


both corrected cotter uvfits

```bash
cotter \
  -m /mnt/data/1254670392_vis/1254670392.metafits \
  -o /mnt/data/1254670392_vis/1254670392.cotter.corrected.uvfits \
  -allowmissing \
  -edgewidth 0 \
  -endflag 0 \
  -initflag 0 \
  -noantennapruning \
  -noflagautos \
  -noflagdcchannels \
  -nosbgains \
  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
  -nostats \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  /mnt/data/1254670392_vis/1254670392_*gpubox*.fits \
  | tee cotter-1254670392-uvfits-corrected.log
```

non-corrected birli uvfits

```bash
RUST_LOG=birli=trace birli aoflagger \
    -m /mnt/data/1254670392_vis/1254670392.metafits \
    -u /mnt/data/1254670392_vis/1254670392.birli.none.uvfits \
    --no-geometric-delay \
    --no-cable-delay \
    --emulate-cotter \
    /mnt/data/1254670392_vis/1254670392_*gpubox*.fits \
    | tee birli-1254670392-uvfits-none.log
```

cable corrected birli uvfits

```bash
RUST_LOG=birli=trace birli aoflagger \
    -m /mnt/data/1254670392_vis/1254670392.metafits \
    -u /mnt/data/1254670392_vis/1254670392.birli.cable.uvfits \
    --no-geometric-delay \
    --emulate-cotter \
    /mnt/data/1254670392_vis/1254670392_*gpubox*.fits \
    | tee birli-1254670392-uvfits-cable.log
```

geom corrected birli uvfits

```bash
RUST_LOG=birli=trace birli aoflagger \
    -m /mnt/data/1254670392_vis/1254670392.metafits \
    -u /mnt/data/1254670392_vis/1254670392.birli.geom.uvfits \
    --no-cable-delay \
    --emulate-cotter \
    /mnt/data/1254670392_vis/1254670392_*gpubox*.fits \
    | tee birli-1254670392-uvfits-geom.log
```

both corrected birli uvfits

```bash
RUST_LOG=birli=trace birli aoflagger \
    -m /mnt/data/1254670392_vis/1254670392.metafits \
    -u /mnt/data/1254670392_vis/1254670392.birli.corrected.uvfits \
    --emulate-cotter \
    /mnt/data/1254670392_vis/1254670392_*gpubox*.fits \
    | tee birli-1254670392-uvfits-corrected.log
```

# dump uvfits files

```bash
for i in \
  1254670392_vis/1254670392.cotter.uvfits \
  1254670392_vis/1254670392.cotter.cable.uvfits \
  1254670392_vis/1254670392.cotter.geom.uvfits \
  1254670392_vis/1254670392.cotter.corrected.uvfits \
  1254670392_vis/1254670392.birli.uvfits \
  1254670392_vis/1254670392.birli.cable.uvfits \
  1254670392_vis/1254670392.birli.geom.uvfits \
  1254670392_vis/1254670392.birli.corrected.uvfits 
do 
  python3 tests/data/dump_uvfits.py  --timestep-limit=24 --baseline-limit=24 --dump-mode vis-only "/mnt/data/$i" --dump-csv "/mnt/data/$i.csv" | tee "/mnt/data/$i.txt";
done
```