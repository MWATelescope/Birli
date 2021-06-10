# Test Data

Full MWA files can get quite large, so this script generates truncated versions of these files for testing. You will not need to run the scripts in this directory to test Birli, the purpose here is to capture how the files were generated.

The test data files in this directory are generated with `generate_test_data.py`. This uses the [`astropy`](https://pypi.org/project/astropy/) Python module to generate mock visibility data.

## 1297526432 - MWAX HydA files

This set of fits is a modified version of the MWAX [HydA Gold Observation](https://wiki.mwatelescope.org/display/MP/HydA+Gold+Observation), limited to two coarse channels, two batches, two scans, two antennas and two fine channels to allow for the smallest possible set of files which could be representative of an observation for testing purposes:

- Coarse channels:
  - `(117) 149.760 Mhz`
  - `(118) 151.040 Mhz`
- Timesteps (first two timesteps from first two batches):
  - `unix=1613491214.000` (`batch=000, scan=0`)
  - `unix=1613491214.500` (`batch=000, scan=1`)
  - `unix=1613491294.000` (`batch=001, scan=0`)
  - `unix=1613491294.500` (`batch=001, scan=1`)
- Baselines:
  - `Tile051` vs `Tile051`
  - `Tile051` vs `Tile052`
  - `Tile052` vs `Tile052`

The observation values are filled with unique 32-bit float values, increasing monotonically in order of `coarse_channel`, `batch`, `timestep`, `baseline`, `fine_channel`, `correlator_polarization`, and `complex_component` (real or imaginary). The floats are made of 3 byte-aligned components to be able to track the exact lineage of each float value:

- the bytes `'A' = 0x41`
- a global ImageHDU index, unique to each ImageHDU (`coarse_channel`, `batch`, `timestep`)
- a float offset within the ImageHDU, where each ImageHDU has the dimensions:
  - axis 2: `num_baselines = 3`
  - axis 1: `num_fine_channels (2) * num_correlator_pols (4) * floats_per_complex (2) = 16`

example:

```txt
0000 0000 0100 0001 0000 0111 0010 1010 -> 0x41032a = 4260650
          |       |       |||   || ||||
          |       |       |||   || |||| ImageHDU offset(0x2a):
          |       |       |||   || |||. -> complex_component = 1 (Re)
          |       |       |||   || |..  -> correlator_polarization = 1 (XY)
          |       |       |||   || .    -> fine_channel = 1
          |       |       |||   ..      -> baseline = 2 (Tile052 vs Tile052)
          |       |       |||
          |       |       |||           global Scan index (0x03):
          |       |       ||.           ->
          |       |       |..           -> scan=1, batch=1 (timestep:1613491294.500)
          |       |       .             -> coarse channel = 1 (118) 151.040 Mhz
          |       |
          |       |                     'A' Byte marker:
          |<----->|                     -> 0x41 = 'A'
```

## 1196175296 - MWA Ord

This is an observation from the legacy MWA Ord Correlator in 2017.

## Generating test files

```bash
python3 tests/data/generate_test_data.py | tee generate.log
cotter \
  -m tests/data/1196175296_mwa_ord/1196175296_cotter-friendly.metafits \
  -o tests/data/1196175296_mwa_ord/FlagfileCotter%%.mwaf \
  -allowmissing \
  -edgewidth 0 \
  -endflag 0 \
  -initflag 0 \
  -noantennapruning \
  -nocablelength \
  -noflagautos \
  -noflagdcchannels \
  -nogeom \
  -norfi \
  -nosbgains \
  -nostats \
  -sbcount 2 \
  -sbstart 1 \
  -flag-strategy /usr/local/share/aoflagger/strategies/generic-minimal.lua \
  tests/data/1196175296_mwa_ord/1196175296_*gpubox*_00.fits \
  | tee cotter-1196175296-flag.log
cotter \
  -m tests/data/1196175296_mwa_ord/1196175296_cotter-friendly.metafits \
  -o tests/data/1196175296_mwa_ord/1196175296.uvfits \
  -allowmissing \
  -noantennapruning \
  -nocablelength \
  -noflagautos \
  -noflagdcchannels \
  -nogeom \
  -norfi \
  -nosbgains \
  -nostats \
  -sbcount 2 \
  -sbstart 1 \
  tests/data/1196175296_mwa_ord/1196175296_*gpubox*_00.fits \
  | tee cotter-1196175296-uvfits.log

for i in \
  1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits \
  1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits \
  1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits \
  1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits \
  1196175296_mwa_ord/1196175296.metafits \
  1196175296_mwa_ord/FlagfileCotter01.mwaf \
  1196175296_mwa_ord/FlagfileCotter02.mwaf \
  1196175296_mwa_ord/1196175296.uvfits \
  1297526432_mwax/1297526432_20210216160014_ch117_000.fits \
  1297526432_mwax/1297526432_20210216160014_ch117_001.fits \
  1297526432_mwax/1297526432_20210216160014_ch118_000.fits \
  1297526432_mwax/1297526432_20210216160014_ch118_001.fits \
  1297526432_mwax/1297526432.metafits
do 
  fitsheader "tests/data/$i" | tee "tests/data/$i.header.txt";
done
```

then with [mwa-scratchpad](https://github.com/derwentx/mwa-scratchpad)

```bash
cargo run dump-all-data \
  --dump-filename=../Birli/tests/data/1297526432_mwax/1297526432_dump.csv \
  --metafits=../Birli/tests/data/1297526432_mwax/1297526432.metafits \
  ../Birli/tests/data/1297526432_mwax/1297526432_20210216160014_ch*.fits \
  | tee dump-1297526432.log
cargo run dump-all-data --vis-radix=16 --absolute \
  --dump-filename=../Birli/tests/data/1297526432_mwax/1297526432_dump_hex.csv \
  --metafits=../Birli/tests/data/1297526432_mwax/1297526432.metafits \
  ../Birli/tests/data/1297526432_mwax/1297526432_20210216160014_ch*.fits \
  | tee dump-1297526432-hex.log
cargo run dump-all-data --vis-radix=2 --absolute \
  --dump-filename=../Birli/tests/data/1297526432_mwax/1297526432_dump_bin.csv \
  --metafits=../Birli/tests/data/1297526432_mwax/1297526432.metafits \
  ../Birli/tests/data/1297526432_mwax/1297526432_20210216160014_ch*.fits \
  | tee dump-1297526432-bin.log

cargo run dump-all-data \
  --dump-filename=../Birli/tests/data/1196175296_mwa_ord/1196175296_dump.csv \
  --metafits=../Birli/tests/data/1196175296_mwa_ord/1196175296.metafits \
  ../Birli/tests/data/1196175296_mwa_ord/1196175296_*.fits \
  | tee dump-1196175296.log
cargo run dump-all-data --vis-radix=16 --absolute \
  --dump-filename=../Birli/tests/data/1196175296_mwa_ord/1196175296_dump_hex.csv \
  --metafits=../Birli/tests/data/1196175296_mwa_ord/1196175296.metafits \
  ../Birli/tests/data/1196175296_mwa_ord/1196175296_*.fits \
  | tee dump-1196175296-hex.log
```