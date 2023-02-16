# Test Data

Full MWA files can get quite large, so this script generates truncated versions of these files for testing. You will not need to run the scripts in this directory to test Birli, the purpose here is to capture how the files were generated.

The test data files in this directory are generated with `generate_test_data.py`. This uses the [`astropy`](https://pypi.org/project/astropy/) Python module to generate mock visibility data.

## 1297526432 - MWAX HydA files

This is a synthetic observation based off the MWAX [HydA Gold Observation](https://wiki.mwatelescope.org/display/MP/HydA+Gold+Observation), limited to two coarse channels, two batches, two scans per batch, two antennas and two fine channels per coarse to allow for the smallest possible set of files which could be representative of an mwax observation for testing purposes:

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

This is a synthetic observation based off [obs id 1196175296](http://ws.mwatelescope.org/observation/obs/?obs_id=1196175296), limited to two coarse channels, two batches, two scans per batch and two fine channels per coarse to allow for the smallest possible set of files which could be representative of an observation for testing purposes.

The observation values are filled with unique positive 32-bit float values, increasing monotonically in order of their position in the file. When reading the visibilities, the order and sign of these values is dificult to predict because of the quirks of legacy Ord correlator gpubox file format. so a CSV dump of the visibilities is provided in hex and decimal.

Because this observation only contains a subset of the channels of a real observation, an additional "cotter friendly" metafits file is provided, with `FINECHAN`, `BANDWDTH` and `NCHANS` set to values that will prevent cotter from segfaulting when reading the files, however Cotter does not correctly calculate the sky frequencies for this partial observation.

Additional "adjusted" gpufits files were generated from this set in order to test edge cases in start time offsets. See [cotter_timesteps.md](../../doc/cotter_timesteps.md)

## 1247842824 - Flags

This is a subset of [obs id 1247842824](http://ws.mwatelescope.org/observation/obs/?obs_id=1247842824), limited to one coarse channel, 1 batch, two scans per batch and two fine channels per coarse to allow for the smallest possible set of files which could be representative of an observation for testing RFI flagging.

Because this observation only contains a subset of the channels of a real observation, an additional "cotter friendly" metafits file is provided, with `FINECHAN`, `BANDWDTH` and `NCHANS` set to values that will prevent cotter from segfaulting when reading the files, however Cotter does not correctly calculate the sky frequencies for this partial observation.

## 1254670392 - Averaged

This is the full [obs id 1254670392](http://ws.mwatelescope.org/observation/obs/?obs_id=1254670392), which is averaged to 2s / 40kHz. Because it is nice and small, all channels fit inside the repository, and a partial CSV dump of the uvfits visibilities is used to test UVFits and corrections, however the uvfits files themselves are too large.

Because the file cuts off 2 scans in, and has all antennae flagged, a fixed version of the metafits file was made for better comparison with cotter.

```python
>>> import os
>>> os.getcwd()
'/home/dev/Birli/tests/data/1254670392_avg'
>>> from astropy.io import fits
>>> meta_fits = fits.open("1254670392.metafits")
>>> meta_fits[0].header["NSCANS"]
4
>>> meta_fits[0].header["NSCANS"] = 2
>>> for row in meta_fits[1].data:
...   row['Flag'] = 0
>>> meta_fits.writeto("1254670392.fixed.metafits", overwrite=True)
```

## 1119683928 - Picket

this is the smallest calibrator observation with 24 channels which has non-contiguous coarse channels. It was found
using the following ADQL query:

```sql
SELECT TOP 100 obs_id, obsname, total_archived_data_bytes, gpubox_files_archived, calibrators, channels_are_contiguous
FROM mwa.observation
WHERE total_archived_data_bytes > 0
  AND gpubox_files_archived > 22
ORDER BY total_archived_data_bytes ASC
```

## Generating test files

warning: the following line needs to be patched out in Cotter for this to work. https://github.com/MWATelescope/cotter/blob/b60b633b6b4b11bc5c84162bd0058a24d98f650c/gpufilereader.cpp#L195

```bash
python3 tests/data/generate_test_data.py | tee generate.log

python3 tests/data/adjust_gpufits.py \
  --in-file=tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits \
  --out-file=tests/data/1196175296_mwa_ord/adjusted_-1/1196175296_20171201145440_gpubox01_00.fits \
  --corr-type=MWA_ORD --timestep-offset=-1
# python3 tests/data/adjust_gpufits.py \
#   --in-file=tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits \
#   --out-file=tests/data/1196175296_mwa_ord/adjusted_-1/1196175296_20171201145540_gpubox01_01.fits \
#   --corr-type=MWA_ORD --timestep-offset=-1
# python3 tests/data/adjust_gpufits.py \
#   --in-file=tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits \
#   --out-file=tests/data/1196175296_mwa_ord/adjusted_+1/1196175296_20171201145440_gpubox01_00.fits \
#   --corr-type=MWA_ORD --timestep-offset=1
# python3 tests/data/adjust_gpufits.py \
#   --in-file=tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits \
#   --out-file=tests/data/1196175296_mwa_ord/adjusted_+1/1196175296_20171201145540_gpubox01_01.fits \
#   --corr-type=MWA_ORD --timestep-offset=1
# python3 tests/data/adjust_gpufits.py \
#   --in-file=tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits \
#   --out-file=tests/data/1196175296_mwa_ord/limited_1/1196175296_20171201145440_gpubox01_00.fits \
#   --corr-type=MWA_ORD --timestep-limit=1
python3 tests/data/adjust_gpufits.py \
  --in-file=tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits \
  --out-file=tests/data/1196175296_mwa_ord/limited_1/1196175296_20171201145540_gpubox01_01.fits \
  --corr-type=MWA_ORD --timestep-limit=1
# python3 tests/data/adjust_gpufits.py \
#   --in-file=tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits \
#   --out-file=tests/data/1196175296_mwa_ord/limited_0/1196175296_20171201145440_gpubox01_00.fits \
#   --corr-type=MWA_ORD --timestep-limit=0
# python3 tests/data/adjust_gpufits.py \
#   --in-file=tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits \
#   --out-file=tests/data/1196175296_mwa_ord/limited_0/1196175296_20171201145540_gpubox01_01.fits \
#   --corr-type=MWA_ORD --timestep-limit=0
# python3 tests/data/adjust_gpufits.py \
#   --in-file=tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits \
#   --out-file=tests/data/1196175296_mwa_ord/empty/1196175296_20171201145440_gpubox01_00.fits \
#   --corr-type=MWA_ORD --timestep-limit=1 --empty-data
# python3 tests/data/adjust_gpufits.py \
#   --in-file=tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits \
#   --out-file=tests/data/1196175296_mwa_ord/empty/1196175296_20171201145540_gpubox01_01.fits \
#   --corr-type=MWA_ORD --timestep-limit=1 --empty-data


# Cotter flags on 1196175296_mwa_ord with generic flagging
cotter \
  -m tests/data/1196175296_mwa_ord/1196175296.cotter.metafits \
  -o "tests/data/1196175296_mwa_ord/FlagfileCotterGeneric%%.mwaf" \
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
  -sbpassband tests/data/subband-passband-2ch-unitary.txt \
  -nostats \
  -sbcount 2 \
  -sbstart 1 \
  -flag-strategy /usr/local/share/aoflagger/strategies/generic-minimal.lua \
  tests/data/1196175296_mwa_ord/1196175296_*gpubox*.fits \
  | tee cotter-1196175296-generic.log
# Cotter flags on 1196175296_mwa_ord with MWA flagging
cotter \
  -m tests/data/1196175296_mwa_ord/1196175296.cotter.metafits \
  -o "tests/data/1196175296_mwa_ord/FlagfileCotterMWA%%.mwaf" \
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
  -sbpassband tests/data/subband-passband-2ch-unitary.txt \
  -nostats \
  -sbcount 2 \
  -sbstart 1 \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1196175296_mwa_ord/1196175296_*gpubox*.fits \
  | tee cotter-1196175296-mwa.log
# Cotter flags on 1196175296_mwa_ord with MWA flagging, 01_00 offset by -1
cotter \
  -m tests/data/1196175296_mwa_ord/1196175296.cotter.metafits \
  -o "tests/data/1196175296_mwa_ord/adjusted_-1/FlagfileCotterMWA%%.mwaf" \
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
  -sbpassband tests/data/subband-passband-2ch-unitary.txt \
  -nostats \
  -sbcount 2 \
  -sbstart 1 \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1196175296_mwa_ord/adjusted_-1/1196175296_20171201145440_gpubox01_00.fits \
  tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits \
  tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits \
  tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits \
  | tee cotter-1196175296-mwa_-1.log
# Cotter flags on 1196175296_mwa_ord with MWA flagging, 01_00 and 01_01 offset by -1
cotter \
  -m tests/data/1196175296_mwa_ord/1196175296.cotter.metafits \
  -o "tests/data/1196175296_mwa_ord/adjusted_-1/FlagfileCotterMWA_both%%.mwaf" \
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
  -sbpassband tests/data/subband-passband-2ch-unitary.txt \
  -nostats \
  -sbcount 2 \
  -sbstart 1 \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1196175296_mwa_ord/adjusted_-1/1196175296_20171201145440_gpubox01_00.fits \
  tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits \
  tests/data/1196175296_mwa_ord/adjusted_-1/1196175296_20171201145540_gpubox01_01.fits \
  tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits \
  | tee cotter-1196175296-mwa_-1_both.log
# Cotter flags on 1196175296_mwa_ord with MWA flagging, 01_00 offset by +1
cotter \
  -m tests/data/1196175296_mwa_ord/1196175296.cotter.metafits \
  -o "tests/data/1196175296_mwa_ord/adjusted_+1/FlagfileCotterMWA%%.mwaf" \
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
  -sbpassband tests/data/subband-passband-2ch-unitary.txt \
  -nostats \
  -sbcount 2 \
  -sbstart 1 \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits \
  tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits \
  tests/data/1196175296_mwa_ord/adjusted_+1/1196175296_20171201145540_gpubox01_01.fits \
  tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits \
  | tee cotter-1196175296-mwa_+1.log
# Cotter flags on 1196175296_mwa_ord with MWA flagging, 01_00 and 01_01 offset by +1
cotter \
  -m tests/data/1196175296_mwa_ord/1196175296.cotter.metafits \
  -o "tests/data/1196175296_mwa_ord/adjusted_+1/FlagfileCotterMWA_both%%.mwaf" \
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
  -sbpassband tests/data/subband-passband-2ch-unitary.txt \
  -nostats \
  -sbcount 2 \
  -sbstart 1 \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1196175296_mwa_ord/adjusted_+1/1196175296_20171201145440_gpubox01_00.fits \
  tests/data/1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits \
  tests/data/1196175296_mwa_ord/adjusted_+1/1196175296_20171201145540_gpubox01_01.fits \
  tests/data/1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits \
  | tee cotter-1196175296-mwa_+1_both.log

# Cotter uvfits on 1196175296_mwa_ord with MWA flagging
cotter \
  -m tests/data/1196175296_mwa_ord/1196175296.cotter.metafits \
  -o tests/data/1196175296_mwa_ord/1196175296.uvfits \
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
  -sbpassband tests/data/subband-passband-2ch-unitary.txt \
  -sbstart 1 \
  -sbcount 2 \
  -nostats \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1196175296_mwa_ord/1196175296_*gpubox*.fits \
  | tee cotter-1196175296-uvfits.log
# # Cotter uvfits on 1196175296_mwa_ord with MWA flagging 4s averaging
# cotter \
#   -m tests/data/1196175296_mwa_ord/1196175296.cotter.metafits \
#   -o tests/data/1196175296_mwa_ord/1196175296_4s.uvfits \
#   -allowmissing \
#   -edgewidth 0 \
#   -endflag 0 \
#   -initflag 0 \
#   -noantennapruning \
#   -nocablelength \
#   -noflagautos \
#   -noflagdcchannels \
#   -nogeom \
#   -nosbgains \
#   -sbpassband tests/data/subband-passband-2ch-unitary.txt \
#   -nostats \
#   -timeres 4 \
#   -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
#   tests/data/1196175296_mwa_ord/1196175296_*gpubox*.fits \
#   | tee cotter-1196175296-uvfits.log

# Cotter flags on 1247842824_flags with generic flagging
cotter \
  -m tests/data/1247842824_flags/1247842824.cotter.metafits \
  -o "tests/data/1247842824_flags/FlagfileCotterGeneric%%.mwaf" \
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
  -sbpassband tests/data/subband-passband-128ch-unitary.txt \
  -nostats \
  -sbcount 1 \
  -sbstart 1 \
  -flag-strategy /usr/local/share/aoflagger/strategies/generic-minimal.lua \
  tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits \
  | tee cotter-1247842824-generic.log
# Cotter flags on 1247842824_flags with MWA flagging
cotter \
  -m tests/data/1247842824_flags/1247842824.cotter.metafits \
  -o "tests/data/1247842824_flags/FlagfileCotterMWA%%.mwaf" \
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
  -sbpassband tests/data/subband-passband-128ch-unitary.txt \
  -nostats \
  -sbcount 1 \
  -sbstart 1 \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits \
  | tee cotter-1247842824-mwa.log
# Cotter uvfits on 1247842824_flags with MWA flagging
cotter \
  -m tests/data/1247842824_flags/1247842824.cotter.metafits \
  -o tests/data/1247842824_flags/1247842824.uvfits \
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
  -sbpassband tests/data/subband-passband-128ch-unitary.txt \
  -nostats \
  -sbcount 1 \
  -sbstart 1 \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits \
  | tee cotter-1247842824-uvfits.log
# Cotter uvfits on 1254670392_avg with MWA flagging, no corrections
cotter \
  -m tests/data/1254670392_avg/1254670392.fixed.metafits \
  -o tests/data/1254670392_avg/1254670392.cotter.none.uvfits \
  -allowmissing \
  -edgewidth 0 \
  -endflag 0 \
  -initflag 0 \
  -noantennapruning \
  -nocablelength \
  -nogeom \
  -noflagautos \
  -noflagdcchannels \
  -nosbgains \
  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
  -nostats \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1254670392_avg/1254670392*gpubox*.fits \
  | tee cotter-1254670392-uvfits-none.log
# Cotter uvfits on 1254670392_avg with MWA flagging, cable corrections only
cotter \
  -m tests/data/1254670392_avg/1254670392.fixed.metafits \
  -o tests/data/1254670392_avg/1254670392.cotter.cable.uvfits \
  -allowmissing \
  -edgewidth 0 \
  -endflag 0 \
  -initflag 0 \
  -noantennapruning \
  -nogeom \
  -noflagautos \
  -noflagdcchannels \
  -nosbgains \
  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
  -nostats \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1254670392_avg/1254670392*gpubox*.fits \
  | tee cotter-1254670392-uvfits-cable.log
# Cotter uvfits on 1254670392_avg with MWA flagging, geometric corrections only
cotter \
  -m tests/data/1254670392_avg/1254670392.fixed.metafits \
  -o tests/data/1254670392_avg/1254670392.cotter.geom.uvfits \
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
  tests/data/1254670392_avg/1254670392*gpubox*.fits \
  | tee cotter-1254670392-uvfits-geom.log
# Cotter uvfits on 1254670392_avg with MWA flagging, both geometric and cable corrections
cotter \
  -m tests/data/1254670392_avg/1254670392.fixed.metafits \
  -o tests/data/1254670392_avg/1254670392.cotter.corrected.uvfits \
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
  -flag-strategy /usr/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1254670392_avg/1254670392*gpubox*.fits \
  | tee cotter-1254670392-uvfits-corrected.log
# Cotter uvfits on 1254670392_avg with MWA flagging, no corrections and 4s/160kHz averaging
cotter \
  -m tests/data/1254670392_avg/1254670392.fixed.metafits \
  -o tests/data/1254670392_avg/1254670392.cotter.none.avg_4s_160khz.uvfits \
  -allowmissing \
  -edgewidth 0 \
  -endflag 0 \
  -initflag 0 \
  -noantennapruning \
  -nocablelength \
  -nogeom \
  -noflagautos \
  -noflagdcchannels \
  -nosbgains \
  -sbpassband tests/data/subband-passband-32ch-unitary.txt \
  -nostats \
  -flag-strategy /usr/share/aoflagger/strategies/mwa-default.lua \
  -timeres 4 \
  -freqres 160 \
  tests/data/1254670392_avg/1254670392*gpubox*.fits \
  | tee cotter-1254670392-uvfits-none-avg_4s_160khz.log

# dump metafits files
# python -m pip install -r tests/data/requirements.txt
for i in \
  1196175296_mwa_ord/1196175296.metafits \
  1247842824_flags/1247842824.metafits \
  1297526432_mwax/1297526432.metafits \
  1254670392_avg/1254670392.metafits \
  1254670392_avg/1254670392.fixed.metafits
do
  python3 tests/data/dump_metafits.py "tests/data/$i" | tee "tests/data/$i.txt"
done
# dump gpufits files
for line in \
  "1196175296_mwa_ord/1196175296_20171201145440_gpubox01_00.fits|MWA_ORD" \
  "1196175296_mwa_ord/1196175296_20171201145440_gpubox02_00.fits|MWA_ORD" \
  "1196175296_mwa_ord/1196175296_20171201145540_gpubox01_01.fits|MWA_ORD" \
  "1196175296_mwa_ord/1196175296_20171201145540_gpubox02_01.fits|MWA_ORD" \
  "1196175296_mwa_ord/adjusted_-1/1196175296_20171201145440_gpubox01_00.fits|MWA_ORD" \
  "1196175296_mwa_ord/limited_1/1196175296_20171201145540_gpubox01_01.fits|MWA_ORD" \
  "1297526432_mwax/1297526432_20210216160014_ch117_000.fits|MWAX" \
  "1297526432_mwax/1297526432_20210216160014_ch117_001.fits|MWAX" \
  "1297526432_mwax/1297526432_20210216160014_ch118_000.fits|MWAX" \
  "1297526432_mwax/1297526432_20210216160014_ch118_001.fits|MWAX" \
  "1247842824_flags/1247842824_20190722150008_gpubox01_00.fits|MWA_ORD"
do
  python3 tests/data/dump_gpufits.py tests/data/${line%|*} --corr-type=${line#*|} | tee tests/data/${line%|*}.txt
done
# dump mwaf files
for i in \
  1196175296_mwa_ord/FlagfileCotterGeneric01.mwaf \
  1196175296_mwa_ord/FlagfileCotterGeneric02.mwaf \
  1196175296_mwa_ord/FlagfileCotterMWA01.mwaf \
  1196175296_mwa_ord/FlagfileCotterMWA02.mwaf \
  1247842824_flags/FlagfileCotterMWA01.mwaf \
  1247842824_flags/FlagfileCotterGeneric01.mwaf
do
  python3 tests/data/dump_mwaf.py --timestep-limit=12 --baseline-limit=12 "tests/data/$i" | tee "tests/data/$i.txt"
done
# dump uvfits files
for i in \
  1196175296_mwa_ord/1196175296.uvfits \
  1247842824_flags/1247842824.uvfits \
  1254670392_avg/1254670392.cotter.none.uvfits \
  1254670392_avg/1254670392.cotter.none.uvfits \
  1254670392_avg/1254670392.cotter.none.avg_4s_160khz.uvfits \
  1254670392_avg/1254670392.cotter.cable.uvfits \
  1254670392_avg/1254670392.cotter.geom.uvfits \
  1254670392_avg/1254670392.cotter.corrected.uvfits
do
  python3 tests/data/dump_uvfits.py "tests/data/$i" \
    --timestep-limit=12 --baseline-limit=12 \
    --dump-csv "tests/data/$i.csv" --dump-mode vis-weight \
    | tee "tests/data/$i.txt";
done

python3 tests/data/dump_uvfits.py \
  tests/data/1196175296_mwa_ord/1196175296.uvfits \
  --timestep-limit=12 --baseline-limit=12 \
  --dump-mode=vis-only \
  --dump-csv=tests/data/1196175296_mwa_ord/1196175296.uvfits.vis.csv \
  | tee tests/data/1196175296_mwa_ord/1196175296.uvfits.txt
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