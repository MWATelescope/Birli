# Quality-assurance metrics

`birli --metrics-out <PATH>` writes a FITS file of cheap summaries of the visibilities, which are
useful for finding radio-frequency interference (RFI) and misbehaving tiles without opening the
visibilities themselves. They are computed for each chunk of visibilities after the corrections
and flagging, and before averaging, so they are at the time and frequency resolution of the raw
data.

```bash
birli \
  -m 1254670392.metafits \
  --metrics-out 1254670392_metrics.fits \
  --uvfits-out 1254670392.uvfits \
  1254670392*gpubox*.fits
```

The same metrics can be computed from Rust with the [`birli::metrics`](../src/metrics.rs) module,
for visibilities from any source.

## What is in the file

| metric | HDUs | what it is |
| --- | --- | --- |
| autocorrelations | `AUTO_*` | the bandpass of each tile, a polynomial fit to it, its delay transform, and the autocorrelations with the bandpass subtracted |
| SSINS | `SSINS_*` | z-scores of the baseline-averaged amplitude of the visibilities differenced in time |
| EAVILS | `EAVILS_*` | z-scores of the baseline-averaged amplitude of the visibilities |
| flag occupancy | `AO_FLAG_METRICS` | fraction of baselines flagged at each time and frequency (`aoflagger` feature only) |
| short baseline delays | `CROSS_DELAY_*` | delay transform of the cross-correlations on short baselines, for each tile (`aoflagger` feature only) |

Things to know about all of them:

- Timesteps that Birli flags (for example the quack time, `--flag-times`), and channels that Birli
  flags (`--flag-coarse-chans`, `--flag-fine-chans`, `--flag-raw-chans`, `--flag-edge-width`, the DC
  channel) are left out of every metric, and are `NaN` in the images.
- Only SSINS and the flag occupancy look at the flags of individual visibilities, like those from
  AOFlagger.
- Every image is written as 64-bit floats.
- Image axes are described with `CTYPEn`, `CRVALn`, `CRPIXn`, and where it is constant, `CDELTn`.
  Times are GPS seconds, frequencies are in Hz and delays are in ns. FITS axis 1 is the fastest
  axis, which is the last axis of the array in Python.
- Tiles are in order of antenna index (their order in the metafits), and only selected tiles
  (`--sel-ants`, `--no-sel-flagged-ants`) are included.

### Which file

- When the observation is processed in more than one chunk of timesteps (`--time-chunk` or
  `--max-memory`), each chunk has its own metrics, in its own file, with the timestep indices of
  the chunk before the extension: `metrics_ts0-29.fits`, `metrics_ts30-59.fits`.
- When the coarse channels are not contiguous, or are selected with `--sel-chan-ranges`, each range
  of channels has its own file, named like the visibility files: `metrics_ch62-63.fits`.
- No file is written for a chunk where every timestep or every channel is flagged.

### Autocorrelations

These need autocorrelations, so they are empty with `--no-sel-autos`. The amplitude of the
autocorrelations is used, for all four polarisations.

| HDU | type | FITS axes | contents |
| --- | --- | --- | --- |
| `AUTO_SUB_ANT=<tile name>`, one for each tile | image | frequency, time, Stokes | the autocorrelation amplitude, with the bandpass of the tile subtracted |
| `AUTO_POL=<pol>` for `XX`, `YY`, `XY`, `YX` | image | frequency, tile | the bandpass: the autocorrelation amplitude averaged over unflagged timesteps |
| `AUTO_COEFFS` | table | one row for each tile | the metadata of the tile, and the coefficients of a polynomial fit to the bandpass |
| `AUTO_DELAY_POL=<pol>` for `XX`, `YY`, `XY`, `YX` | image | delay, tile | the delay transform of the bandpass |

- The Stokes axis of `AUTO_SUB_ANT=*` is in the FITS order: `XX`, `YY`, `XY`, `YX` (`CRVAL3 = -5`,
  `CDELT3 = -1`).
- The frequency axis has no `CDELT1`, because CARTA does not display the images well with it. The
  channel width is in `CHAN_WIDTH`.
- The header of each `AUTO_SUB_ANT=*` HDU, and each row of `AUTO_COEFFS` has the metadata of the
  tile:

  | key or column | contents |
  | --- | --- |
  | `ANTNAME` (`ANT_NAME` in the table) | tile name |
  | `ANT_ID` | antenna index: the order of the tile in the metafits |
  | `ANT_NUM` | tile ID |
  | `RX_NUMBER`, `RX_SLOT`, `RX_TYPE` | the receiver of the tile, its slot in the receiver, and the type of the receiver |
  | `ANT_TYPE` | the same as `RX_TYPE` |
  | `CABLE_FLAVOUR` | the flavour of the cable, like `RG6_150` |
  | `WHITENING_FILTER` | 1 if the tile has a whitening filter |
  | `OBSGEO-X`, `OBSGEO-Y`, `OBSGEO-Z` | the position of the tile relative to the centre of the array in metres: north, east and height. These are *not* geocentric coordinates. |

  The receiver and cable are those of the X polarisation.
- `AUTO_COEFFS` has the columns `COEFFS_XX`, `COEFFS_YY`, `COEFFS_XY` and `COEFFS_YX`, each with
  `POLY_ORD + 1` coefficients `c`, lowest order first, so that the fit to the bandpass at frequency
  `f` in Hz is `c[0] + c[1] f + c[2] f²`. The coefficients are `NaN` if there are too few unflagged
  channels to fit.
- The delay transform is the amplitude of the Fourier transform of the bandpass along frequency,
  with a Blackman window, zero-padded to a resolution of about 1 ns, between 100 and 3000 ns. A
  reflection in a cable shows up as a peak at the delay of the reflection. Flagged channels are
  treated as zero. The delay of the first bin is `CRVAL1`, and the width of the bins is
  `DELAY_RES`, both in ns. With channels wider than about 170 kHz, delays above half the reciprocal
  of the channel width are aliases.

### SSINS

The amplitude of the difference between the visibilities at adjacent timesteps, averaged over
cross-correlation baselines, as a z-score. The sky changes slowly so it cancels in the difference,
leaving noise and anything that changes quickly, like RFI. This is the incoherent noise spectrum of
[SSINS](https://ssins.readthedocs.io/en/latest/) ([Wilensky et al.
2019](https://arxiv.org/abs/1906.01093)), and matches
`SSINS.INS(SS.read(..., diff=True), spectrum_type="cross", order=0).metric_ms`:

- only cross-correlations are averaged,
- a baseline is left out of a cell if either of the visibilities that are differenced is flagged,
- the mean over time is weighted by the number of baselines in each cell.

With `N` unflagged baselines in a cell, and `C = 4/π - 1`:

```txt
zscore(t, f) = sqrt(N(t, f) / C) * (mean_amp(t, f) / mean_amp(f) - 1)
```

There is one fewer timestep than in the visibilities, and the time axis starts halfway between the
first two timesteps. The metrics are left out if a chunk has only one timestep.

| HDU | FITS axes | contents |
| --- | --- | --- |
| `SSINS_POL=<pol>` for `XX`, `XY`, `YX`, `YY` | frequency, time | the z-scores |
| `SSINS_FLAGS` | frequency, time | 1 where the z-scores are flagged. With the `aoflagger` feature, these are from AOFlagger's generic minimal strategy; otherwise they are all 0 |
| `SSINS_N_BL` | frequency, time | `N`, the number of unflagged cross-correlation baselines in each cell |
| `SSINS_DIFF_MEAN_AMP_FP` | polarisation, frequency | `mean_amp(f)`, the weighted mean over time of the baseline-averaged amplitude of the difference |

`N_BL` in the header of `SSINS_POL=*` is the number of cross-correlation baselines selected.

### EAVILS

The amplitude of the visibilities, averaged over all selected baselines including autocorrelations,
as a z-score relative to its mean over time. This is like SSINS without differencing in time, so it
also sees things that don't change between timesteps. With the amplitudes `amp` as an array with
dimensions (time, baseline, frequency, polarisation):

```python
blmean = np.mean(amp, axis=1)
stdv = np.sqrt(np.mean(np.var(amp, axis=0), axis=0))
zscore = (blmean - np.mean(blmean, axis=0)) / stdv * np.sqrt(Nbls)
```

where the means and variances over time only include unflagged timesteps.

| HDU | FITS axes | contents |
| --- | --- | --- |
| `EAVILS_POL=<pol>` for `XX`, `XY`, `YX`, `YY` | frequency, time | the z-scores |
| `EAVILS_FLAGS` | frequency, time | all 0 |
| `EAVILS_MEAN_AMP_FP` | polarisation, frequency | `np.mean(blmean, axis=0)` |
| `EAVILS_SQRT_MEAN_VAR_AMP_FP` | polarisation, frequency | `stdv` |

### Flag occupancy

`AO_FLAG_METRICS` is an image with FITS axes (frequency, time) of the fraction of the selected
cross-correlation baselines which are flagged, for example by AOFlagger. It is only written by a
Birli with the `aoflagger` feature.

### Short baseline delays

`CROSS_DELAY_POL=<pol>` for `XX`, `YY`, `XY`, `YX` are images with FITS axes (delay, tile). The
visibilities of each cross-correlation baseline no longer than `--cross-metrics-cutoff` metres
(default 30) are averaged coherently over unflagged timesteps, and the amplitude of that spectrum
is delay transformed like the autocorrelations. The metric for a tile is the mean of the delay
spectra of the short baselines it is in, or `NaN` if it has none. They are only written by a Birli
with the `aoflagger` feature.

## Reading the file

```python
from astropy.io import fits
import numpy as np

with fits.open("1254670392_metrics.fits") as hdus:
    hdus.info()

    # SSINS z-scores: (time, frequency)
    hdu = hdus["SSINS_POL=XX"]
    zscore = hdu.data
    freqs_hz = hdu.header["CRVAL1"] + hdu.header["CDELT1"] * np.arange(zscore.shape[1])
    gps_times = hdu.header["CRVAL2"] + hdu.header["CDELT2"] * np.arange(zscore.shape[0])

    # bandpass of each tile: (tile, frequency), in the order of the rows of AUTO_COEFFS
    bandpass = hdus["AUTO_POL=XX"].data
    tiles = hdus["AUTO_COEFFS"].data
    for name, flavour, spectrum in zip(tiles["ANT_NAME"], tiles["CABLE_FLAVOUR"], bandpass):
        print(name, flavour, np.nanmedian(spectrum))

    # delay spectrum of each tile: (tile, delay)
    hdu = hdus["AUTO_DELAY_POL=XX"]
    delay_spectrum = hdu.data
    delays_ns = hdu.header["CRVAL1"] + hdu.header["DELAY_RES"] * np.arange(delay_spectrum.shape[1])
```
