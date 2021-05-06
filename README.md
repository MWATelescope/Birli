# Birli

<!-- markdownlint-disable MD033 -->
<div class="bg-gray-dark" align="center" style="background-color:#24292e">
<img src="img/birli.png" height="200px" alt="Birli logo">
<br/>
<img src="https://github.com/MWATelescope/Birli/workflows/Cross-Platform%20Tests/badge.svg" alt="Cross-Platform Tests">
<a href="https://codecov.io/gh/MWATelescope/Birli"><img src="https://codecov.io/gh/MWATelescope/Birli/branch/main/graph/badge.svg?token=PK2KYEZOW9" alt="codecov"></a>
</div>

A Rust library for common preprocessing tasks performed in the data pipeline of the Murchison
Widefield Array (MWA), located on the land of the Wajarri Yamatji people in Murchison Shire, Western
Australia.

Birl reads MWA correlator visibilities in the gpufits file format using
[mwalib](https://github.com/MWATelescope/mwalib), which supports the existing "legacy" MWA
correlator, as well as the in-development "MWAX" correlator.

**Birli** is the Wajarri word for lightning, a common cause of outages at the MWA, and a great
descriptor for the speed which this library intends to deliver.

## Installation

### Prerequisites

- A Rust compiler with a version >= 1.51.0 - <https://www.rust-lang.org/tools/install>
- [AOFlagger](https://gitlab.com/aroffringa/aoflagger) >= 3.0
- [CFitsIO](https://heasarc.gsfc.nasa.gov/fitsio/) >= 3.49

for OS-specific instructions, check out the [linux](https://github.com/MWATelescope/Birli/blob/main/.github/workflows/linux_test.yml) and [macOS](https://github.com/MWATelescope/Birli/blob/main/.github/workflows/macos_test.yml) CI Scripts, as these are tested regularly.

### Installing the binary

```bash
cargo install --path .
```

## Usage

`birli -h`

```txt
USAGE:
    birli [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    aoflagger    flag visibilities with aoFlagger
    help         Prints this message or the help of the given subcommand(s)
```

`birli aoflagger -h`

```txt
    birli aoflagger <fits-files>... -f <flag-template> -m <metafits>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f <flag-template>        Sets the template used to name flag files. Percents are substituted for the zero-prefixed
                              GPUBox ID, which can be up to 3 characters log. Similar to -o in Cotter. Example:
                              FlagFile%%%.mwaf
    -m <metafits>             Sets the metafits file.

ARGS:
    <fits-files>...
```

### Examples

A direct comparison between Birli and Cotter might look like this.

```bash
birli aoflagger \
  -m /mnt/data/1247842824_vis/1247842824.metafits \
  -f "/mnt/data/1247842824_vis/flags_birli/Flagfile_Birli_%%.mwaf" \
  /mnt/data/1247842824_vis/1247842824*gpubox*.fits
```

```bash
cotter \
  -m /mnt/data/1247842824_vis/1247842824.metafits \
  -o "/mnt/data/1247842824_vis/flags_cotter/Flagfile_Cotter_%%.mwaf" \
  -nostats \
  -nogeom \
  -noantennapruning \
  -nosbgains \
  -noflagautos \
  -noflagdcchannels \
  -nocablelength \
  -edgewidth 0 \
  -initflag 0 \
  -endflag 0 \
  -flag-strategy  /usr/local/share/aoflagger/strategies/mwa-default.lua \
  /mnt/data/1247842824_vis/1247842824*gpubox*.fits
```

## Acknowledgement

This scientific work uses data obtained from the Murchison Radio-astronomy Observatory. We
acknowledge the Wajarri Yamatji people as the traditional owners of the Observatory site.

## Coverage

<a href="https://codecov.io/gh/MWATelescope/Birli"><img src="https://codecov.io/gh/MWATelescope/Birli/branch/main/graphs/sunburst.svg" height="200px" alt="Birli logo"></a>

This repo is approved by...

<img src="https://github.com/MWATelescope/Birli/raw/main/img/CIRA_Rust_Evangelism_Strike_Force.png" height="200px" alt="CIRA Rust Evangelism Strike Force logo">
