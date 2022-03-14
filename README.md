# Birli

<!-- markdownlint-disable MD033 -->
<div class="bg-gray-dark" align="center" style="background-color:#24292e">
<img src="img/birli.png" height="200px" alt="Birli logo">
<br/>
<a href="https://github.com/MWATelescope/Birli/actions/workflows/linux_test.yml">
  <img src="https://github.com/MWATelescope/Birli/actions/workflows/linux_test.yml/badge.svg" alt="MacOS Tests"></a>
<a href="https://github.com/MWATelescope/Birli/actions/workflows/macos_test.yml">
  <img src="https://github.com/MWATelescope/Birli/actions/workflows/macos_test.yml/badge.svg" alt="Linix Tests"></a>
<a href="https://crates.io/crates/birli">
  <img alt="Crates.io" src="https://img.shields.io/crates/d/birli?label=crates.io%20%E2%AC%87%EF%B8%8F"></a>
<a href="https://docs.rs/crate/birli/">
  <img src="https://docs.rs/birli/badge.svg" alt="codecov"></a>
<a href="https://codecov.io/gh/MWATelescope/Birli">
  <img src="https://codecov.io/gh/MWATelescope/Birli/branch/main/graph/badge.svg?token=PK2KYEZOW9" alt="codecov"></a>
<a href="https://rust-reportcard.xuri.me/report/github.com/mwatelescope/birli">
  <img src="https://rust-reportcard.xuri.me/badge/github.com/mwatelescope/birli" alt="rust-reportcard"></a>
<a href="https://github.com/MWATelescope/Birli/blob/main/LICENSE">
  <img alt="Crates.io" src="https://img.shields.io/crates/l/birli"></a>

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

- A Rust compiler with a version >= 1.58.0 - <https://www.rust-lang.org/tools/install>
- [AOFlagger](https://gitlab.com/aroffringa/aoflagger) >= 3.0
  (Ubuntu > 21.04: apt install aoflagger-dev)
- [CFitsIO](https://heasarc.gsfc.nasa.gov/fitsio/) >= 3.49
  (Ubuntu > 20.10: apt install libcfitsio-dev)
- [LibERFA](https://github.com/liberfa/erfa) >= 1.7.1
  (Ubuntu > 20.04: apt install liberfa-dev)

for OS-specific instructions, check out the [linux](https://github.com/MWATelescope/Birli/blob/main/.github/workflows/linux_test.yml) and [macOS](https://github.com/MWATelescope/Birli/blob/main/.github/workflows/macos_test.yml) CI Scripts; the [Makefile.toml](https://github.com/MWATelescope/Birli/blob/main/Makefile.toml); and the [Dockerfile](https://github.com/MWATelescope/Birli/blob/main/Dockerfile) as these are tested regularly. The instructions below may be updated less frequently, but are better documented.

### (Debian/Ubuntu) Linux Setup

```bash
# Prerequisites for rustup, cargo and cargo-make
sudo apt install -y gcc libssl-dev pkg-config curl unzip wget
# Run the Rustup install script, profile=default, toolchain=stable
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs -sSf | sh -s -- -y
# Cargo make uses Makefile.toml to automate development tasks
cargo install --force cargo-make
# Use multiple cores when compiling C/C++ libraries
export MAKEFLAGS="-j $MAKEFLAGS"
# Install prerequisite C/C++ libraries
cargo make install_deps
# Ensure that rust can find the C/C++ libraries.
# AOFlagger and CFitsIO default to /usr/local/lib,
# however packages installed with apt (LibERFA) end up in /usr/lib/x86_64-linux-gnu/,
# so we need both.
export LD_LIBRARY_PATH="/usr/local/lib/:/usr/lib/x86_64-linux-gnu/"
```

### MacOS Setup

```bash
# Install homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
# Run the Rustup install script, profile=default, toolchain=stable
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs -sSf | sh -s -- -y
# Cargo make uses Makefile.toml to automate development tasks
cargo install --force cargo-make
# Add the MWATelescope homebrew tap
brew tap mwaTelescope/tap
# Install prerequisite libraries
brew cask install casacore-data casacore aoflagger erfa
```

### Windows Setup

Unfortunately most of the prerequisites aren't available on Windows. However, WSL is great, and there is a docker image! You could use VSCode remote for WSL or Docker. Your best best is Ubuntu LTS

### Installing the binary

```bash
cargo --features=aoflagger install --path .
```

This creates a `birli` binary with the `aoflagger` feature enabled in `$HOME/.cargo/bin`

## Troubleshooting

### Test suite

Having issues with Birli? run the test suite to narrow down your issue.

```bash
cargo test
```

### Dependencies

Experiencing segfaults? I can guarantee it's because of one of the C library dependencies.
Make sure you have the right versions of all the libraries. These are specified in [Prerequisites](#Prerequisites).

Get library versions on linux with:

```bash
pkg-config --modversion cfitsio
pkg-config --modversion erfa
aoflagger --version
```

If you have something like CASA installed from apt, it's going to put an
ancient cfitsio library version in `/usr/lib/x86_64-linux-gnu/`, to get around
this, you must export `LD_LIBRARY_PATH=/usr/local/lib/` in the shell so that Birli can find the correct library version.

### Logging

You can enable additional logging on individual Rust modules by setting the `RUST_LOG` environment variable. For example:

```bash
RUST_LOG=trace birli ... # set log level to trace for all module (including dependencies)
RUST_LOG=birli=debug birli ... # set log level to debug for birli only
RUST_LOG=birli::io=error birli ... # only show warnings for birli's io module
```

For more examples, see the [env_logger docs](https://docs.rs/env_logger/latest/env_logger/)

The default log level in `info`

### Docker

Couldn't get it working on your environment? You can always run Birli in Docker

```bash
docker run mwatelescope/birli:latest -h
```

Want to open a shell within a fully provisioned Birli development environment? Easy!

```bash
docker run -it --entrypoint /bin/bash --volume $PWD:/app mwatelescope/birli:latest
```

Note: This mounts the current directory to `/app` in the Docker image, meaning both of these systems share the same
`target` folder. so if your host system is a different
architecture than Docker, you may need to `cargo clean` each time you switch between these environments. You may also want to temporarily disable any linters or language servers that use

### Singularity on HPC

```bash
# - load the singularity module
module load singularity
# - cd into your preferred sif file location, e.g. /pawsey/mwa/singularity/birli
# - create a .sif file from the latest mwatelescope/birli docker image
singularity pull --dir . docker://mwatelescope/birli:latest
# - run birli within the singularity image
singularity exec  /pawsey/mwa/singularity/birli/birli_latest.sif /app/target/release/birli ${YOUR_BIRLI_ARGS}
```

### Singularity on HPC (debug mode)

This will give you much more information about any problem you're having with Birli, however the
debug build is not optimised, and is much slower.

```bash
# - request an interactive HPC session
salloc --partition workq --time 1:00:00 --nodes 1 -c 38 --mem=350G
# - load the singularity module
module load singularity
# - cd into your preferred sif file location, e.g. /pawsey/mwa/singularity/birli
# - create a .sif file from the latest mwatelescope/birli docker image
singularity pull --dir . docker://mwatelescope/birli:debug
# - run birli within the singularity image
singularity exec  /pawsey/mwa/singularity/birli/birli_debug.sif /bin/bash
```

then within this shell

```bash
# - enable lots of logs
export RUST_LOG=trace
# - run birli in debug mode with GDB
gdb --args /app/target/debug/birli ${YOUR_BIRLI_ARGS}
# > run
```

## Usage

`birli -h`

```txt
USAGE:
    birli [OPTIONS] --metafits <PATH> <PATHS>...

OPTIONS:
        --apply-di-cal <PATH>        Apply DI calibration solutions before averaging
        --dry-run                    Just print the summary and exit
        --emulate-cotter             Use Cotter's array position, not MWAlib's
    -h, --help                       Print help information
        --no-draw-progress           do not show progress bars
        --phase-centre <RA> <DEC>    Override Phase centre from metafits (degrees)
        --pointing-centre            Use pointing instead phase centre
    -V, --version                    Print version information

INPUT:
    -m, --metafits <PATH>    Metadata file for the observation
    <PATHS>...           GPUBox files to process

SELECTION:
        --no-sel-autos            [WIP] Deselect autocorrelations
        --no-sel-flagged-ants     [WIP] Deselect flagged antennas
        --sel-ants <ANTS>...      [WIP] Antenna to select
        --sel-time <MIN> <MAX>    Timestep index range (inclusive) to select

RESOURCE LIMITS:
        --max-memory <GIBIBYTES>    [WIP] Estimate --time-chunk with <GIBIBYTES> GiB each chunk.
        --time-chunk <STEPS>        [WIP] Process observation in chunks of <STEPS> timesteps.

FLAGGING:
        --flag-antennas <ANTS>...         [WIP] Flag antenna indices
        --flag-autos                      [WIP] Flag auto correlations
        --flag-coarse-chans <CHANS>...    [WIP] Flag additional coarse chan indices
        --flag-dc                         [WIP] Force flagging of DC centre chans
        --flag-edge-chans <COUNT>         [WIP] Flag <COUNT> fine chans on the ends of each coarse
        --flag-edge-width <KHZ>           [WIP] Flag bandwidth [kHz] at the ends of each coarse chan
        --flag-end <SECONDS>              [WIP] Flag seconds before the last provided time
        --flag-end-steps <COUNT>          [WIP] Flag <COUNT> steps before the last provided
        --flag-fine-chans <CHANS>...      [WIP] Flag fine chan indices in each coarse chan
        --flag-init <SECONDS>             [WIP] Flag <SECONDS> after first common time (quack time)
        --flag-init-steps <COUNT>         [WIP] Flag <COUNT> steps after first common time
        --flag-times <STEPS>...           [WIP] Flag additional time steps
        --no-flag-dc                      [WIP] Do not flag DC centre chans
        --no-flag-metafits                [WIP] Ignore antenna flags in metafits

CORRECTION:
        --no-cable-delay           Do not perform cable length corrections
        --no-digital-gains         Do not perform digital gains corrections
        --no-geometric-delay       Do not perform geometric corrections
        --passband-gains <TYPE>    Type of PFB passband filter gains correction to apply [default:
                                   jake] [possible values: none, cotter, jake]

AVERAGING:
        --avg-freq-factor <FACTOR>    Average <FACTOR> channels per averaged channel
        --avg-freq-res <KHZ>          Frequency resolution of averaged data
        --avg-time-factor <FACTOR>    Average <FACTOR> timesteps per averaged timestep
        --avg-time-res <SECONDS>      Time resolution of averaged data

OUTPUT:
    -f, --flag-template <TEMPLATE>    The template used to name flag files. Percents are substituted
                                      for the zero-prefixed GPUBox ID, which can be up to 3
                                      characters long. Example: FlagFile%%%.mwaf
    -M, --ms-out <PATH>               Path for measurement set output
    -u, --uvfits-out <PATH>           Path for uvfits output

AOFLAGGER:
        --aoflagger-strategy <PATH>    Strategy to use for RFI Flagging
        --no-rfi                       Do not perform RFI Flagging with aoflagger
```

Note: the aoflagged options are only available when the aoflagger feature is enabled.

Operations are performed in the order described by the following sections.

### Cable Delay Corrections

Cable delay correction involves adjusting visibility phases to correct for the differences in electrical length of the cable between each tile and it's receiver.

Legacy MWA correlator observations do not typically have cable delays applied, however MWAX observations can. The [`CABLEDEL`](https://wiki.mwatelescope.org/display/MP/MWAX+Metafits+Changes) key in the metafits describes what geometric delays have been applied.

By default, Birli will apply cable length corrections. You can use `--no-cable-delay` to disable this.

A baseline's cable lengths are determined by the difference between a baseline's rfInput electrical lengths, as specified the the `TILEDATA` HDU of the metafits. Complex visibilities are phase-shifted by an angle determined by the electrical length, and the channel's frequency.

```rust
let angle = -2.0 * PI * electrical_length_m * freq_hz / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
```

### Digital Gain Corrections

Each input in the raw data is scaled by a factor for each coarse channel. This is defined in the metafits primary hdu in the Gains column. Birli corrects these digital gains by default, you can disable this with `--no-digital-gains`

### Coarse PFB Passband Corrections

Adjust each coarse channel within a fine channel to correct for the shape of the pfb passband curve. Birli will apply the gains defined in the mwa wiki [on pfb gains](https://wiki.mwatelescope.org/display/MP/RRI+Receiver+PFB+Filter)   by default. They can be disabled with `--passband-gains none`. Another option is to emulate Cotter's `_sb128ChannelSubbandValue2014FromMemo` from `subbandpassband.cpp`, sometimes referred to as Levine Gains. Since these gains were computed at the base legacy correlator resolution of 10KHz, they will not work on all MWAX resolutions. Cotter's implementation of this functionality is slightly different, in that it does not include the channel from the gains when scaling. It's not clear if this is a bug or a feature.

When applying pfb gains to an observation that is not at the same resolution as the gains, the gains need to be averaged to fit the data, and the exact details of this averaging depends on the correlator type. For more dtails, see the mwa wiki on [averaging fine channels](https://wiki.mwatelescope.org/display/MP/MWA+Fine+Channel+Centre+Frequencies)

### RFI Flagging.

By default, Birli will flag the data using the default MWA strategy in AOFlagger. You can use the
`--no-rfi` option to disable this, or the `--aoflagger-strategy` option to proived your own strategy
file.

### Geometric Delay Corrections (AKA Phase Tracking)

Geometric correction involves adjusting visibility phases to correct for the differences in distance that light from the phase center has to travel to reach each tile.

Legacy MWA correlator observations are not typically phase tracked, however MWAX observations can have phase tracking applied. The [`GEODEL`](https://wiki.mwatelescope.org/display/MP/MWAX+Metafits+Changes) card in the metafits describes what geometric delays have been applied.

By default, Birli will apply geometric corrections at the phase center if they have not already been applied. It determines the observations phase center from the [`RAPHASE` and `DECPHASE`](https://wiki.mwatelescope.org/display/MP/Metafits+files) cards in the metafits. If these are not available, the pointing center cards ([`RA` and `DEC`](https://wiki.mwatelescope.org/display/MP/Metafits+files)) from the metafits are used. You can use `--no-geometric-delay` to disable geometric corrections, as well as the `--phase-centre` and `--pointing-centre` options to override the phase center.

A baseline's geometric length is determined by the w component of it's UVW fourier-space vector, after applying precession and nutation to it's tiles' positions and the phase center to the J2000 epoch, accounting for stellar aberration. Complex visibilities are phase-shifted by an angle determined by the w-component, and the channel's frequency.

```rust
let angle = -2.0 * PI * uvw.w * freq_hz / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
```

### Calibration

Birli can apply direction independent calibration solutions using the `--apply-di-cal` flag. Solutions are applied before averaging. The number of channels in the un-averaged visibilities must be an integer multiple of the number of channels in the calibration solutions file. Unlike Cotter, Birli will handle calibration solutions where a `NaN` value is present by flagging any visibilities where a NaN is present.

Currently, only the MWA aocal format (.bin), historically generated by the `calibrate` binary in the `mwa-reduce` package is supported. This format is described [here](https://github.com/MWATelescope/cotter/blob/master/solutionfile.h), however due to the ambiguous definition of the startTime and endTime fields, their values are ignored and so only a single timeblock of solutions can be applied.

### Cotter Emulation

The `--emulate-cotter` flag ensures that outputs match Cotter as much as possible. You should only use this flag if you need to perform a direct comparison with Cotter.

By default, Birli will use the MWA array position from MWALib in order to calculate UVWs and geometric corrections. This is more accurate than the one that Cotter uses, and is the main source of error when doing direct comparisons.

This flag is used as part of the tests in `src/main.rs` to validate that Birli's output matches that of Cotter to within an acceptable margin.

### Averaging

To average the data in time or frequency by a given whole number factor, you can provide the `--avg-time-factor`
or `--avg-freq-factor` options. This can also be achieved with the `--avg-time-res` and
`--avg-freq-res` options which take a duration \[seconds\] or ammount of bandwidth \[kHz\]
respectively. This second group of options will choose the closest whole number averaging factor
based on the resolution of the input data.

### Output

Birli can output visibility data to uvfits or measurement set with `--ms-out` (`-M`) or
`--uvfits-out` (`-u`). It can also output flags for each coarse channel in .mwaf format with
`--flag-template` (`-f`), where the `%` characters in the template argument are replaced with
the same zero-prefixed coarse channel identifiers that are used to identify the coarse channel
GPUBox files that the coarse channel data came from. For legacy data, use two percentage characters,
since the coarse channel identifier is the GPUBox number. However, for MWAX data, the coarse channel
identifier is the channel number, which needs three digits.

### Comparison with Cotter

The following table shows how Birli options map onto Cotter options:

| **Birli**                           | **Cotter**              | **Cotter Description**
| ----------------------------------- | ----------------------- | ------
| `--version`                         | `-version`              | Output version and exit.
| `-m <PATH>`                         | `-m <filename>`         | Read meta data from given fits filename.
| `-f`,`-u`,`-M`                      | `-o <filename>`         | Save output to given filename
| `--no-rfi`                          | `-norfi`                | Disable RFI detection.
| `--aoflagger-strategy <PATH>`       | `-flag-strategy <file>` | Use the specified aoflagger strategy.
| `--no-cable-delay`                  | `-nocablelength`        | Do not perform cable length corrections.
| `--no-geom`                         | `-nogeom`               | Disable geometric corrections.
| `--phase-centre <RA> <DEC>`         | `-centre <ra> <dec>`    | Set alternative phase centre, e.g. -centre 00h00m00.0s 00d00m00.0s.
| `--pointing-centre`                 | `-usepcentre`           | Centre on pointing centre.
| `--avg-time-res <SECONDS>`          | `-timeres <s>`          | Average nr of sec of timesteps together before writing to measurement set.
| `--avg-freq-res <KHZ>`              | `-freqres <kHz>`        | Average kHz bandwidth of channels together before writing to measurement set.
| `--apply-di-cal <PATH>`             | `-full-apply <file>`    | Apply a solution file before averaging.
| `--no-digital-gains`                | `-nosbgains`            | Do not correct for the digital gains.
| `--max-memory` (WIP)                | `-absmem <gb>`          | Use at most the given amount of memory, specified in gigabytes.
| `--flag-edge-width <kHz>` (WIP)     | `-edgewidth <kHz>`      | Flag the given width of edge channels of each sub-band (default: 80 kHz).
| `--flag-init <sec>` (WIP)           | `-initflag <sec>`       | Specify number of seconds to flag at beginning of observation (default: QUACK)
| `--flag-end <sec>` (WIP)            | `-endflag <sec>`        | Specify number of seconds to flag extra at end of observation (default: 0s).
| `--flag-dc` (WIP)                   | `-flagdcchannels`       | Flag the centre channel of each sub-band (currently the default).
| `--no-flag-dc` (WIP)                | `-noflagdcchannels`     | Do not flag the centre channel of each sub-band.
| `--flag-antennae <ANTS>...` (WIP)   | `-flagantenna <lst>`    | Mark the comma-separated list of zero-indexed antennae as flagged antennae.
| `--flag-coarse-chans <CHANS>` (WIP) | `-flagsubband <lst>`    | Flag the comma-separated list of zero-indexed sub-bands.
| `--no-sel-autos` (WIP)              | `-noautos`              | Do not output auto-correlations.
| (not `--flag-autos`)                | `-noflagautos`          | Do not flag auto-correlations (default for uvfits file output).
| (default)                           | `-nostats`              | Disable collecting statistics (default for uvfits file output).
| (not `--no-sel-flagged-ants`, WIP)  | `-noantennapruning`     | Do not remove the flagged antennae.
| (default)                           | `-allowmissing`         | Do not abort when not all GPU box files are available (default is to abort).

Birli will eventually perform all the same default preprocessing steps as Cotter when no flags are provided. The exceptions are that we have not yet implemented flagging of edge / centre fine channels / quack timesteps / auto-correlations, pruning of flagged antennas. This means that `birli <in/out args>` is equivalent to:

```bash
 cotter \
  -allowmissing \
  -edgewidth 0 \
  -endflag 0 \
  -initflag 0 \
  -noantennapruning \
  -noflagautos \
  -noflagdcchannels \
  -nostats \
  -flag-strategy <mwa default aoflagger strategy>
  <in/out args>
```

There is no intention of replicating the following options Birli at this point, so please open an issue if these are important to you:

- Coarse channel selection (`-sbcount`, `-sbstart`): This can be done by simply changing which coarse channel files are given in the CLI arguments)
- Dysco compression (`-use-dysco`, `-dysco-config`)
- Manual metadata specification (`-a`, `-h`, `-i`): This information is readily available from metafits.
- `-offline-gpubox-format`
- Quality statistics (`-saveqs`, `-histograms`, `-skipwrite`, `-nostats`)
- `-noflagmissings`: If an HDU is missing, it should always be flagged.
- `-apply`: only `-full-apply` is supported.
- `-noalign`: gpuboxes are always aligned.
- CPU limit (`-j`): Birli uses crossbeam for concurrency which intelligency uses the compute resources available. Strict resource limits can be achieved with cgroups.
- Memory percentage limit (`-mem`): Only `-absmem` is supported. Determining memory limits on HPC systems is unreliable, so we recommend manually specifying a memory limit instead.
- `-sbpassband <file>`
- `-flagfiles <name>` apply existing flags

### Example: RFI Flagging, corrections, averaging, output

In this example, we use the aoflagger subcommand to:

- Perform RFI flagging using the MWA-default flagging strategy
- Perform geometric and cable length corrections
- average the data to 4 seconds, 160khz
- disable passband correction
- Output visibilities to .uvfits (`-u`)

```bash
birli \
  -m tests/data/1254670392_avg/1254670392.metafits \
  -f "/tmp/Flagfile.Birli.MWA.%%.mwaf" \
  -u "/tmp/1254670392.birli.uvfits" \
  --avg-time-res 4 --avg-freq-res 160 \
  tests/data/1254670392_avg/1254670392_*gpubox*.fits
```

The equivalent Cotter
commands would be:

```bash
# output uvfits
cotter \
  -m tests/data/1254670392_avg/1254670392.metafits \
  -o "tests/data/1254670392_avg/1254670392.cotter.uvfits" \
  -allowmissing \
  -edgewidth 0 \
  -endflag 0 \
  -initflag 0 \
  -noantennapruning \
  -noflagautos \
  -noflagdcchannels \
  -nostats \
  -timeres 4 \
  -freqres 160 \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1254670392_avg/1254670392_20191009153257_gpubox*.fits
```

## Contributing

Pull requests are welcome! Please do your best to ensure that the high standards
of test coverage are maintained.

Before each commit, use `cargo make ci` to ensure your code is formatted correctly.

## Acknowledgement

This scientific work uses data obtained from the Murchison Radio-astronomy Observatory. We
acknowledge the Wajarri Yamatji people as the traditional owners of the Observatory site.

## Coverage

<a href="https://codecov.io/gh/MWATelescope/Birli"><img src="https://codecov.io/gh/MWATelescope/Birli/branch/main/graphs/sunburst.svg" height="200px" alt="Birli logo"></a>

This repo is approved by...

<img src="https://github.com/MWATelescope/Birli/raw/main/img/CIRA_Rust_Evangelism_Strike_Force.png" height="200px" alt="CIRA Rust Evangelism Strike Force logo">
