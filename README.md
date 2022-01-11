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
<a href="https://deps.rs/crate/birli/">
  <img alt="Libraries.io dependency status for GitHub repo" src="https://img.shields.io/librariesio/github/mwatelescope/birli"></a>
<a href="">
  <img alt="Lines of code" src="https://img.shields.io/tokei/lines/github/mwatelescope/birli"></a>

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

### the trait bound `Jones<f32>: AbsDiffEq<_>` is not satisfied

if you see an error that looks like this:

```txt
error[E0277]: the trait bound `Jones<f32>: AbsDiffEq<_>` is not satisfied
    --> src/corrections.rs:1029:9
     |
1029 | /         assert_abs_diff_eq!(
1030 | |             *jones_array.get((3, 3, 1)).unwrap(),
1031 | |             &Jones::from([
1032 | |                 Complex::new(rot_1_xx_3_3_re, rot_1_xx_3_3_im),
...    |
1036 | |             ])
1037 | |         );
     | |__________^ the trait `AbsDiffEq<_>` is not implemented for `Jones<f32>`
     |
     = note: this error originates in the macro `abs_diff_eq` (in Nightly builds, run with -Z macro-backtrace for more info)
```

try

```bash
cargo update
cargo update -p approx:0.5.0 --precise 0.4.0
cargo update -p ndarray:0.15.3 --precise 0.14.0
```

## Usage

`birli -h`

```txt
USAGE:
    birli [OPTIONS] --metafits <PATH> <PATHS>...

OPTIONS:
        --emulate-cotter             Use Cotter's array position, not MWAlib's
    -h, --help                       Print help information
        --no-cable-delay             Do not perform cable length corrections
        --no-geometric-delay         Do not perform geometric corrections
        --phase-centre <RA> <DEC>    Override Phase centre from metafits (degrees)
        --pointing-centre            Use pointing instead phase centre
    -V, --version                    Print version information

INPUT:
    -m, --metafits <PATH>    Metadata file for the observation
    <PATHS>...           GPUBox files to process

FLAGGING:
        --flag-antennae <ANTS>...         Flag antenna indices
        --flag-coarse-chans <CHANS>...    Flag additional coarse channel indices
        --flag-timesteps <STEPS>...       Flag additional timestep indices
        --no-flag-metafits                Ignore antenna flags in metafits

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

### Cable Delay Corrections

Cable delay correction involves adjusting visibility phases to correct for the differences in electrical length of the cable between each tile and it's receiver.

Legacy MWA correlator observations do not typically have cable delays applied, however MWAX observations can. The [`CABLEDEL`](https://wiki.mwatelescope.org/display/MP/MWAX+Metafits+Changes) key in the metafits describes what geometric delays have been applied.

By default, Birli will apply cable length corrections. You can use `--no-cable-delay` to disable this.

A baseline's cable lengths are determined by the difference between a baseline's rfInput electrical lengths, as specified the the `TILEDATA` HDU of the metafits. Complex visibilities are phase-shifted by an angle determined by the electrical length, and the channel's frequency.

```rust
let angle = -2.0 * PI * electrical_length_m * freq_hz / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
```

### Geometric Delay Corrections (AKA Phase Tracking)

Geometric correction involves adjusting visibility phases to correct for the differences in distance that light from the phase center has to travel to reach each tile.

Legacy MWA correlator observations are not typically phase tracked, however MWAX observations can have phase tracking applied. The [`GEODEL`](https://wiki.mwatelescope.org/display/MP/MWAX+Metafits+Changes) card in the metafits describes what geometric delays have been applied.

By default, Birli will apply geometric corrections at the phase center if they have not already been applied. It determines the observations phase center from the [`RAPHASE` and `DECPHASE`](https://wiki.mwatelescope.org/display/MP/Metafits+files) cards in the metafits. If these are not available, the pointing center cards ([`RA` and `DEC`](https://wiki.mwatelescope.org/display/MP/Metafits+files)) from the metafits are used. You can use `--no-geometric-delay` to disable this.

A baseline's geometric length is determined by the w component of it's UVW fourier-space vector, after applying precession and nutation to it's tiles' positions and the phase center to the J2000 epoch, accounting for stellar aberration. Complex visibilities are phase-shifted by an angle determined by the w-component, and the channel's frequency.

```rust
let angle = -2.0 * PI * uvw.w * freq_hz / SPEED_OF_LIGHT_IN_VACUUM_M_PER_S;
```

### Cotter Emulation

The `--emulate-cotter` flag ensures that outputs match Cotter as much as possible. You should only use this flag if you need to perform a direct comparison with Cotter.

By default, Birli will use the MWA array position from MWALib in order to calculate UVWs and geometric corrections. This is more accurate than the one that Cotter uses, and is the main source of error when doing direct comparisons.

This flag is used as part of the tests in `src/main.rs` to validate that Birli's output matches that of Cotter to within an acceptable margin.

### Example: RFI Flagging, corrections and UVFits/mwaf output

In this example, we use the aoflagger subcommand to:

- Perform RFI flagging using the MWA-default flagging strategy
- Perform geometric and cable length corrections
- Output flags to .mwaf (`-f`)
- Output visibilities to .uvfits (`-u`)

```bash
birli \
  -m tests/data/1254670392_avg/1254670392.metafits \
  -f "/tmp/Flagfile.Birli.MWA.%%.mwaf" \
  -u "/tmp/1254670392.birli.uvfits" \
  tests/data/1254670392_avg/1254670392_*gpubox*.fits
```

Since Cotter can't output flags and uvfits at the same time, the equivalent Cotter commands would be:

```bash
# output flags
cotter \
  -m tests/data/1254670392_avg/1254670392.metafits \
  -o "tests/data/1247842824_flags/FlagfileCotterMWA%%.mwaf" \
  -allowmissing \
  -edgewidth 0 \
  -endflag 0 \
  -initflag 0 \
  -noantennapruning \
  -noflagautos \
  -noflagdcchannels \
  -nosbgains \
  -sbpassband tests/data/subband-passband-128ch-unitary.txt \
  -nostats \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1254670392_avg/1254670392_20191009153257_gpubox*.fits
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
  -nosbgains \
  -sbpassband tests/data/subband-passband-128ch-unitary.txt \
  -nostats \
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

## Release Checklist

- [ ] pipeline is green
- [ ] update `RELEASES.md`
- [ ] update `package.version` in `Cargo.toml`
- [ ] `cargo make pre_commit`
- [ ] commit (include Cargo.toml)
- [ ] `git tag -a $tag -m $tag`
- [ ] `git push`
- [ ] `git push --tags`