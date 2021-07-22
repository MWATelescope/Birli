# Birli

<!-- markdownlint-disable MD033 -->
<div class="bg-gray-dark" align="center" style="background-color:#24292e">
<img src="img/birli.png" height="200px" alt="Birli logo">
<br/>
<a href="https://github.com/MWATelescope/Birli/actions/workflows/linux_test.yml"><img src="https://github.com/MWATelescope/Birli/actions/workflows/linux_test.yml/badge.svg" alt="Cross-Platform Tests"></a>
<a href="https://github.com/MWATelescope/Birli/actions/workflows/macos_test.yml"><img src="https://github.com/MWATelescope/Birli/actions/workflows/macos_test.yml/badge.svg" alt="Cross-Platform Tests"></a>
<a href="https://crates.io/birli"><img alt="Crates.io" src="https://img.shields.io/crates/d/birli?label=crates.io%20%E2%AC%87%EF%B8%8F"></a>
<a href="https://docs.rs/crate/birli/"><img src="https://docs.rs/birli/badge.svg" alt="codecov"></a>
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
```

### Windows Setup

Unfortunately most of the prerequisites aren't available on Windows. However, WSL is great, and there is a docker image! You could use VSCode remote for WSL or Docker. Your best best is Ubuntu LTS

### Installing the binary

```bash
cargo install --path .
```

This creates a `birli` binary in `$HOME/.cargo/bin`

## Troubleshooting

Having issues with Birli? run the test suite to narrow down your issue.

```bash
cargo test
```

Experiencing segfaults? I can guarantee it's because of one of the C library dependencies.
Make sure you have the right versions of all the libraries. These are specified in [Prerequisites](#Prerequisites). 

Get library versions with:

```bash
apt show liberfa-dev
aoflagger --version
# cfitsio: ???
```


If you have something like CASA installed from apt, it's going to put an 
ancient cfitsio library version in `/usr/lib/x86_64-linux-gnu/`, to get around
this, you must export `LD_LIBRARY_PATH=/usr/local/lib/:/usr/lib/x86_64-linux-gnu/` in the shell so that Birli can find the correct library version.
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
architecture than Docker, you may need to `cargo clean` each time you switch between these environments. You
may also want to temporarily disable any linters or language servers that use 

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
flag visibilities with aoFlagger

USAGE:
    birli aoflagger [FLAGS] [OPTIONS] <fits-files>... -m <metafits>

FLAGS:
    -h, --help              Prints help information
        --no-cable-delay    Do not perform cable length corrections.
    -V, --version           Prints version information

OPTIONS:
    -f <flag-template>        Sets the template used to name flag files. Percents are substituted for the zero-prefixed
                              GPUBox ID, which can be up to 3 characters log. Similar to -o in Cotter. Example:
                              FlagFile%%%.mwaf
    -m <metafits>             Sets the metafits file.
    -u <uvfits-out>           Filename for uvfits output. Similar to -o in Cotter. Example: 1196175296.uvfits

ARGS:
    <fits-files>...  
```

### Examples

A direct comparison between Birli and Cotter might look like this.

```bash
birli aoflagger \
  -m tests/data/1247842824_flags/1247842824.metafits \
  -f "tests/data/1247842824_flags/flags_birli/FlagfileBirliMWA%%.mwaf" \
  tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits
```

```bash
cotter \
  -m tests/data/1247842824_flags/1247842824cotter-friendly.metafits \
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
  -sbpassband tests/data/subband-passband-128ch-unitary.txt \
  -nostats \
  -sbcount 1 \
  -sbstart 1 \
  -flag-strategy /usr/local/share/aoflagger/strategies/mwa-default.lua \
  tests/data/1247842824_flags/1247842824_20190722150008_gpubox01_00.fits
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
