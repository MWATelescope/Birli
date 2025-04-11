<!-- markdownlint-disable=MD025 -->

# Version 0.17.1 (2025-04-11)

- ➕ marlu 0.16.1

# Version 0.17.0 (2025-04-10)

- ➕ marlu 0.16.0 (with mwalib 1.8.7 to support observations with large numbers of tiles/antennas)
- 🙏 quality of life:
  - Fixed clippy lints

# Version 0.16.1 (2025-04-02)

- ✨ new features:
  - `--no-sel-flagged-ants`
  - `--sel-ants`
  - `--no-sel-autos`
- ➕ marlu 0.15.1

# Version 0.16.0 (2024-11-13)

- ✨ new features:
  - add `--van-vleck` to correct for legacy correlator quantization
- ➕ update mwalib 1.8.2 via Marlu 0.15
- ➕ fitsio 0.21.6

# Version 0.15.1 (2024-09-29)

- add read and write rate

# Version 0.15.0 (2024-09-28)

- add `--provided-chan-ranges` to fix #149

# Version 0.14.0 (2024-09-23)

- 🐛 fix linking issues on macOS with aoflagger_sys 0.1.2
- ➕ update mwalib 1.5.0 via Marlu 0.14.0
- ➕ update shlex

# Version 0.13.0 (2024-08-14)

- 🐛 fix issues compiling on arm64:
  - update Marlu 0.13.0

# Version 0.12.0 (2024-06-27)

- 🐛 bug fixes:
  - #152 default pfb gains should be cotter2014 for legacy correlator, not jake
- 🙏 quality of life:
  - #150 implement antenna (tile) selection (`--sel-ants`)
- ➕ dependencies:
  - update default aoflagger version in makefile and dockerfile

# Version 0.11.0 (2024-06-19)

- ➕ dependencies:
  - Marlu 0.11.0
  - built ~0.7.3

# Version 0.10.0 (2023-08-11)

- 🙏 quality of life:
  - uvfits are now written out with a second DATE and INTTIM
- ➕ dependencies:
  - use Marlu 0.10.1
  - use rust version 1.64

# Version 0.9.2 (2023-07-18)

- 🙏 quality of life:
  - update modtime when writing ms

# Version 0.9.1 (2023-03-10)

- 🙏 quality of life:
  - Fix #140 - observations which were scheduled with a contiguous band of coarse channels are
    written to a single file (with any missing channels flagged).
  - `--sel-chan-ranges` can select individual coarse channels in addition to ranges
  - Fix #133 - print more details in impl Display for errors.
  - Fix display of selected and flagged coarse channel table

# Version 0.9.0 (2023-02-17)

- ✨ new features:
  - support for non-contiguous coarse channel selections (picket fence!) 🎉
- 🙏 quality of life:
  - Remove the `--ignore-dut1` flag
    - The DUT1 of an observation is now written as the `UT1UTC` key in a
      measurement set. It continues to be written to `UT1UTC` in uvfits files.
      As this value has no effect on the rest of the visibility file, the
      `--ignore-dut1` flag is now useless.
  - Clippy lints
  - Update Marlu
    - The ERFA C library is no longer needed, so all mentions of it and
      installations of it have been removed.
  - Use a newer aoflagger_sys crate
    - This allows for AOFLAGGER_INCLUDE_DIR and AOFLAGGER_LIB to be set and
      point to where aoflagger is installed, rather than searching in a default
      location
  - report error messages properly, not as Rust-internal representation

# Version 0.8.0 (2022-08-24)

- ✨ new features:
  - use DUT1 functionality as part of Marlu 0.8.0
- ➕ dependencies:
  - use Marlu 0.8.0
  - MSRV 1.63

# Version 0.7.0 (2022-08-04)

- ⚡ performance:
  - Fix unnecessary allocation for smaller visibility chunks to improve memory performance
  - better uvfits performance from marlu 0.7.0
  - mwaf files are written in parallel
- 🙏 quality of life:
  - write `history` metadata in ms and uvfits
  - mwaf files
    - are overwritten if they already exist (previous versions simply panicked)
    - can be written out when chunking the observation
    - are now version 2.0 of the format
    - `GPSTIME` -> `OBSID` (less ambiguous)
    - add `GPSSTART` (the centroid timestep of the first scan of flags)
    - `COTVER` -> `SOFTWARE`
    - remove `COTVDATE`
    - write `AO_VER` (AOFlagger version used)
    - write `AO_STRAT` (AOFlagger strategy file used)
    - write `CMDLINE`
    - add `CH_OCC` (channel occupancy), `BL_OCC` (baseline occupancy) and
      `TILES` (unflagged tile information) HDUs
  - PFB gains are now borrowed instead of cloned
  - log library version info, and applied delay types
- 🏗 api changes:
  - use array views for flagging
  - mwaf files have a new interface
- ➕ dependencies:
  - use Marlu 0.7.2
  - use rust version 1.60
  - use shlex for command line argument escaping

# Version 0.6.4 (2022-05-04)

- ✨ new features:
  - Automatic quack time flagging
  - CLI options for timestep flagging

# Version 0.6.3 (2022-04-20)

- ✨ new features:
  - --flag-edge-chans and --flag-edge-width for edge channel flagging

# Version 0.6.2 (2022-04-11)

- ✨ new features:
  - Automatic DC channel flagging for old correlator inputs
  - --flag-dc and --no-flag-dc control DC flagging explicitly

# Version 0.6.1 (2022-04-04)

- ➕ dependencies:
  - Use Marlu 0.6.1

# Version 0.6.0 (2022-03-24)

- 🏗 api changes:
  - refactor preprocessing into PreprocessingContext
  - refactor context_to_jones_array into VisSelection::read_mwalib
  - refactor flags into FlagContext
  - refactor io ino IOContext
  - refactor main into cli::BirliContext::{from_args, run}
  - refactor pfb gains to take ScrunchType enum
  - improve error handling in cli and main
  - separate allocation and filling of vis, weights, flags
  - bake flags into weights before writing for better performance
  - move uvfits and selection into Marlu
  - rename `--flag-antennae` to `--flag-antennas`
- ➕ dependencies:
  - Use Marlu 0.6.0

# Version 0.5.2 (2022-03-17)

- 🐛 bug fixes:
  - fix <https://github.com/MWATelescope/Birli/issues/58>
- ➕ dependencies:
  - Bump clap from 3.0.14 to 3.1.5

# Version 0.5.1 (2022-03-03)

- ✨ new features:
  - [EOR-71] Correction of coarse pfb passband gains is enabled by default, use `--passband-gains none` to disable it.
- 🙏 quality of life:
  - shows the total time elapsed at the end.

# Version 0.5.0 (2022-02-15)

- ✨ new features:
  - [EOR-40] Correction of digital gains is enabled by default, use `--no-digital-gains` to disable it.
- 🐛 bug fixes:
  - fix [https://github.com/MWATelescope/Birli/issues/44](#44 - Can't apply calibration solutions to partial observation)
- ➕ dependencies:
  - Bump Marlu from 0.4.0 to 0.5.0

# Version 0.4.1 (2022-02-09)

- ✨ new features:
  - --apply-di-cal applies direction independent calibration solutions
  - read mwa ao cal calibration solution files
- 🙏 quality of life:
  - add timing info to logs
  - default to info log level
- 🏗 api changes:
  - deprecate expand_flag_array, use add_dimension instead
- ➕ dependencies
  - Bump crossbeam-utils from 0.8.6 to 0.8.7
  - Bump clap from 3.0.13 to 3.0.14
  - Bump approx from 0.5.0 to 0.5.1
  - Bump crossbeam-utils from 0.8.5 to 0.8.6
  - Bump tempfile from 3.2.0 to 3.3.0

# Version 0.4.0 (2022-02-01)

- Preprocess data in timestep chunks based on maximum memory or number of timesteps.
- --time-sel selects only specific timestep indices.
- Log info about chunks, memory, build version.
- Use Rust 1.58 and 2021 edition

# Version 0.3.1 (2022-01-21)

- expose optional Marlu features
- fix thread 'main' panicked at '`flag-timesteps`.

# Version 0.3.0 (2022-01-20)

- [EOR-5] Spectral and temporal averaging.
- add untested options for basic flagging and phase centre configuration.
- display human-readable information how the input channels and timesteps will be processed.
- Update Marlu

# Version 0.2.1 (2021-12-15)

- upgrade ndarray and marlu for better ms writing performance.

# Version 0.2.0 (2021-11-19)

- Breaking CLI changes, most aoflagger subcommands moved to root
- [EOR-4](https://github.com/MWATelescope/Birli/projects/1#card-73186083) use Marlu to write measurement sets
- write some missing mandatory AIPS 117 fields to uvfits

# Version 0.1.10 (2021-10-22)

- add CLI flag to skip flagging
- calculate the correct reference frequency for uvfits (issue #6)
- extra logging for uvfits

# Version 0.1.9 (2021-09-01)

- [EOR-34] put aoflagger behind feature
- [EOR-33] vis io traits
- slightly faster uvfits

# Version 0.1.8 (2021-08-26)

- [EOR-32] replace aoflagger::imageset with `ndarray::Array3<Jones>` wherever possible.
- use `mwa_rust_core` for positional and Jones
- use mwalib::get_fine_chan_freqs_hz_array

# Version 0.1.7 (2021-08-09)

- Bug fixes for uvfits output and geometric correction
- Cotter emulation for validation of uvfits against cotter output

# Version 0.1.6 (2021-07-22)

- basically just testing that release automation is working
- update all dependencies
- more optimized docker image

# Version 0.1.5 (2021-07-22)

- [EOR-3] generate uvfits file:
  - weights, uvw coordinates and visibilities cross-matched with cotter in automated tests
  - weights aren't yet affected by averaging
- Added Docker image

# Version 0.1.4 (2021-07-09)

- [EOR-22] correct_cable_lengths
- [EOR-18] Implement updated mwalib timestep / coarse channel interface
- [EOR-28] flag antennas flagged in TILEDATA HDU
- [EOR-29] handle misaligned coarse channels, flag missing hdus
- Updates:
  - mwalib = "0.8.3"
  - bindgen = "0.58.1"
- testing:
  - validate flagging output matches cotter exactly
  - utilities for dumping fits files
  - cotter-friendly metafits files
- other:
  - refactor lib,flag_io into flags, corrections, io

# Version 0.1.3 (2021-05-13)

- Release was created to test a fix for docs.rs.

# Version 0.1.2 (2021-05-05)

- Release was created to test release automation.
- [EOR-19] optimizations:
  - reducing `pin()`s.
  - parallelize allocating, loading and flagging
  - use chunks_exact instead of chunks
  - use Vec instead of BTreeMap where possible
- [EOR-21] quality of life features:
  - progress bars!
  - basic timing info for each stage
- [EOR-20] automation and [documentation](doc/benchmark_results.md) of benchmarks
- Public API is fully documented.
- Use Rust 1.51

# Version 0.1.1 (2021-05-04)

- Release was created to test docs.rs. Flagging with AOFlagger mostly works and seems to resemble
  output from Cotter with some caveats:
  - does not match for observations where not all coarse channels start at the same time because of
    <https://github.com/MWATelescope/mwalib/issues/22>
  - no optimization, so significantly slower.
- [EOR-6] CI/CD runs multi-platform tests and tracks coverage automatically
- [EOR-13] 95% test coverage using synthetic test data
- [EOR-14] Created CXX Bindings for AOFlagger
- [EOR-11] GPUFits files can be read into AOFlagger ImageSet objects
- [EOR-15] Writes to .mwaf flag file sets.
- [EOR-17] Implemented `birli aoflagger` command line interface
- [EOR-20] CI/CD tracks benchmark results automatically
