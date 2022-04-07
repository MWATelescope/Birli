<!-- markdownlint-disable=MD025 -->

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
