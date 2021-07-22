<!-- markdownlint-disable=MD025 -->

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

* Release was created to test a fix for docs.rs.

# Version 0.1.2 (2021-05-05)

* Release was created to test release automation.
* [EOR-19] optimizations:
  * reducing `pin()`s.
  * parallelize allocating, loading and flagging
  * use chunks_exact instead of chunks
  * use Vec instead of BTreeMap where possible
* [EOR-21] quality of life features:
  * progress bars!
  * basic timing info for each stage
* [EOR-20] automation and [documentation](doc/benchmark_results.md) of benchmarks
* Public API is fully documented.
* Use Rust 1.51

# Version 0.1.1 (2021-05-04)

* Release was created to test docs.rs. Flagging with AOFlagger mostly works and seems to resemble
  output from Cotter with some caveats:
  * does not match for observations where not all coarse channels start at the same time because of
    <https://github.com/MWATelescope/mwalib/issues/22>
  * no optimization, so significantly slower.
* [EOR-6] CI/CD runs multi-platform tests and tracks coverage automatically
* [EOR-13] 95% test coverage using synthetic test data
* [EOR-14] Created CXX Bindings for AOFlagger
* [EOR-11] GPUFits files can be read into AOFlagger ImageSet objects
* [EOR-15] Writes to .mwaf flag file sets.
* [EOR-17] Implemented `birli aoflagger` command line interface
* [EOR-20] CI/CD tracks benchmark results automatically
