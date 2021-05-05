<!-- markdownlint-disable=MD025 -->

# Version 0.1.2 (2021-05-05)

* Release was created to test release automation.
* [EOR-19] initial optimization of visibility loading by reducing `pin()`s.
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
