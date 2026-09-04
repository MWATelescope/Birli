#!/usr/bin/env bash
set -e

echo Cleaning
cargo clean

echo Upadting
cargo update --verbose

echo Cargo check...
cargo check --all-features --all-targets

echo Cargo clippy...
cargo clippy --all-features --all-targets -- -D warnings

echo Cargo fmt...
cargo fmt --check

echo Done
