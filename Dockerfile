FROM ubuntu:21.04

ENV DEBIAN_FRONTEND=noninteractive
ARG DEBUG
RUN apt-get update \
    && apt-get install -y \
        aoflagger-dev \
        build-essential \
        clang \
        curl \
        git \
        jq \
        lcov \
        libcfitsio-dev \
        liberfa-dev \
        libssl-dev \
        pkg-config \
        unzip \
        zip \

RUN test -z "$DEBUG" || ( \
        apt-get install -y vim gdb \
    )
RUN apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Get Rust
RUN mkdir -m755 /opt/rust /opt/cargo
ENV RUSTUP_HOME=/opt/rust CARGO_HOME=/opt/cargo PATH=/opt/cargo/bin:$PATH
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y

# Get cargo make
RUN cargo install --force cargo-make cargo-binutils

ADD . /app
WORKDIR /app

RUN cargo clean \
    && cargo install --path . --features aoflagger --locked
RUN test -z "$DEBUG" || (\
        mkdir benches \
        && touch benches/bench.rs \
        && cargo build --features aoflagger \
    )

# setup the toolchain used for coverage analysis
RUN rustup toolchain install nightly-2021-05-09 --component llvm-tools-preview --profile minimal \
    && cargo +nightly-2021-05-09 update -p syn --precise 1.0.80 \
    && cargo +nightly-2021-05-09 update -p proc-macro2 --precise 1.0.28 \
    && cargo +nightly-2021-05-09 install --force cargo-make --locked --version '=0.32' \
    && cargo +nightly-2021-05-09 install --force cargo-binutils --locked --version '=0.3.3'

ENTRYPOINT [ "/app/target/release/birli" ]