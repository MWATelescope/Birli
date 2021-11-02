FROM ubuntu:21.04

ENV DEBIAN_FRONTEND=noninteractive
ARG DEBUG
RUN apt-get update \
    && apt-get install -y \
        aoflagger-dev \
        build-essential \
        curl \
        git \
        jq \
        libcfitsio-dev \
        liberfa-dev \
        libssl-dev \
        pkg-config \
        unzip \
        zip
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

ENTRYPOINT [ "/app/target/release/birli" ]