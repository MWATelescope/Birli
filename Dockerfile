FROM python:3.11-slim-bookworm

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    automake \
    build-essential \
    casacore-data \
    casacore-dev \
    clang \
    cmake \
    curl \
    git \
    jq \
    lcov \
    libboost-date-time-dev \
    libboost-filesystem-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-test-dev \
    libcfitsio-dev \
    liberfa-dev \
    libexpat1-dev \
    libfftw3-dev \
    libhdf5-dev \
    liblapack-dev \
    liblua5.3-dev \
    libpng-dev \
    libssl-dev \
    libtool \
    pkg-config \
    unzip \
    zip \
    && rm -rf /var/lib/apt/lists/*

# # Get Rust
ARG RUST_VERSION=stable
ENV RUSTUP_HOME=/opt/rust CARGO_HOME=/opt/cargo
ENV PATH="${CARGO_HOME}/bin:${PATH}"
RUN mkdir -m755 $RUSTUP_HOME $CARGO_HOME && ( \
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | env RUSTUP_HOME=$RUSTUP_HOME CARGO_HOME=$CARGO_HOME sh -s -- -y \
    --profile=minimal \
    --component llvm-tools \
    --default-toolchain=${RUST_VERSION} \
    )

RUN python -m pip install --force-reinstall --no-cache-dir \
    mwalib

# installing aoflagger with `apt install aoflagger-dev` gives weird errors
# looks like -j`nproc` also breaks arm64
ARG AOFLAGGER_BRANCH=v3.4.0
RUN git clone --depth 1 --branch=${AOFLAGGER_BRANCH} --recurse-submodules https://gitlab.com/aroffringa/aoflagger.git /aoflagger && \
    cd /aoflagger && \
    mkdir build && \
    cd build && \
    cmake $CMAKE_ARGS \
    -DENABLE_GUI=OFF \
    .. && \
    make install && \
    ldconfig && \
    cd / && \
    rm -rf /aoflagger
# set up aoflagger python library
ENV PYTHONPATH="/usr/local/lib/"

ADD . /birli
WORKDIR /birli

# e.g. docker build . --build-arg=TEST_SHIM=cargo\ test\ --release
ARG TEST_SHIM=""
RUN ${TEST_SHIM}

RUN cargo install --path . --locked --features=all-static,aoflagger && \
    cargo clean

ENTRYPOINT [ "/opt/cargo/bin/birli" ]
