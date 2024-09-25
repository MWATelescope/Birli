FROM python:3.11-bookworm as base

# suppress perl locale errors
ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8
# suppress apt-get prompts
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    automake \
    build-essential \
    casacore-data \
    casacore-dev \
    clang \
    cmake \
    curl \
    jq \
    lcov \
    libblas-dev \
    libboost-date-time-dev \
    libboost-system-dev \
    libcfitsio-dev \
    liberfa-dev \
    libfftw3-dev \
    libgsl-dev \
    libgtkmm-3.0-dev \
    libhdf5-serial-dev \
    liblapack-dev \
    liblua5.3-dev \
    libpng-dev \
    libssl-dev \
    libtool \
    pkg-config \
    unzip \
    zip \
    2>&1 | tee /birli-apt.txt \
    && \
    apt-get clean all && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get -y autoremove

# Get Rust
ARG RUST_VERSION=stable
ENV RUSTUP_HOME=/opt/rust CARGO_HOME=/opt/cargo
ENV PATH="${CARGO_HOME}/bin:${PATH}"
RUN mkdir -m755 $RUSTUP_HOME $CARGO_HOME && ( \
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | env RUSTUP_HOME=$RUSTUP_HOME CARGO_HOME=$CARGO_HOME sh -s -- -y \
    --profile=minimal \
    --component llvm-tools \
    --default-toolchain=${RUST_VERSION} \
    )
# # Get cargo make, llvm-cov
RUN cargo install --force cargo-make cargo-llvm-cov && \
    rm -rf ${CARGO_HOME}/registry

RUN python -m pip install --no-cache-dir maturin[patchelf]==1.7.0

ARG MWALIB_BRANCH=v1.5.0
RUN git clone --depth 1 --branch=${MWALIB_BRANCH} https://github.com/MWATelescope/mwalib.git /mwalib && \
    cd /mwalib && \
    maturin build --release --features=python && \
    python -m pip install $(ls -1 target/wheels/*.whl | tail -n 1) && \
    cd / && \
    rm -rf /mwalib ${CARGO_HOME}/registry

# for example, CMAKE_ARGS="-D CMAKE_CXX_FLAGS='-march=native -mtune=native -O3 -fomit-frame-pointer'"
ARG CMAKE_ARGS="-D PORTABLE=ON"

ARG AOFLAGGER_BRANCH=v3.4.0
RUN git clone --depth 1 --branch=${AOFLAGGER_BRANCH} --recurse-submodules https://gitlab.com/aroffringa/aoflagger.git /aoflagger && \
    cd /aoflagger && \
    mkdir build && \
    cd build && \
    cmake $CMAKE_ARGS \
    -DENABLE_GUI=OFF \
    .. && \
    make install -j`nproc` && \
    ldconfig && \
    cd / && \
    rm -rf /aoflagger
# set up aoflagger python library
ENV PYTHONPATH="/usr/local/lib/:$PYTHONPATH"

ADD . /app
WORKDIR /app

ARG DEBUG
ARG BIRLI_FEATURES=cli,aoflagger
RUN cargo install --path . --no-default-features --features=${BIRLI_FEATURES} --locked $(test -z "$DEBUG" || echo "--debug") && \
    rm -rf /app/target ${CARGO_HOME}/registry

RUN test -z "$DEBUG" || ( \
    apt-get install -y vim gdb && \
    apt-get clean all && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get -y autoremove \
    )

ARG TEST_SHIM
RUN ${TEST_SHIM}

ENTRYPOINT [ "/opt/cargo/bin/birli" ]
