FROM mwatelescope/mwalib:latest-python3.11-slim-bookworm

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8
ARG DEBUG
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
    liberfa-dev \
    libexpat1-dev \
    libfftw3-dev \
    libhdf5-dev \
    liblua5.3-dev \
    liblapack-dev \
    libpng-dev \
    libssl-dev \
    libtool \
    pkg-config \
    unzip \
    zip \
    && rm -rf /var/lib/apt/lists/*

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
ENV PYTHONPATH="/usr/local/lib/"

ADD . /birli
WORKDIR /birli

ARG TEST_SHIM=""
RUN ${TEST_SHIM}

RUN cargo install --path . --features aoflagger --locked $(test -z "$DEBUG" || echo "--debug") \
    && cargo clean

ENTRYPOINT [ "/opt/cargo/bin/birli" ]
