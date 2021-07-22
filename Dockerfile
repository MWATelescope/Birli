FROM ubuntu:20.04

RUN export DEBIAN_FRONTEND=noninteractive \
    && apt-get update \
    && apt-get install -y \
        build-essential \
        casacore-data \
        casacore-dev \
        cmake \
        curl \
        git \
        libblas-dev \
        libboost-date-time-dev \
        libboost-filesystem-dev \
        libboost-system-dev \
        libboost-test-dev \
        liberfa-dev \
        libfftw3-dev \
        libgsl-dev \
        libgtkmm-3.0-dev \
        liblapack-dev \
        liblua5.3-dev \
        libpng-dev \
        libpython3-dev \
        libssl-dev \
        libxml2-dev \
        pkg-config \
        python3 \
        unzip \
        wget \
        zip

# Install CFitsIO

RUN cd /tmp \
    && wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-3.49.tar.gz \
    && tar -zxvf cfitsio-3.49.tar.gz \
    && cd cfitsio-3.49/ \
    && CFLAGS="-O3" ./configure --prefix=/usr/local --enable-reentrant --enable-ssse3 --enable-sse2 \
    && make -j \
    && make install

# Install AOFlagger

RUN cd /tmp \
    && git clone --recurse-submodules https://gitlab.com/aroffringa/aoflagger.git --branch v3.1.0 \
    && cd aoflagger \
    && chmod a+rwx . \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make -j \
    && make install

# Get Rust
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y

# Add .cargo/bin to PATH
ENV PATH="/root/.cargo/bin:${PATH}"

# Get cargo make
RUN cargo install --force cargo-make

ADD . /app
WORKDIR /app

ENV LD_LIBRARY_PATH=/usr/local/lib/:/usr/lib/x86_64-linux-gnu/

RUN cargo install --path .

ENTRYPOINT [ "/root/.cargo/bin/birli" ]