FROM ubuntu:21.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y \
        aoflagger-dev \
        build-essential \
        curl \
        git \
        libcfitsio-dev \
        liberfa-dev \
        libssl-dev \
        pkg-config \
        unzip \
        zip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Get Rust
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y

# Add .cargo/bin to PATH
ENV PATH="/root/.cargo/bin:${PATH}"

# Get cargo make
RUN cargo install --force cargo-make

ADD . /app
WORKDIR /app

RUN cargo clean
RUN cargo install --path .

ENTRYPOINT [ "/root/.cargo/bin/birli" ]
