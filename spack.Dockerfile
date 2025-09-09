# ---------------------------
# Builder Stage: Build environment, install dependencies, and generate Spack view
# ---------------------------
FROM spack/ubuntu-jammy:0.23.0 AS builder
# note to replicate some of these steps in your own container, you should first:
# . /opt/spack/share/spack/setup-env.sh

# some packages must be installed by apt in addition to spack.
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt update && apt-get --no-install-recommends install -y \
    python3 \
    # i think aoflagger-sys needs these?
    pkg-config \
    libssl-dev \
    pkg-config \
    # aoflagger needs libpython3-dev?
    # libpython3-dev \
    && apt-get clean

# Clone the custom Spack repository from GitLab
RUN git clone https://github.com/PawseySC/pawsey-spack-config.git /opt/pawsey-spack && \
    spack repo add /opt/pawsey-spack/repo

# Create a new Spack environment which writes to /opt
RUN --mount=type=cache,target=/opt/buildcache \
    mkdir -p /opt/{software,spack_env,view} && \
    spack env create --dir /opt/spack_env && \
    spack env activate /opt/spack_env && \
    spack mirror add --autopush --unsigned mycache file:///opt/buildcache && \
    spack config add "config:install_tree:root:/opt/software" && \
    spack config add "view:/opt/view" && \
    spack add 'rust' && \
    spack install --no-check-signature --fail-fast --no-checksum --jobs=$(nproc) && \
    spack add \
    'birli' \
    && \
    spack spec -I && \
    spack install --no-check-signature --fail-fast --no-checksum --only dependencies --jobs=$(nproc)

ADD . /birli
WORKDIR /birli

RUN --mount=type=cache,target=/opt/buildcache \
    --mount=type=cache,target=/birli/target/ \
    spack env activate /opt/spack_env && \
    cargo test --release --features=all-static,aoflagger && \
    cargo install --path . --locked --features=all-static,aoflagger && \
    cargo clean


# ---------------------------
# Final Stage: Create a lean runtime image
# ---------------------------
# FROM ubuntu:jammy AS runtime

# # Copy necessary files from builder
# COPY --from=builder /opt/software /opt/software
# COPY --from=builder /opt/view /opt/view
# COPY --from=builder /opt/spack_env /opt/spack_env
# COPY --from=builder /opt/spack /opt/spack
# COPY --from=builder /opt/pawsey-spack /opt/pawsey-spack

# # Setup Spack environment
# ENV SPACK_ROOT=/opt/spack \
#     PATH=/opt/view/bin:/opt/software/bin:/usr/local/bin:/usr/bin:/bin
# RUN . /opt/spack/share/spack/setup-env.sh && \
#     spack repo add /opt/pawsey-spack && \
#     spack env activate /opt/spack_env && \
#     echo ". /opt/spack/share/spack/setup-env.sh" >> /etc/profile.d/spack.sh && \
#     echo "spack env activate /opt/spack_env" >> /etc/profile.d/spack.sh && \
#     . /etc/profile.d/spack.sh

# # Create a startup script that activates the environment
# RUN echo '#!/bin/bash' > /usr/local/bin/entrypoint.sh && \
#     echo 'source /opt/spack/share/spack/setup-env.sh' >> /usr/local/bin/entrypoint.sh && \
#     echo 'spack env activate /opt/spack_env' >> /usr/local/bin/entrypoint.sh && \
#     echo 'exec "$@"' >> /usr/local/bin/entrypoint.sh && \
#     chmod +x /usr/local/bin/entrypoint.sh

# # Set the entrypoint to our custom script
# ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
# CMD ["/bin/bash", "-l"]

# # TEST ME :
# # docker build -f spack.Dockerfile . --tag runtime && docker run --rm -it $_ python -m pytest