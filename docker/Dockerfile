# syntax=docker/dockerfile:1
#
# vpfkit application image.
#
# Base image
# ----------
# bioconductor/bioconductor_docker already carries the system libraries every
# Bioconductor package needs, and its versioned tags publish both linux/amd64
# and linux/arm64 manifests, so this image builds natively on Apple silicon and
# on ARM64 Linux hosts. The previous base, rocker/shiny-verse, publishes
# linux/amd64 only.
#
# Pin the versioned tag, never `latest`: at the time of writing the `latest`
# tag was still amd64-only while RELEASE_3_23 was multi-arch.
#
# Architecture note
# -----------------
# On amd64 the Bioconductor container binary repository is enabled and the
# dependency install takes a few minutes. On arm64 that repository is not
# available and every Bioconductor package is compiled from source, so a clean
# arm64 build takes considerably longer. It works; it is not fast. Use
# BuildKit's cache mount (already configured below) when iterating.
#
# Build:
#   docker build -f docker/Dockerfile -t vpfkit:latest .
# Build for one architecture explicitly:
#   docker build --platform linux/arm64 -f docker/Dockerfile -t vpfkit:latest .
# Publish both architectures:
#   docker buildx build --platform linux/amd64,linux/arm64 \
#     -f docker/Dockerfile -t denglab/vpfkit:latest --push .
# Run:
#   docker run --rm -p 3838:3838 vpfkit:latest
#   # then open http://localhost:3838

FROM bioconductor/bioconductor_docker:RELEASE_3_23-r-4.6.1

LABEL org.opencontainers.image.title="vpfkit" \
      org.opencontainers.image.description="Shiny viewer and toolkit for ViroProfiler results" \
      org.opencontainers.image.source="https://github.com/deng-lab/vpfkit" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

# Quarto CLI, needed by generate_report(). The .deb is published for both
# amd64 and arm64; dpkg --print-architecture picks the right one.
ARG QUARTO_VERSION=1.9.38
RUN set -eux; \
    arch="$(dpkg --print-architecture)"; \
    curl -fsSL -o /tmp/quarto.deb \
      "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${arch}.deb"; \
    dpkg -i /tmp/quarto.deb; \
    rm -f /tmp/quarto.deb; \
    quarto --version

WORKDIR /opt/vpfkit

# Install dependencies from DESCRIPTION first, so that editing package sources
# does not invalidate the expensive dependency layer.
COPY DESCRIPTION ./DESCRIPTION
RUN --mount=type=cache,target=/root/.cache/R/pkgcache \
    R --vanilla -q -e "\
      options(Ncpus = max(1L, parallel::detectCores())); \
      install.packages('pak', repos = sprintf('https://r-lib.github.io/p/pak/stable/%s/%s/%s', .Platform\$pkgType, R.Version()\$os, R.Version()\$arch)); \
      pak::local_install_deps('.', dependencies = TRUE, ask = FALSE)"

COPY . .

RUN R CMD INSTALL --no-multiarch --with-keep.source . && \
    R --vanilla -q -e "library(vpfkit); stopifnot(is.function(vpfkit::run_app))"

# Do not serve the application as root.
RUN useradd --create-home --shell /bin/bash vpfkit
USER vpfkit

EXPOSE 3838

CMD ["R", "--vanilla", "-q", "-e", \
     "options(shiny.host='0.0.0.0', shiny.port=3838); vpfkit::run_app()"]
