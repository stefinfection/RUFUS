# RUFUS — canonical container definition (single source of truth).
# SIFs are built from this image via `apptainer pull/build docker://...`; there is no separate .def file.
FROM ubuntu:22.04

# Build-time version, passed by CI from resources/globals.txt (defaults to "dev" for local builds).
ARG RUFUS_VERSION=dev

LABEL author="Stephanie Georges"
LABEL org.opencontainers.image.title="RUFUS"
LABEL org.opencontainers.image.version="${RUFUS_VERSION}"
LABEL org.opencontainers.image.source="https://github.com/stefinfection/RUFUS"

ENV DEBIAN_FRONTEND=noninteractive

# Pinned external tool versions (last reviewed Jan 2025).
ARG HTSLIB_VERSION="1.21"
ARG BAMTOOLS_VERSION="v2.5.2"
ARG BEDTOOLS_VERSION="2.31.1"

# System + build dependencies. Most of these (parallel, gawk, vt, bc, libgsl, the lib*-dev
# packages) are required at RUNTIME by runRufus.sh and the samtools/bcftools stack, so this
# stays a single stage rather than a slimmed multi-stage build.
# Note: the historical Dockerfiles added ppa:ubuntu-toolchain-r/test but never installed a
# newer g++ from it -- the default ubuntu 22.04 g++ 11 builds RUFUS. The PPA was vestigial and
# (under --no-install-recommends) broke the build on a missing gnupg, so it is dropped.
RUN apt-get update && \
    apt-get install -y \
      git cmake wget g++ build-essential zlib1g-dev libbz2-dev bc \
      libgsl0-dev libncurses5-dev autoconf automake make liblzma-dev \
      libcurl4-gnutls-dev libssl-dev vt parallel gawk libjsoncpp-dev \
      libjsoncpp25 curl unzip ca-certificates file && \
    rm -rf /var/lib/apt/lists/*

# AWS CLI — required at runtime: resource_helpers/download_hash.sh fetches region-specific
# 1000G/control exclusion hashes from S3 via `aws s3 ... --no-sign-request`.
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip -q awscliv2.zip && \
    ./aws/install && \
    rm -rf awscliv2.zip aws

# htslib (provides bgzip/tabix + the libhts that samtools/bcftools link against) -> /usr/local.
# Keep the /opt/htslib source tree so samtools/bcftools can build against it via --with-htslib;
# run ldconfig so the runtime linker picks up the freshly installed /usr/local/lib/libhts.
RUN cd /opt && \
    git clone --recurse-submodules https://github.com/samtools/htslib.git --depth 1 --branch "${HTSLIB_VERSION}" && \
    cd htslib && autoreconf -i && ./configure && make && make install && ldconfig

# samtools -> /usr/local (linked against the /opt/htslib source above)
RUN cd /opt && \
    git clone https://github.com/samtools/samtools.git --depth 1 --branch "${HTSLIB_VERSION}" && \
    cd samtools && autoheader && autoconf -Wno-syntax && \
    ./configure --with-htslib=/opt/htslib && make && make install && \
    cd /opt && rm -rf samtools

# bcftools -> /usr/local (last htslib consumer; drop the htslib source tree afterward)
RUN cd /opt && \
    git clone https://github.com/samtools/bcftools.git --depth 1 --branch "${HTSLIB_VERSION}" && \
    cd bcftools && autoheader && autoconf && \
    ./configure --with-htslib=/opt/htslib --enable-libgsl && make && make install && \
    cd /opt && rm -rf bcftools htslib

# bamtools -> /usr/local
RUN cd /opt && \
    git clone https://github.com/pezmaster31/bamtools.git --depth 1 --branch "${BAMTOOLS_VERSION}" && \
    cd bamtools && mkdir build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=/usr/local .. && make && make install && \
    cd /opt && rm -rf bamtools

# bedtools -> /usr/local/bin
RUN cd /opt && \
    wget -q "https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz" && \
    tar -zxf "bedtools-${BEDTOOLS_VERSION}.tar.gz" && \
    cd bedtools2 && make && cp bin/* /usr/local/bin/ && \
    cd /opt && rm -rf bedtools2 "bedtools-${BEDTOOLS_VERSION}.tar.gz"

# RUFUS — built from the checked-out source tree (not a pinned git clone), so each branch/tag
# builds its own code. CMake fetches and builds the bundled externals (modified jellyfish, etc.).
COPY . /opt/RUFUS
RUN cd /opt/RUFUS && mkdir -p bin && cd bin && cmake ../ && make

# Drop the largest build-only packages; keep the rest since runtime depends on them.
RUN apt-get purge -y --auto-remove git wget && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

ENV DEBIAN_FRONTEND=
ENV PATH=/opt/RUFUS/bin:${PATH}
ENV RUFUS_ROOT=/opt/RUFUS
ENV LC_CTYPE=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US.UTF-8

# Smoke test — fails the build if any expected binary or env var is missing.
RUN bash /opt/RUFUS/tests/smoke_test.sh

WORKDIR /data
CMD ["/bin/bash"]
