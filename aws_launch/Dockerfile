FROM ubuntu:22.04

# Metadata labels
LABEL author="Stephanie Georges"
LABEL version="v0.1.1-epsilon"

# Internal install location
ARG INSTALL_PREFIX=/opt
ENV INSTALL_PREFIX=${INSTALL_PREFIX}
ENV DEBIAN_FRONTEND=noninteractive

# Versions
# Note: the following versions were updated Jan 2025
ENV HTSLIB_VERSION="1.21"
ENV BAMTOOLS_VERSION="v2.5.2"
ENV RUFUS_VERSION=""

WORKDIR ${INSTALL_PREFIX}

# Install dependencies
RUN apt-get update && apt-get install -y software-properties-common && \
    add-apt-repository ppa:ubuntu-toolchain-r/test && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
      git cmake wget g++ build-essential zlib1g-dev libbz2-dev bc \
      libgsl0-dev libncurses5-dev autoconf automake make liblzma-dev \
      libcurl4-gnutls-dev libssl-dev vt parallel gawk libjsoncpp-dev \
      libjsoncpp25 curl unzip && \
    rm -rf /var/lib/apt/lists/*

# Install AWS CLI
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm -rf awscliv2.zip aws

# Clean up apt cache
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install htslib
RUN git clone --recurse-submodules https://github.com/samtools/htslib.git --depth 1 --branch "${HTSLIB_VERSION}" && \
    cd htslib && \
    autoreconf -i && \
    ./configure --prefix=${INSTALL_PREFIX}/htslib && \
    make && \
    make install && \
    cd .. && rm -rf htslib


# Install samtools
RUN git clone https://github.com/samtools/samtools.git --depth 1 --branch "${HTSLIB_VERSION}" && \
    cd samtools && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure --prefix=${INSTALL_PREFIX}/samtools --with-htslib=${INSTALL_PREFIX}/htslib && \
    make && make install && \
    cd .. && rm -rf samtools

# Install bcftools
RUN git clone https://github.com/samtools/bcftools.git --depth 1 --branch "${HTSLIB_VERSION}" && \
    cd bcftools && \
    autoheader && \
    autoconf && \
    ./configure --prefix=${INSTALL_PREFIX}/bcftools --with-htslib=${INSTALL_PREFIX}/htslib --enable-libgsl && \
    make && make install && \
    cd .. && rm -rf bcftools

# Install bamtools
RUN git clone https://github.com/pezmaster31/bamtools.git --depth 1 --branch "${BAMTOOLS_VERSION}" && \
    cd bamtools && \
    mkdir -p build && cd build && \
    cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX}/bamtools .. && \
    make && make install && \
    cd ../.. && rm -rf bamtools

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar -zxvf bedtools-2.31.1.tar.gz && \
    cd bedtools2 && \
    make && \
    mkdir -p ${INSTALL_PREFIX}/bedtools2 && \
    cp -r bin ${INSTALL_PREFIX}/bedtools2/ && \
    cd .. && rm -rf bedtools-2.31.1.tar.gz bedtools2

# Install RUFUS
RUN git clone https://github.com/stefinfection/RUFUS.git && \
    cd RUFUS && \
    git checkout docker --branch ${RUFUS_VERSION} --depth 1 && \
    mkdir bin && \
    cd bin && \
    cmake ../ && \
    make

# Clean up for run time
RUN apt-get purge -y --auto-remove git wget build-essential cmake && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Unset DEBIAN_FRONTEND
ENV DEBIAN_FRONTEND=

# Set runtime environment variables
ENV PATH=${INSTALL_PREFIX}/RUFUS/bin:${INSTALL_PREFIX}/samtools/bin:${INSTALL_PREFIX}/bamtools/bin:${INSTALL_PREFIX}/bedtools2/bin:${INSTALL_PREFIX}/bcftools/bin:${PATH}
ENV RUFUS_ROOT=${INSTALL_PREFIX}/RUFUS
ENV LC_CTYPE=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US.UTF-8

# Health check / test commands
RUN samtools || true && \
    bamtools || true && \
    bedtools || true && \
    bcftools || true && \
    bgzip || true

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]