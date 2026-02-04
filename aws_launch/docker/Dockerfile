FROM ubuntu:22.04 AS BUILD

# Metadata labels
LABEL author="Stephanie Georges"
LABEL version="v0.1.0-epsilon"

# Set environment variables for installation
ENV DEBIAN_FRONTEND=noninteractive
ENV HTSLIB_VERSION="1.21"

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
        git \
        cmake \
        wget \
        g++ \
        build-essential \
        software-properties-common \
        zlib1g-dev \
        libbz2-dev \
        bc \
        libgsl0-dev \
        libncurses5-dev \
        autoconf \
        automake \
        make \
        liblzma-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        vt \
        parallel \
        gawk \
        libjsoncpp-dev \
        libjsoncpp25 \
        curl \
        unzip && \
    add-apt-repository ppa:ubuntu-toolchain-r/test

# Install AWS CLI
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm -rf awscliv2.zip aws

# Clean up apt cache
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Note: the following versions were updated Jan 2025

# Install htslib
RUN cd /opt && \
    git clone --recurse-submodules https://github.com/samtools/htslib.git --branch "$HTSLIB_VERSION" && \
    cd htslib && \
    autoreconf -i && \
    ./configure && \
    make && \
    make install

# Install samtools
RUN cd /opt && \
    git clone https://github.com/samtools/samtools.git --branch "$HTSLIB_VERSION" && \
    cd samtools && \
    autoheader && \
    autoconf -Wno-syntax && \
    ./configure --with-htslib=/opt/htslib/ && \
    make && \
    make install

# Install bcftools
RUN cd /opt && \
    git clone https://github.com/samtools/bcftools.git --branch "$HTSLIB_VERSION" && \
    cd bcftools && \
    autoheader && \
    autoconf && \
    ./configure --with-htslib=/opt/htslib/ --enable-libgsl && \
    make && \
    make install

# Install bamtools
RUN cd /opt && \
    git clone https://github.com/pezmaster31/bamtools.git --branch v2.5.2 && \
    cd bamtools/ && \
    mkdir build && \
    cd build/ && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/bamtools/build .. && \
    make && \
    make install

# Install bedtools
RUN cd /opt && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar -zxvf bedtools-2.31.1.tar.gz && \
    cd bedtools2 && \
    make && \
    make install

# Install RUFUS
RUN cd /opt && \
    git clone https://github.com/stefinfection/RUFUS.git -b singularity && \
    cd RUFUS && \
    mkdir -p bin && \
    cd bin && \
    cmake ../ && \
    make

# Clean up unnecessary packages
RUN apt-get purge -y --auto-remove git wget

# Unset DEBIAN_FRONTEND
ENV DEBIAN_FRONTEND=

# Set runtime environment variables
ENV PATH=$PATH:/opt/RUFUS/bin
ENV RUFUS_ROOT=/opt/RUFUS
ENV LC_CTYPE=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US.UTF-8
ENV CC=gcc
ENV CXX=g++

# Health check / test commands
RUN /opt/samtools/samtools || true && \
    /opt/bamtools/build/src/bamtools || true && \
    /opt/bedtools2/bin/bedtools || true && \
    bgzip || true

# Default command
CMD ["/bin/bash"]
