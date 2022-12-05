FROM debian:bullseye-20220328-slim AS downloader

ARG SENTIEON_VERSION
RUN test -n "$SENTIEON_VERSION"

LABEL container.base.image="debian:bullseye-20220328-slim" \
      software.version="${VERSION}" \
      software.website="https://www.sentieon.com/"

# Download the software from the permalink
RUN apt-get update && apt-get install -y curl && \
    mkdir -p /opt/sentieon/ && \
    curl -L "https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-${SENTIEON_VERSION}.tar.gz" | \
      tar -zxf - -C /opt/sentieon/

# Install igzip for fastq decompression
RUN apt-get update && apt-get install -y curl autoconf libtool binutils make nasm yasm && \
    mkdir -p /opt/isa-l && \
    curl -L "https://github.com/intel/isa-l/archive/refs/tags/v2.30.0.tar.gz" | \
    tar -zxf - -C /opt/isa-l/ && \
    cd /opt/isa-l/isa-l-2.30.0/ && \
    ./autogen.sh && \
    ./configure --prefix=/usr --libdir=/usr/lib && \
    make install

# Install samtools
RUN apt-get update && apt-get install -y bzip2 curl autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev && \
    mkdir -p /opt/samtools && \
    curl -L "https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2" | \
    tar -jxf - -C /opt/samtools && \
    cd /opt/samtools/samtools-1.16.1/ && \
    ./configure && \
    make install

# Build the container
FROM debian:bullseye-20220328-slim
ARG SENTIEON_VERSION
ENV SENTIEON_VERSION=$SENTIEON_VERSION

# Copy dependencies from the first stage
COPY --from=downloader /opt/sentieon/sentieon-genomics-${SENTIEON_VERSION} /opt/sentieon/sentieon-genomics-${SENTIEON_VERSION}
COPY --from=downloader /usr/bin/igzip /usr/bin/igzip
COPY --from=downloader /usr/lib/libisal.a /usr/lib/libisal.a
COPY --from=downloader /usr/lib/libisal.so.2.0.30 /usr/lib/libisal.so.2.0.30
COPY --from=downloader /usr/lib/libisal.so.2 /usr/lib/libisal.so.2
COPY --from=downloader /usr/lib/libisal.so /usr/lib/libisal.so
COPY --from=downloader /usr/lib/libisal.la /usr/lib/libisal.la
COPY --from=downloader /usr/local/bin/samtools /usr/local/bin/samtools

CMD ["/bin/bash"]

# Install jemalloc as the recommended memory allocation tool, see https://support.sentieon.com/appnotes/jemalloc/
# Install procps for process monitoring
# Install samtools shared libraries
RUN apt-get update && apt-get install -y \
    libjemalloc2 \
    procps \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev \
    tini \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

ENV SENTIEON_INSTALL_DIR=/opt/sentieon/sentieon-genomics-$SENTIEON_VERSION
ENV PATH $SENTIEON_INSTALL_DIR/bin/:$PATH
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2

# A default jemalloc configuration that should work well for most use-cases, see http://jemalloc.net/jemalloc.3.html
ENV MALLOC_CONF=metadata_thp:auto,background_thread:true,dirty_decay_ms:30000,muzzy_decay_ms:30000

# Test the container
RUN sentieon driver --help && \
    igzip --help && \
    samtools --help

ENTRYPOINT ["/usr/bin/tini" , "--"]
