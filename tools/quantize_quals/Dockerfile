FROM debian:bookworm-slim AS build

ARG HTSLIB_VERSION=1.22.1

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential ca-certificates curl pkg-config \
        zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev \
        samtools \
    && rm -rf /var/lib/apt/lists/*

RUN curl -fsSL "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" \
        | tar -xj -C /tmp \
    && cd /tmp/htslib-${HTSLIB_VERSION} \
    && ./configure --prefix=/opt/htslib --disable-libcurl --disable-gcs --disable-s3 \
    && make -j"$(nproc)" \
    && make install

WORKDIR /src
COPY Makefile ./
COPY src ./src
COPY tests ./tests

ENV PKG_CONFIG_PATH=/opt/htslib/lib/pkgconfig \
    LD_LIBRARY_PATH=/opt/htslib/lib
RUN make test && make install PREFIX=/opt/quantize_quals

FROM debian:bookworm-slim

LABEL org.opencontainers.image.source="https://github.com/break-through-cancer/parabricks-bqsr" \
      org.opencontainers.image.description="quantize_quals: GATK-compatible static quality-score quantization (HTSlib)" \
      org.opencontainers.image.version="0.1.0"

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps zlib1g libbz2-1.0 liblzma5 libdeflate0 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/htslib/lib/ /usr/local/lib/
COPY --from=build /opt/quantize_quals/bin/quantize_quals /usr/local/bin/quantize_quals
RUN ldconfig && quantize_quals --version
