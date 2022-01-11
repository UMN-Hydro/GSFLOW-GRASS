FROM ubuntu:16.04
USER root

RUN useradd -m -s /bin/bash dgketchum

RUN apt-get update && apt-get install tzdata
RUN apt-get install -y apt-utils

ARG CURL_VERSION=7.59.0
ARG GDAL_VERSION=2.3.2
ARG LIBJPEG_TURBO_VERSION=1.5.90
ARG NGHTTP2_VERSION=1.32.0
ARG OPENJPEG_VERSION=2.3.0
ARG ZSTD_VERSION=1.3.5

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update \
  && apt-get upgrade -y \
  && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    debhelper \
    dh-autoreconf \
    autotools-dev \
    zlib1g-dev \
    libcrypto++-dev \
    libjasper-dev \
    libpng-dev \
    libgif-dev \
    libwebp-dev \
    libhdf5-dev \
    libpcre3-dev \
    libxerces-c-dev \
    d-shlibs \
    libgeos-dev \
    python-all-dev \
    python-numpy \
    libsqlite3-dev \
    libexpat1-dev \
    libproj-dev \
    libxml2-dev \
    libspatialite-dev \
    liblzma-dev \
    libarmadillo-dev \
    liburiparser-dev \
    pkg-config \
    libgnutls-dev \
    cmake \
    nasm \
  && mkdir /tmp/nghttp2 \
  && curl -sfL https://github.com/nghttp2/nghttp2/releases/download/v${NGHTTP2_VERSION}/nghttp2-${NGHTTP2_VERSION}.tar.gz | tar zxf - -C /tmp/nghttp2 --strip-components=1 \
  && cd /tmp/nghttp2 \
  && ./configure --enable-lib-only \
  && make -j $(nproc) install \
  && mkdir /tmp/curl \
  && curl -sfL https://curl.haxx.se/download/curl-${CURL_VERSION}.tar.gz | tar zxf - -C /tmp/curl --strip-components=1 \
  && cd /tmp/curl \
  && ./configure --prefix=/opt/curl --disable-manual --disable-cookies --with-gnutls \
  && make -j $(nproc) install \
  && mkdir /tmp/zstd \
  && curl -sfL https://github.com/facebook/zstd/archive/v${ZSTD_VERSION}.tar.gz | tar zxf - -C /tmp/zstd --strip-components=1 \
  && cd /tmp/zstd \
  && make -j $(nproc) install \
  && mkdir -p /tmp/libjpeg-turbo \
  && curl -sfL https://github.com/libjpeg-turbo/libjpeg-turbo/archive/${LIBJPEG_TURBO_VERSION}.tar.gz | tar zxf - -C /tmp/libjpeg-turbo --strip-components=1 \
  && cd /tmp/libjpeg-turbo \
  && cmake -G"Unix Makefiles" -DCMAKE_INSTALL_PREFIX=/usr . \
  && make -j $(nproc) install \
  && cd / \
  && rm -rf /tmp/curl /tmp/libjpeg-turbo /tmp/nghttp2 /tmp/zstd \
  && curl -sfL https://github.com/uclouvain/openjpeg/releases/download/v${OPENJPEG_VERSION}/openjpeg-v${OPENJPEG_VERSION}-linux-x86_64.tar.gz | tar zxf - -C /usr/local --strip-components=1 \
  && ldconfig \
  && mkdir -p /tmp/gdal \
  && curl -sfL https://github.com/OSGeo/gdal/archive/v${GDAL_VERSION}.tar.gz | tar zxf - -C /tmp/gdal --strip-components=2 \
  && cd /tmp/gdal \
  && ./configure \
    --prefix=/usr \
    --mandir=/usr/share/man \
    --includedir=/usr/include/gdal \
    --with-threads \
    --with-grass=no \
    --with-hide-internal-symbols=yes \
    --with-rename-internal-libtiff-symbols=yes \
    --with-rename-internal-libgeotiff-symbols=yes \
    --with-libtiff=internal \
    --with-geotiff=internal \
    --with-webp \
    --with-jasper \
    --with-jpeg=/usr \
    --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial \
    --with-xerces \
    --with-geos \
    --with-sqlite3 \
    --with-curl=/opt/curl/bin/curl-config \
    --with-proj=yes \
    --with-spatialite=/usr \
    --with-cfitsio=no \
    --with-ecw=no \
    --with-mrsid=no \
    --with-openjpeg=yes \
    --with-armadillo=yes \
    --with-liblzma=yes \
    --with-zstd \
    --with-cryptopp=yes \
  && make -j $(nproc) \
  && make -j $(nproc) install \
  && cd / \
  && rm -rf /tmp/gdal \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

ENV CURL_CA_BUNDLE /etc/ssl/certs/ca-certificates.crt
ENV GDAL_DATA=/usr/share/gdal
ENV GDAL_HTTP_VERSION 2

ADD install_gsflow_grass /install_gsflow_grass

RUN apt update && \
    chmod +x /install_gsflow_grass && /install_gsflow_grass

ADD install_grass_package /install_grass_package
RUN chmod +x /install_grass_package && /install_grass_package

#ADD install_binaries_grass /install_binaries_grass
#RUN chmod +x /install_binaries_grass && /install_binaries_grass

ADD install_kit /home/dgketchum/install_kit
RUN chmod +x /home/dgketchum/install_kit

RUN apt install -y libgl1-mesa-glx libgl1-mesa-dri
USER dgketchum
WORKDIR /home/dgketchum
RUN ./install_kit


ENV TERM xterm-256color
CMD sh -c "echo You need to run grass from a shell created by a docker exec invocation; sleep 5d"
