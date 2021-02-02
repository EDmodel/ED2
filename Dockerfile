# ----------------------------------------------------------------------
# BUILD MODEL BINARY
# ----------------------------------------------------------------------
FROM debian:buster-slim AS builder

# Some variables that can be used to set control the docker build
ARG MODEL_VERSION="2.2.0"
ARG BINARY_VERSION="2.2"

# specify fortran compiler
ENV FC_TYPE=GNU

# install dependencies
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential \
       curl \
       gfortran \
       git \
       libhdf5-dev \
       libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

# download, unzip and build ed2
COPY . /src/
WORKDIR /src/ED/build
RUN ./install.sh -g -p docker && mv ed_*-opt ed

########################################################################

# ----------------------------------------------------------------------
# BUILD PECAN FOR MODEL
# ----------------------------------------------------------------------
FROM debian:buster-slim

# ----------------------------------------------------------------------
# INSTALL MODEL SPECIFIC PIECES
# ----------------------------------------------------------------------

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       libhdf5-103 \
       libopenmpi3 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /data

# Some variables that can be used to set control the docker build
ARG MODEL_VERSION="2.2.0"
ENV MODEL_VERSION="${MODEL_VERSION}"

# COPY model binary
COPY --from=builder /src/ED/build/ed /usr/local/bin/ed


