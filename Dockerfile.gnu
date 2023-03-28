# ----------------------------------------------------------------------
# Build ED2 model
# ----------------------------------------------------------------------
#FROM debian:bullseye-slim AS build
FROM ubuntu:20.04 AS build

# Some variables that can be used to set control the docker build
ARG ED2_KIND=E

# environment variables controlling the build
ENV FC_TYPE=GNU \
    ED2_KIND=${ED2_KIND} \
    TZ=America/Chicago \
    DEBIAN_FRONTEND=noninteractive

# install dependencies
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential \
       gfortran \
       libhdf5-openmpi-dev \
    && rm -rf /var/lib/apt/lists/*

# copy the source only, this prevents a full rebuild in case of changes to Dockerfile or non source files
COPY BRAMS    /ED2/BRAMS
COPY ED       /ED2/ED
COPY EDR      /ED2/EDR
COPY RAPP     /ED2/RAPP
COPY Ramspost /ED2/Ramspost

WORKDIR /ED2/ED/build
RUN ./install.sh -k ${ED2_KIND} -g -p docker.gnu
RUN if [ -e ed_*-opt ]; then mv ed_*-opt ed2; else mv ed_*-dbg ed2; fi

########################################################################

# ----------------------------------------------------------------------
# Minimal image with just ED2 model
# ----------------------------------------------------------------------
#FROM debian:bullseye-slim
FROM ubuntu:20.04

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       libhdf5-openmpi-103 \
       libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /data
COPY --from=build /ED2/ED/build/ed2 /usr/bin

CMD /usr/bin/ed2
