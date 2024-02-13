ARG IMAGE_SOURCE
FROM ${IMAGE_SOURCE}ubuntu:22.04 AS systemdependencies
# FROM ubuntu:kinetic AS systemdependencies
LABEL maintainer: "robin.buratti@magellium.fr"

ENV LANG C.UTF-8
ENV LC_ C.UTF-8

RUN ulimit -s unlimited

# Proxy from secret volumes
RUN if [ -f "/kaniko/run/secrets/http_proxy" ]; then export http_proxy=$(cat /kaniko/run/secrets/http_proxy); export https_proxy=$(cat /kaniko/run/secrets/https_proxy); fi && \
    apt-get update -y && \
    apt-get install -y ca-certificates

# Ajout des certificats
COPY cert[s]/* /usr/local/share/ca-certificates/
RUN update-ca-certificates

# Install libraries
RUN if [ -f "/kaniko/run/secrets/http_proxy" ]; then export http_proxy=$(cat /kaniko/run/secrets/http_proxy); export https_proxy=$(cat /kaniko/run/secrets/https_proxy); fi \
    && apt-get -qq update \
    && DEBIAN_FRONTEND=noninteractive apt-get -qq install -y --no-install-recommends \
        software-properties-common \
        gcc \
        gfortran \
        python3.10 \
        python3-dev \
        build-essential \
        gdal-bin \
        libgdal-dev \
    && rm -rf /var/lib/apt/lists/*
  
ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
ENV C_INCLUDE_PATH=/usr/include/gdal

# GRS INSTALL
WORKDIR /home/
COPY ecmwf ./grs2/ecmwf
COPY exe ./grs2/exe
COPY grs ./grs2/grs
COPY grsdata ./grs2/grsdata
COPY Makefile ./grs2/
COPY setup.py ./grs2/
WORKDIR /home/grs2

RUN make

#########################
FROM ${IMAGE_SOURCE}ubuntu:22.04
LABEL maintainer: "robin.buratti@magellium.fr"

ENV LANG C.UTF-8
ENV LC_ C.UTF-8

RUN ulimit -s unlimited

RUN if [ -f "/kaniko/run/secrets/http_proxy" ]; then export http_proxy=$(cat /kaniko/run/secrets/http_proxy); export https_proxy=$(cat /kaniko/run/secrets/https_proxy); fi \
    && apt-get -qq update \
    && DEBIAN_FRONTEND=noninteractive apt-get -qq install -y --no-install-recommends \
        python-is-python3 \
        python3.10 \
        python3-dev \
        python3-pip \
        python3-affine \
        python3-gdal \
        python3-lxml \
        python3-xmltodict \
        gdal-bin \
    && rm -rf /var/lib/apt/lists/*
  
# get GRS from systemdependencies
COPY --from=systemdependencies /home/grs2 /home/grs2

WORKDIR /home/grs2

# Add additionnal dependencies + GRS
RUN if [ -f "/kaniko/run/secrets/http_proxy" ]; then export http_proxy=$(cat /kaniko/run/secrets/http_proxy); export https_proxy=$(cat /kaniko/run/secrets/https_proxy); fi \
    && pip3 install \
        --trusted-host pypi.org --trusted-host pypi.python.org --trusted-host files.pythonhosted.org \
        --no-cache-dir \
        # requirements
        cdsapi \
        dask \
        dask[array] \
        docopt \
        geopandas \
        matplotlib \
        netCDF4 \
        numpy \
        pandas \
        pyproj \
        python-dateutil \
        scipy \
        "xarray<=2023.4.2" \
        # other dependencies        
        "eoreader<=0.19.4" \
         numba \
    && pip3 install --trusted-host pypi.org --trusted-host pypi.python.org --trusted-host files.pythonhosted.org .

RUN mkdir -p /datalake/watcal/GRS \
    && cp -r grsdata /datalake/watcal/GRS/

WORKDIR /home/
#ENTRYPOINT ["tail", "-f", "/dev/null"]
