ARG IMAGE_SOURCE
FROM ${IMAGE_SOURCE}ubuntu:22.04 AS systemdependencies
LABEL maintainer="robin.buratti@magellium.fr"

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

RUN ulimit -s unlimited

# Proxy depuis les secrets Podman + installation ca-certificates
RUN --mount=type=secret,id=http_proxy --mount=type=secret,id=https_proxy \
    if [ -f "/run/secrets/http_proxy" ]; then \
        export http_proxy=$(cat /run/secrets/http_proxy); \
        export https_proxy=$(cat /run/secrets/https_proxy); \
    fi && \
    apt-get update -y && \
    apt-get install -y ca-certificates

# Ajout des certificats CNES / Entreprise
COPY cert[s]/* /usr/local/share/ca-certificates/
RUN update-ca-certificates
COPY certs/ca-bundle.crt /usr/local/share/ca-certificates/ca-bundle.crt

# Install libraries système
RUN --mount=type=secret,id=http_proxy --mount=type=secret,id=https_proxy \
    if [ -f "/run/secrets/http_proxy" ]; then \
        export http_proxy=$(cat /run/secrets/http_proxy); \
        export https_proxy=$(cat /run/secrets/https_proxy); \
    fi \
    && apt-get -qq update \
    && DEBIAN_FRONTEND=noninteractive apt-get -qq install -y --no-install-recommends \
        software-properties-common \
        gcc \
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
COPY pyproject.toml ./grs2/
COPY requirements.txt ./grs2/
WORKDIR /home/grs2


#########################
# STAGE 2: Image finale
#########################
FROM ${IMAGE_SOURCE}ubuntu:22.04
LABEL maintainer="robin.buratti@magellium.fr"

ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8

RUN ulimit -s unlimited

# Récupération et activation immédiate des certificats du stage 1 dans le stage 2
COPY --from=systemdependencies /usr/local/share/ca-certificates/ /usr/local/share/ca-certificates/
COPY --from=systemdependencies /etc/ssl/certs/ca-certificates.crt /etc/ssl/certs/ca-certificates.crt

RUN --mount=type=secret,id=http_proxy --mount=type=secret,id=https_proxy \
    if [ -f "/run/secrets/http_proxy" ]; then \
        export http_proxy=$(cat /run/secrets/http_proxy); \
        export https_proxy=$(cat /run/secrets/https_proxy); \
    fi \
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
  
# Récupération de GRS depuis systemdependencies
COPY --from=systemdependencies /home/grs2 /home/grs2

WORKDIR /home/grs2

ARG ARTIFACTORY_HOST="artifactory.cnes.fr"

ENV PIP_CERT=/etc/ssl/certs/ca-certificates.crt

# Installation Python robuste (Les certificats système étant là, pip fera confiance au proxy !)
RUN --mount=type=secret,id=http_proxy \
    --mount=type=secret,id=https_proxy \
    --mount=type=secret,id=artifactory_url \
    if [ -f /run/secrets/http_proxy ]; then \
        export http_proxy=$(cat /run/secrets/http_proxy); \
        export https_proxy=$(cat /run/secrets/https_proxy); \
    fi \
    && pip3 install \
        --timeout 3000 \
        --upgrade pip \
    && ARTIFACTORY_SECRET_URL=$(cat /run/secrets/artifactory_url) \
    && pip3 install \
        --timeout 3000 \
        --retries 10 \
        --progress-bar off \
        --no-cache-dir \
        --index-url "$ARTIFACTORY_SECRET_URL" \
        --trusted-host ${ARTIFACTORY_HOST} \
        -r requirements.txt \
    && pip3 install \
        --timeout 3000 .

RUN mkdir -p /datalake/watcal/GRS \
    && cp -r grsdata /datalake/watcal/GRS/

WORKDIR /home/