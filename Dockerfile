FROM ubuntu:kinetic AS systemdependencies
LABEL maintainer: "robin.buratti@magellium.fr"

ENV LANG C.UTF-8
ENV LC_ C.UTF-8

RUN ulimit -s unlimited

# Install libraries
RUN apt-get -qq update \
    && apt-get -qq install -y --no-install-recommends software-properties-common curl g++ gcc gfortran \
	                                                  python3.10 build-essential cmake sqlite3 libsqlite3-dev \
    && rm -rf /var/lib/apt/lists/*

#MICROMAMBA
WORKDIR /home
RUN curl -sL micro.mamba.pm/install.sh | bash
ENV MAMBA_EXE="/root/.local/bin/micromamba"
ENV MAMBA_ROOT_PREFIX="/root/micromamba"
ENV PATH="/root/micromamba/bin:$PATH"
ENV PATH="/root/.local/bin/:$PATH"
RUN micromamba update --yes --name base --channel conda-forge micromamba

# MICROMAMBA INSTALLS
RUN micromamba install --yes --name base --channel conda-forge affine cartopy descartes eoreader gdal geopandas lxml matplotlib numba numpy pandas \
                                                               pyproj rasterio scipy setuptools shapely xarray xmltodict \
	&& micromamba clean --all --yes

# GDAL INSTALL
SHELL ["/bin/sh","-c"]
RUN python3 -m pip install --upgrade pip
RUN apt-get update \
  && apt-get -y -q install gdal-bin libgdal-dev
ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
ENV C_INCLUDE_PATH=/usr/include/gdal
RUN python3 -m pip install GDAL

# S2DRIVER INSTALL
WORKDIR /home/
ADD s2driver /home/s2driver
WORKDIR /home/s2driver
RUN python3 setup.py install

# GRS INSTALL
WORKDIR /home/
ADD grs2 /home/grs2
WORKDIR /home/grs2
RUN make \
	&& python3 setup.py install \
	&& ls -sf grsdata /root/micromamba/lib/python3.10/site-packages/grs-2.0.1-py3.10.egg/

WORKDIR /home/
ENTRYPOINT ["tail", "-f", "/dev/null"]

