FROM ubuntu:kinetic AS systemdependencies
LABEL maintainer: "robin.buratti@magellium.fr"

ENV LANG C.UTF-8
ENV LC_ C.UTF-8

RUN ulimit -s unlimited

# Install libraries
RUN apt-get update -q
RUN apt-get install -q -y --no-install-recommends apt-utils
RUN apt-get install -q -y software-properties-common
RUN apt-get -q update
RUN apt install -y  curl dash g++  gcc  gettext-base  python-apt-common \
                    python-babel-localedata  python3-apt python3-attr python3-babel  python3-blinker \
                    python3-certifi python3-cffi-backend python3-chardet python3-commandnotfound python3-configobj\
                    python3-cryptography python3-dbus python3-debconf python3-distro python3-httplib2 \
                    python3-idna python3-importlib-metadata python3-jeepney python3-jinja2 python3-json-pointer \
                    python3-jsonpatch python3-jsonschema python3-keyring python3-launchpadlib python3-lazr.restfulclient \
                    python3-lazr.uri python3-markupsafe python3-more-itertools python3-netifaces \
                    python3-pkg-resources  python3-pyparsing python3-pyrsistent python3-requests  \
                    python3-secretstorage python3-serial python3-setuptools python3-six python3-software-properties \
                    python3-tz python3-urllib3 python3-wadllib python3-yaml python3-zipp python3.10-minimal  \
                    python3.10
RUN apt-get -q -y install build-essential
RUN rm -rf /var/lib/apt/lists/*

#CONDA
#remove defaults channel from conda
RUN curl -sL "https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh" > "Miniconda3.sh"
RUN bash Miniconda3.sh -b -p /miniconda3
ENV PATH="/miniconda3/bin/":$PATH
RUN conda update -y conda
RUN rm Miniconda3.sh
RUN echo "conda activate" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]
RUN conda config --remove channels defaults; conda config --append channels conda-forge/label/dev

# CONDA INSTALLS
RUN conda install -y -c conda-forge affine cartopy descartes eoreader gdal geopandas lxml matplotlib numba numpy pandas pyproj rasterio scipy setuptools shapely xarray xmltodict

# PIP INSTALLS
SHELL ["/bin/sh","-c"]
RUN apt-get update -q && apt-get install -y --no-install-recommends apt-utils \
    && apt-get install -q -y software-properties-common \
    && apt-get -q update \
    && apt-get -q -y install build-essential tar file apt-utils pkg-config cmake unzip git wget curl sqlite3 libsqlite3-dev software-properties-common unzip \
    && apt-get -q clean

RUN python3 -m pip install --upgrade pip
RUN apt-get update \
    && apt-get install gdal-bin --yes \
    && apt-get install libgdal-dev --yes
ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
ENV C_INCLUDE_PATH=/usr/include/gdal
RUN python3 -m pip install GDAL

RUN python3 -m pip install  asciitree==0.3.3 azure-core==1.26.4 azure-storage-blob==12.16.0 cdsapi==0.6.1 \
                            click-plugins==1.1.1 colorlog==6.7.0 descartes==1.1.0 docopt==0.6.2 h5py==3.7.0 isodate==0.6.1 jmespath==1.0.1 \
                            munch==2.5.0 netCDF4==1.6.2 pandas==1.5.2 ply==3.11 protobuf==4.21.12 pyarrow==11.0.0 pyasn1==0.4.8 pyasn1-modules==0.2.7 \
                            PyQt5==5.15.7 PyQt5-sip==12.11.0 s3transfer==0.6.0 scipy==1.10.1 snuggs==1.4.7 azure-core==1.26.4 azure-storage-blob==12.16.0 boto3==1.26.122 \
                            botocore==1.29.122 certifi==2022.12.7 cftime==1.6.2 h5py==3.7.0 isodate==0.6.1 rasterio==1.3.6 \
                            setuptools==67.7.2

RUN apt-get install -y -q libstdc++6
RUN ln -sf /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /miniconda3/bin/../lib/libstdc++.so.6

# S2DRIVER INSTALL
WORKDIR /home/
ADD s2driver /home/s2driver
WORKDIR /home/s2driver
RUN python3 setup.py install

# GRS INSTALL
WORKDIR /home/
ADD grs2 /home/grs2
WORKDIR /home/grs2
RUN apt-get update -q
RUN apt-get install -q -y gfortran
RUN make
RUN python3 setup.py install
RUN cp -r grsdata /miniconda3/lib/python3.10/site-packages/grs-2.0.1-py3.10.egg/

WORKDIR /home/
ENTRYPOINT ["tail", "-f", "/dev/null"]
