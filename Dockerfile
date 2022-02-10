ARG IMAGE_SOURCE

FROM ${IMAGE_SOURCE}/snap


# Montage du volume temporaire et utilisation pour apt le site du cnes
# Il faut utiliser le secret dans le même run que le montage sinon cela ne fonctionnera pas

USER root
ARG HTTP_PROXY
ARG ARTI_CONDA

RUN mkdir -p ~/.pip/
RUN touch ~/.pip/pip.conf
RUN echo '[global]\nindex-url = https://${ARTI_CONDA}/api/pypi/pypi/simple'  > ~/.pip/pip.conf

RUN cat ~/.pip/pip.conf
RUN mkdir -p /home/jovyan/grs

COPY * /home/jovyan/grs
DIR /home/jovyan/grs

RUN pip install -r /home/jovyan/grs2/requirements.txt

RUN conda install --override-channels -c ${ARTI_CONDA}/api/conda/conda/ gdal

RUN --mount=type=secret,id=proxy_http_cnes \ 
    --mount=type=secret,id=proxy_https_cnes \
    export http_proxy=$(cat /run/secrets/proxy_http_cnes) && export https_proxy=$(cat /run/secrets/proxy_https_cnes) && \
    apt-get -y update && \
    apt-get -y install ca-certificates && \
    apt-get -y install gfortran

#Ajout des certificats
COPY certs/* /usr/local/share/ca-certificates/
RUN update-ca-certificates

LABEL maintainer="obs2co"

#FROM docker.pkg.github.com/snap-contrib/docker-snap/snap
    
# && \ pip install -r /home/jovyan/grs2/requirements.txt

RUN ln -s /srv/conda/envs/env_snap/lib/python3.9/site-packages/snappy /srv/conda/envs/env_snap/lib/python3.9/site-packages/esasnappy

RUN cd /home/jovyan/grs2 && make clean && make
RUN python /home/jovyan/grs2/setup.py build && python /home/jovyan/grs2/setup.py install

CMD grs
