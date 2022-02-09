ARG IMAGE_SOURCE
ARG HTTP_PROXY

FROM ${IMAGE_SOURCE}/snap


# Montage du volume temporaire et utilisation pour apt le site du cnes
# Il faut utiliser le secret dans le même run que le montage sinon cela ne fonctionnera pas
#RUN --mount=type=secret,id=proxy_http_cnes \ 
#    --mount=type=secret,id=proxy_https_cnes \
#    export http_proxy=$(cat /run/secrets/proxy_http_cnes) && export https_proxy=$(cat /run/secrets/proxy_https_cnes) && \
#   

RUN export http_proxy=${HTTP_PROXY}
RUN export https_proxy=${HTTP_PROXY}
RUN apt-get update && \
    apt install ca-certificates 

#Ajout des certificats
COPY certs/* /usr/local/share/ca-certificates/
RUN update-ca-certificates

LABEL maintainer="obs2co"

#FROM docker.pkg.github.com/snap-contrib/docker-snap/snap

COPY grs /home/jovyan/grs2

RUN cd /home/jovyan/grs2

RUN apt-get update && apt-get -y install ca-certificates
COPY certs/* /usr/local/share/ca-certificates/
RUN update-ca-certificates

RUN conda install gdal
RUN pip install -r /home/jovyan/grs2/requirements.txt

RUN apt-get -y install gfortran

RUN ln -s /srv/conda/envs/env_snap/lib/python3.9/site-packages/snappy /srv/conda/envs/env_snap/lib/python3.9/site-packages/esasnappy

RUN cd /home/jovyan/grs2 && make clean && make
RUN python /home/jovyan/grs2/setup.py build && python /home/jovyan/grs2/setup.py install

CMD grs
