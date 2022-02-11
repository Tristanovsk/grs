ARG IMAGE_SOURCE

FROM ${IMAGE_SOURCE}/snap

USER root
LABEL maintainer="OBS2CO"


# Montage du secret contenant un password pour se connecter au proxy du cnes
# Il faut utiliser le secret dans le même run que le montage sinon cela ne fonctionnera pas. Car les secrets sont montes seulement dans une commande
RUN --mount=type=secret,id=proxy_http_cnes \ 
    export http_proxy=$(cat /run/secrets/proxy_http_cnes) && export https_proxy=$(cat /run/secrets/proxy_https_cnes) && \
    apt-get -y update && \
    apt-get -y install ca-certificates gfortran


#Ajout des certificats
COPY certs/* /usr/local/share/ca-certificates/
RUN update-ca-certificates

# UL : installation Conda apres la mise a jour des certificats pour atteindre Artifactory. 
RUN --mount=type=secret,id=arti_conda_repo \
    CONDA_SSL_VERIFY=/etc/ssl/certs/ca-certificates.crt conda install --override-channels -c $(cat /run/secrets/arti_conda_repo) gdal

COPY . /home/jovyan/grs

RUN --mount=type=secret,id=arti_pip_repo \
    PIP_CERT=/etc/ssl/certs/ca-certificates.crt pip install -i $(cat /run/secrets/arti_pip_repo) -r /home/jovyan/grs/requirements.txt
    
WORKDIR /home/jovyan/grs
RUN ls
RUN make clean
RUN make
RUN python /home/jovyan/grs2/setup.py build && python /home/jovyan/grs2/setup.py install

RUN ln -s /srv/conda/envs/env_snap/lib/python3.9/site-packages/snappy /srv/conda/envs/env_snap/lib/python3.9/site-packages/esasnappy

CMD grs
