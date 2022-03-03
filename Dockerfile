ARG IMAGE_SOURCE

FROM ${IMAGE_SOURCE}/snap

USER root
LABEL maintainer="OBS2CO"

# Montage du secret contenant un password pour se connecter au proxy du cnes
## Il faut utiliser le secret dans le même run que le montage sinon cela ne fonctionnera pas. Car les secrets sont montes seulement dans une commande
RUN --mount=type=secret,id=proxy_http_cnes \ 
    export http_proxy=$(cat /run/secrets/proxy_http_cnes) && export https_proxy=$(cat /run/secrets/proxy_http_cnes) && \
    apt-get -y update && \
    apt-get -y install ca-certificates gfortran

#Ajout des certificats
COPY certs/* /usr/local/share/ca-certificates/
RUN update-ca-certificates

RUN useradd -ms /bin/bash grsuser
WORKDIR /home/grsuser
RUN usermod -aG sudo grsuser

#USER grsuser

# UL : installation Conda apres la mise a jour des certificats pour atteindre Artifactory. 
RUN --mount=type=secret,id=arti_conda_repo \
    CONDA_SSL_VERIFY=/etc/ssl/certs/ca-certificates.crt conda install --override-channels -c $(cat /run/secrets/arti_conda_repo) gdal

COPY . /home/grsuser/grs
RUN chmod -R 777 /home/grsuser/grs

RUN --mount=type=secret,id=arti_pip_repo \
    PIP_CERT=/etc/ssl/certs/ca-certificates.crt pip install -i $(cat /run/secrets/arti_pip_repo) -r /home/grsuser/grs/requirements.txt

RUN ln -s /srv/conda/envs/env_snap/lib/python3.9/site-packages/snappy /srv/conda/envs/env_snap/lib/python3.9/site-packages/esasnappy

#WORKDIR /home/grsuser/grs
#WORKDIR /home/grsuser/grs/grs/landsat_angles/OLI/
#RUN gcc -g -Wall -O2 -march=nocona -mfpmath=sse -msse2  -I./ias_lib/ -I./ -c -o #l8_angles.o l8_angles.c

WORKDIR /home/grsuser/grs
RUN make clean && make

RUN ls
RUN python setup.py build 

RUN python setup.py install

RUN chmod -R 777 /srv/conda/envs/env_snap/snap/.snap

RUN grs -h

CMD grs
