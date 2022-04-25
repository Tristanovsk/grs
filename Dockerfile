ARG IMAGE_SOURCE

FROM ${IMAGE_SOURCE}/snap

USER root
LABEL maintainer="OBS2CO"

RUN mkdir /app
RUN mkdir /app/.snap

# Montage du secret contenant un password pour se connecter au proxy du cnes
## Il faut utiliser le secret dans le même run que le montage sinon cela ne fonctionnera pas. Car les secrets sont montes seulement dans une commande
RUN --mount=type=secret,id=proxy_http_cnes \ 
    export http_proxy=$(cat /run/secrets/proxy_http_cnes) && export https_proxy=$(cat /run/secrets/proxy_http_cnes) && \
    apt-get -y update && \
    apt-get -y install ca-certificates gfortran

#Ajout des certificats
COPY certs/* /usr/local/share/ca-certificates/
RUN update-ca-certificates

# UL : installation Conda apres la mise a jour des certificats pour atteindre Artifactory. 
#RUN --mount=type=secret,id=arti_conda_repo \
#    CONDA_SSL_VERIFY=/etc/ssl/certs/ca-certificates.crt conda install --override-channels -c $(cat /run/secrets/arti_conda_repo) gdal

RUN chmod -R 777 /app
COPY . /app/grs
WORKDIR /app/grs 

RUN ln -s /srv/conda/envs/env_snap/snap/.snap/snap-python/snappy /srv/conda/envs/env_snap/lib/python3.9/site-packages/esasnappy

RUN --mount=type=secret,id=arti_pip_repo \
    PIP_CERT=/etc/ssl/certs/ca-certificates.crt pip install -i $(cat /run/secrets/arti_pip_repo) -r /app/grs/requirements.txt

RUN make clean && make

RUN python setup.py build 
RUN python setup.py install

RUN echo 'snap.versionCheck.interval=NEVER\nsnap.jai.tileCacheSize=1024' > /srv/conda/envs/env_snap/snap/.snap/etc/snap.properties

RUN sed -i 's#/srv/conda/envs/env_snap/snap//.snap/system#/app/.snap/system/#g' /srv/conda/envs/env_snap/snap/etc/snap.conf
RUN sed -i 's#/srv/conda/envs/env_snap/snap/.snap#//app/.snap/#g' /srv/conda/envs/env_snap/snap//etc/snap.properties
RUN echo 'snap.versionCheck.interval=NEVER\nsnap.jai.tileCacheSize=1024' >> /srv/conda/envs/env_snap/snap/etc/snap.properties
RUN sed -i '11 a AuxDataPath = /app/.snap/auxdata/' /srv/conda/envs/env_snap/snap//etc/snap.auxdata.properties

RUN --mount=type=secret,id=proxy_http_cnes \ 
    export http_proxy=$(cat /run/secrets/proxy_http_cnes) && export https_proxy=$(cat /run/secrets/proxy_http_cnes) && \
    timeout 300 snap --nosplash --nogui --modules --update-all || true

#RUN cp /app/grs/snap.auxdata.properties /srv/conda/envs/env_snap/snap/etc/snap.auxdata.properties

RUN chmod -R 777 /app
RUN mkdir /snap && chmod -R 777 /snap

#ENTRYPOINT ['python', '/app/grs/exe/launcher.py', "/app/grs/exe/global_config.yml'] 
