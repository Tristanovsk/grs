ARG IMAGE_SOURCE

FROM ${IMAGE_SOURCE}/ubuntu

USER root
LABEL maintainer="OBS2CO"

ENV LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    LANGUAGE=en_US.UTF-8 \
    SHELL=/bin/bash \
    PATH=/srv/conda/bin:$PATH \
    DEBIAN_FRONTEND=noninteractive \
    NB_USER=jovyan \
    NB_UID=1002 \
    NB_GID=1002 \
    APP_BASE=/home/${NB_USER}/
    
ENV USER=${NB_USER} \ 
    HOME=/home/${NB_USER} \
    CONDA_DIR=${APP_BASE}/conda
    
RUN groupadd --gid ${NB_GID} ${NB_USER}                                                                                             && \
    useradd --comment "Default user" --create-home --gid ${NB_GID} --no-log-init --shell /bin/bash --uid ${NB_UID} ${NB_USER}
    
RUN --mount=type=secret,id=proxy_http_cnes \ 
    export http_proxy=$(cat /run/secrets/proxy_http_cnes) && export https_proxy=$(cat /run/secrets/proxy_http_cnes) && \
    apt-get -qq update                                                                                                              && \
    apt-get -qq install --yes apt-utils                                                                                             && \
    apt-get -qq install --yes --no-install-recommends ttf-dejavu \
            wget make g++ sudo vim less unzip tree file libgfortran5 locales > /dev/null                                            && \
    apt-get -qq purge                                                                                                               && \
    apt-get -qq clean                                                                                                               && \
    rm -rf /var/lib/apt/lists/*                                                                                                     && \
    echo "LC_ALL=en_US.UTF-8" >> /etc/environment                                                                                   && \
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen                                                                                     && \
    echo "LANG=en_US.UTF-8" > /etc/locale.conf                                                                                      && \
    locale-gen en_US.UTF-8

# conda installation via miniforge
ADD install-miniforge.bash /tmp/install-miniforge.bash

RUN chmod 755 /tmp/install-miniforge.bash                                                                                           && \
    /tmp/install-miniforge.bash                                                                                                     && \
    rm -f /tmp/install-*.bash                                                                                                       && \
    chown -R $NB_USER:$NB_GID ${HOME} 

USER ${NB_USER}

ENV PATH=${CONDA_DIR}/envs/env_snap/snap/bin:${CONDA_DIR}/envs/env_snap/bin:${CONDA_DIR}/bin:$PATH \
    KERNEL_PYTHON_PREFIX=${CONDA_DIR}/envs/env_snap  \
    PREFIX=${CONDA_DIR}/envs/env_snap

ADD environment.yml /tmp/environment.yml

RUN mamba env create -f /tmp/environment.yml                                                                 && \
    mamba clean --all -f -y

WORKDIR ${HOME}

# Montage du secret contenant un password pour se connecter au proxy du cnes
## Il faut utiliser le secret dans le même run que le montage sinon cela ne fonctionnera pas. Car les secrets sont montes seulement dans une commande
RUN --mount=type=secret,id=proxy_http_cnes \ 
    export http_proxy=$(cat /run/secrets/proxy_http_cnes) && export https_proxy=$(cat /run/secrets/proxy_http_cnes) && \
    apt-get -y install ca-certificates gfortran

#Ajout des certificats
COPY certs/* /usr/local/share/ca-certificates/
RUN update-ca-certificates

# UL : installation Conda apres la mise a jour des certificats pour atteindre Artifactory. 
#RUN --mount=type=secret,id=arti_conda_repo \
#    CONDA_SSL_VERIFY=/etc/ssl/certs/ca-certificates.crt conda install --override-channels -c $(cat /run/secrets/arti_conda_repo) gdal

COPY . ${HOME}/grs
RUN chmod -R 777 ${HOME}


RUN --mount=type=secret,id=arti_pip_repo \
    PIP_CERT=/etc/ssl/certs/ca-certificates.crt pip install -i $(cat /run/secrets/arti_pip_repo) -r /home/jovyan/grs/requirements.txt

RUN ln -s /srv/conda/envs/env_snap/lib/python3.9/site-packages/snappy /srv/conda/envs/env_snap/lib/python3.9/site-packages/esasnappy

WORKDIR ${HOME}/grs
RUN make clean && make

RUN python setup.py build 
RUN python setup.py install

#RUN sed -i -e '/default\_userdir= =/ s/= .*/= \/home\/jovyan\/.snap/' /srv/conda/envs/env_snap/snap/etc/snap.conf
RUN echo 'snap.versionCheck.interval=NEVER\nsnap.jai.tileCacheSize=1024' > /srv/conda/envs/env_snap/snap/.snap/etc/snap.properties

RUN chmod -R 777 ${HOME}

#RUN grs -h

#CMD grs
