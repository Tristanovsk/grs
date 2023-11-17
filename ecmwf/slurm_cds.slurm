#!/bin/bash

#SBATCH --account=cnes_level2
#SBATCH --job-name=ecmwf_cams
#SBATCH --mem-per-cpu=8G
#SBATCH -N 1 # number of nodes
#SBATCH -n 4 # number of cores
#SBATCH --time=48:00

# Recuperer les donnees ecmwf cams à travers un job PBS  :
# To be launched with -- export option

module load conda

export https_proxy=http://<proxy-login>:<proxy-password>@proxy-tech-web.cnes.fr:8060
export http_proxy=$https_proxy
export REQUESTS_CA_BUNDLE="/etc/pki/tls/certs/ca-bundle.crt"
export SSL_CERT_FILE="/etc/pki/tls/certs/ca-bundle.crt"

# conda activate /work/scratch/harmelt/grs_py3.6
conda activate /work/SWOT/swotdev/aval/OBS2CO/conda_env/grs_py3.8
# using gitlab package
api_call=/work/SWOT/swotdev/aval/OBS2CO/git/grs2/ecmwf/download_cams_cds_day.py
#api_call=/work/scratch/env/zeghoudij/workspace/grs2/ecmwf/download_cams_cds_day.py

$CONDA_PREFIX/bin/python $api_call cams_forecast "/work/datalake/watcal/ECMWF/CAMS/" 

