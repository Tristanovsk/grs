#!/bin/bash

#SBATCH --array=1-30:1
#SBATCH --account=cnes_level2
#SBATCH --job-name=ecmwf_cams
#SBATCH --mem-per-cpu=8G
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of cores
#SBATCH --time=08:00

# Exemple de Job array permettant de récupérer des données pour un mois entier
module load conda

export https_proxy=http://<proxy-login>:<proxy-password>@proxy-tech-web.cnes.fr:8060
export http_proxy=$https_proxy
export REQUESTS_CA_BUNDLE="/etc/pki/tls/certs/ca-bundle.crt"
export SSL_CERT_FILE="/etc/pki/tls/certs/ca-bundle.crt"

# conda activate /work/scratch/harmelt/grs_py3.6
conda activate /work/SWOT/swotdev/aval/OBS2CO/conda_env/grs_py3.8
# using gitlab package
api_call=/work/SWOT/swotdev/aval/OBS2CO/git/grs2/ecmwf/download_cams_day.py

# by default download from ads at J-7 with relative humidity from ECMWF at J-2
$CONDA_PREFIX/bin/python $api_call

# to get data from eac4 reanalysis for November 2023
# $CONDA_PREFIX/bin/python $api_call --mode reanalysis --source copernicus --cstart 2023-11-01 --cend 2023-11-31

# to get the relative humidity from ecmwf for November 2023
# $CONDA_PREFIX/bin/python $api_call --source ecmwf --estart 2023-11-01 --eend 2023-11-31