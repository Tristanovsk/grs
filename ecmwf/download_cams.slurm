#!/bin/bash

#SBATCH --account=cnes_level2
#SBATCH --job-name=ecmwf_cams
#SBATCH --mem-per-cpu=8G
#SBATCH -N 1 # number of nodes
#SBATCH -n 8 # number of cores
#SBATCH --time=08:00:00

# Exemple de Job array permettant de récupérer des données pour un mois entier
module load conda

export https_proxy=http://<proxy-login>:<proxy-password>@proxy-tech-web.cnes.fr:8060
export http_proxy=$https_proxy
export REQUESTS_CA_BUNDLE="/etc/pki/tls/certs/ca-bundle.crt"
export SSL_CERT_FILE="/etc/pki/tls/certs/ca-bundle.crt"

conda activate grs2_3
# using gitlab package
api_call=/work/SWOT/swotdev/aval/OBS2CO/git/grs2/ecmwf/download_cams.py
# api_call=/work/scratch/env/zeghoudij/workspace/grs2/ecmwf/download_cams_day.py
# by default download from ads at J-2 with relative humidity from ECMWF at J-2
$CONDA_PREFIX/bin/python $api_call

# SAMPLE USAGE:

# to get data from ads for the previous month of the current year
#$CONDA_PREFIX/bin/python $api_call --outdir ~/scratch_data/cams_perso --source copernicus --month $(($(date +%m)-1))

# to get data for a custom date range and combine the result into one netcdf file
# CONDA_PREFIX/bin/python $api_call --cstart 2023-10-20 --cend 2023-10-30 --estart 2023-10-20 --eend 2023-10-30 --combine

# to get data from eac4 reanalysis for november 2022
# $CONDA_PREFIX/bin/python $api_call --mode reanalysis --source copernicus --month 11 --year 2022

# to get the relative humidity from ecmwf for November of the current year
# $CONDA_PREFIX/bin/python $api_call --source ecmwf --month 11
