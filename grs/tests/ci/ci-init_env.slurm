#!/bin/bash
#
#SBATCH --job-name=grs-init_env
#SBATCH --time=00:30:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem-per-cpu=4G
#SBATCH --account=campus

#cnes_level2

set -e

cd ${MY_JOB_WORKDIR}/grs2
module load conda
eval "$(conda shell.bash hook)"

ENV_PATH="/work/scratch/data/${USER}/.conda/envs/grs2_tests"
CONDA_BASE="https://${ARTIFACTORY_USER}:${ARTIFACTORY_TOKEN}@artifactory.cnes.fr/artifactory/api/conda"

# GDAL via conda (osgeo bindings need C headers); numpy pinned to 1.25 for numba compat
if ! conda env list 2>/dev/null | grep -q "${ENV_PATH}"; then
    conda create -y -p "${ENV_PATH}" \
        --override-channels \
        -c "${CONDA_BASE}/conda-forge-remote" \
        -c "${CONDA_BASE}/conda-main-remote" \
        python=3.10 'numpy=1.25' gdal
else
    if ! conda run -p "${ENV_PATH}" python -c "
from osgeo import gdal
import numpy as np
assert np.__version__.startswith('1.25'), f'numpy={np.__version__}'
" 2>/dev/null; then
        conda install -y -p "${ENV_PATH}" \
            --override-channels \
            -c "${CONDA_BASE}/conda-forge-remote" \
            -c "${CONDA_BASE}/conda-main-remote" \
            'numpy=1.25' gdal
    fi
fi

conda activate "${ENV_PATH}"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"
pip install -U pip
pip install 'setuptools<80'
pip install -r requirements.txt
pip install .[dev]
pip install coverage

mkdir -p /work/scratch/data/${USER}/OBS2CO/CI/tests/coverage/
