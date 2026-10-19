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

# Python 3.12 from conda (grs requires Python >= 3.12); all the other dependencies,
# including numpy >= 2, are installed with pip from requirements.txt (GDAL is not needed).
# An existing environment with another Python version is recreated.
if conda env list 2>/dev/null | grep -q "${ENV_PATH}" && \
   ! conda run -p "${ENV_PATH}" python -c "import sys; assert sys.version_info[:2] == (3, 12)" 2>/dev/null; then
    conda env remove -y -p "${ENV_PATH}"
fi
if ! conda env list 2>/dev/null | grep -q "${ENV_PATH}"; then
    conda create -y -p "${ENV_PATH}" \
        --override-channels \
        -c "${CONDA_BASE}/conda-forge-remote" \
        -c "${CONDA_BASE}/conda-main-remote" \
        python=3.12
fi

conda activate "${ENV_PATH}"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"
pip install -U pip
pip install 'setuptools<80'
pip install -r requirements.txt
pip install .[dev]
pip install coverage

mkdir -p /work/scratch/data/${USER}/OBS2CO/CI/tests/coverage/
