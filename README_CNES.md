# GRS Installation and Usage on CNES Machines

This document gathers the instructions specific to installing and running GRS on
CNES infrastructure (TREX, the HAL HPC cluster, and the PBS/SLURM clusters). For
general installation instructions (any other machine), see the main
[README](README.md).

## Installation on TREX (CNES)

1. First clone the repository (from https or ssh):
```commandline
git clone https://gitlab.cnes.fr/waterquality/grs2.git
```
or
```commandline
git clone git@gitlab.cnes.fr:waterquality/grs2.git
```

Choose your branch (example grs_cnes_v2.1.6)
```commandline
git checkout grs_cnes_v2.1.6
```

2. Make sure that the grsdata variable is set as follows in the config.yml file:
```commandline
grsdata: '/work/datalake//watcal/GRS/grsdata_v21'
```

3. To complete installation please activate your conda grs_cnes environment as follows:
```commandline
ml conda
conda activate grs_cnes
pip install .
```

You are done, please check [Running GRS on TREX](#running-grs-on-trex-with-a-slurm-interactive-job) below.

## Installation on the HAL CNES HPC

Script for installation on the HAL CNES HPC:
```commandline
# set your grs path here
your_path_to_grs=/work/scratch/$USER/dev/grs

cd $your_path_to_grs
git clone git@gitlab.cnes.fr:waterquality/grs2.git
ml conda/4.12.0
mkdir /work/scratch/$USER/tmp
export TMPDIR=/work/scratch/$USER/tmp
conda create python=3.10 -n grs_v2
conda activate grs_v2
conda install gdal geopandas -c conda-forge-remote
pip install cdsapi netCDF4 matplotlib docopt xarray dask dask[array] toolz>=0.8.2 affine xmltodict bokeh eoreader lxml numba
ml gcc
make
pip install .

grs -h
```

## Installation on the PBS cluster (legacy)

> **Note:** this PBS workflow is kept for reference but is superseded by SLURM on the
> current CNES cluster. See the SLURM scripts under
> [grs/tests/ci](grs/tests/ci) and [ecmwf](ecmwf) (e.g. `ci-func-run.slurm`,
> `download_cams.slurm`) for up-to-date job submission examples, and the
> ["Scheduling and Triggers"](https://cnes.github.io/GRSprocessor/processing_chain.html#scheduling-and-triggers)
> section of the documentation.

Installing from sources with conda on the CNES cluster:

Create the conda environment using the definition file available in the conda folder:
```
conda env create -f conda/grs_conda_3.6.yml -p /work/scratch/$user/grs_py3.6
```
The option -p set the directory where the conda environment will be installed

To install the package grs in conda:

```
source conda/conda_grs.sh -ci
```

To launch GRS on a pbs node:

```
qsub launch_grs_exemple.pbs
```

## Running GRS on TREX with a SLURM interactive job

If you are on TREX CNES you can run the grs example using a SLURM interactive job:
```commandline
unset SLURM_JOB_ID
srun -A cnes_level2 -N 1 -c 8 --time=02:00:00 --mem=64G --x11 --pty bash
ml conda
conda activate grs_cnes
grs /work/datalake/S2-L1C/31TFJ/2023/06/16/S2B_MSIL1C_20230616T103629_N0509_R008_T31TFJ_20230616T111826.SAFE --cams_file /work/datalake/watcal/ECMWF/CAMS/2023/06/16/2023-06-16-cams-global-atmospheric-composition-forecasts.nc --odir /work/datalake/watcal/test --resolution 20 --dem_file /work/datalake/static_aux/MNT/COP-DEM_GLO-30-DGED_S2_tiles/COP-DEM_GLO-30-DGED_31TFJ.tif
```

See the [Testing](README.md#testing) section of the main README for the full list of
`grs` command-line options.

## Running GRS with Docker on CNES [deprecated]

```
qsub -q qdev -I -l walltime=4:00:00

/opt/bin/drunner run -it -v /datalake/watcal:/datalake/watcal artifactory.cnes.fr/obs2co-docker/grs:1.4.0 python /app/grs/exe/launcher.py /app/grs/exe//app/grs/exe/global_config.yml
```

## GitLab CI (internal CNES)

GRSprocessor is developed on GitHub; the `sync` job of the GitHub Actions
[`main.yml`](.github/workflows/main.yml) pipeline strips Git-LFS pointer files from the
whole history and force-pushes the cleaned repository to the internal CNES GitLab mirror
(`gitlab.cnes.fr/waterquality/grs2.git`). This triggers the GitLab pipeline defined in
[`.gitlab-ci.yml`](.gitlab-ci.yml), which runs on the mirrored repository at CNES, on the
`Usine_Logicielle` runners:

1. `init` — `python-init` installs `grs` through the internal Artifactory/JFrog pip
   mirror.
2. `test` — `python-tests` submits SLURM jobs on the `trex.sis.cnes.fr` HPC cluster
   (`ci-init_env.slurm` → `ci-func-run.slurm` → `ci-cleanup_env.slurm`) to run functional
   tests against reference outputs (non-blocking, `allow_failure: true`).
3. `package` — `podman-build` (tags only) builds and pushes the image to the CNES
   Artifactory registry (`obs2co-docker/grs`); `docker-test` (tags only) pulls it back via
   Singularity on the HPC and re-runs the functional tests inside the container;
   `podman-test` is a manual, build-only sanity check.
4. `sonarqube` / `security` — static analysis (SonarQube), SAST, and ClamAV antivirus
   scanning, via shared CNES "Usine Logicielle" pipeline components.

This GitLab mirror pipeline runs the CNES-internal validation (HPC functional tests
against reference data, code quality/security gates) that cannot run outside the CNES
network. See the main [README](README.md#how-the-cicd-pipeline-works) for the public
GitHub Actions side of the pipeline.
