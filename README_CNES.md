# GRS Installation and Usage on CNES Machines

This document gathers the instructions specific to installing and running GRS on the CNES TREX HPC
cluster (SLURM). For general installation instructions (any other machine), see the main
[README](README.md) and the [documentation](https://cnes.github.io/GRSprocessor).

## Installation on TREX

1. Clone the repository from the internal GitLab mirror (https or ssh) and choose the branch or tag to
   install:

   ```bash
   git clone https://gitlab.cnes.fr/waterquality/grs2.git
   # or: git clone git@gitlab.cnes.fr:waterquality/grs2.git
   cd grs2
   git checkout <branch-or-tag>
   ```

2. Create a conda environment (Python 3.11 or newer; conda only provides Python, the packages are
   installed with pip) and install GRS. The TREX nodes have no GPU: install the CPU build of PyTorch
   (used by the OmniCloudMask cloud masks) first, much smaller than the default CUDA one:

   ```bash
   ml conda
   conda create -n grs python=3.12 pip
   conda activate grs
   python -m pip install torch --index-url https://download.pytorch.org/whl/cpu
   python -m pip install .
   ```

3. Set the folder of the LUTs, either in the user configuration file (once):

   ```bash
   grs_config --init /work/datalake/watcal/GRS/grsdata_v21
   ```

   or in the environment of the jobs:

   ```bash
   export GRS_DATA_ROOT=/work/datalake/watcal/GRS/grsdata_v21
   ```

   The compute nodes have no internet access: switch the automatic download of the LUTs off with
   `export GRS_NO_DOWNLOAD=1` (or `zenodo: auto_download: false` in the user configuration). If a LUT
   is missing from the folder, download it from a node with internet access with
   `grs_config --download`.

4. Check the installation and the configuration (`grs_config` lists the LUTs and exits with code 1 if
   one is missing):

   ```bash
   grs -v
   grs_config
   ```

## CAMS data on TREX

The default configuration reads the CAMS forecasts from the TREX archive `/work/datalake/watcal/cams`
(`path: cams`, daily files `YYYY/MM/DD/YYYY-MM-DD-cams-global-atmospheric-composition-forecasts.nc`):
on TREX, `--cams_file` can be omitted. Give `--cams_file` to use another file. The scripts of the
[ecmwf](ecmwf) folder (e.g. `download_cams.slurm`) download the CAMS data.

## Running GRS on TREX with a SLURM interactive job

```bash
unset SLURM_JOB_ID
srun -A cnes_level2 -N 1 -c 8 --time=02:00:00 --mem=64G --pty bash
ml conda
conda activate grs
grs /work/datalake/S2-L1C/31TFJ/2023/06/16/S2B_MSIL1C_20230616T103629_N0509_R008_T31TFJ_20230616T111826.SAFE \
    --odir /work/datalake/watcal/test --resolution 20 \
    --dem_file /work/datalake/static_aux/MNT/COP-DEM_GLO-30-DGED_S2_tiles/COP-DEM_GLO-30-DGED_31TFJ.tif
```

GRS uses `processor: ncpu` threads (8 by default): request as many cores (`-c`) as threads. A Sentinel-2
tile at 20 m needs about 10 GB of memory. See the [Usage](README.md#usage) section of the main README
for the `grs` options (`grs -h` lists them all).

For batch jobs, see the SLURM scripts under [grs/tests/ci](grs/tests/ci) and [ecmwf](ecmwf), and the
["Scheduling and Triggers"](https://cnes.github.io/GRSprocessor/processing_chain.html#scheduling-and-triggers)
section of the documentation.

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
