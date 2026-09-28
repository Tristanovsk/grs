# GRS algorithm package
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue?logo=python&logoColor=white)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/license/Apache-2.0)
[![GRS](https://github.com/CNES/GRSProcessor/actions/workflows/main.yml/badge.svg)](https://github.com/CNES/GRSprocessor/actions)
[![Coverage](badges/coverage.svg)]()
[![Docker Image Version](https://img.shields.io/docker/v/guillaumeeb/grs/2.1.9)](https://hub.docker.com/r/guillaumeeb/grs)
[![PyPI](https://img.shields.io/pypi/v/GRSProcessor)](https://pypi.org/project/GRSprocessor/)

## GRS (Glint Removal for Sentinel-2-like sensors)

Please check [grs documentation](https://cnes.github.io/GRSprocessor/)

The GRS (Glint Removal for Sentinel-2) algorithm [Harmel et al., 2018](https://www.sciencedirect.com/science/article/pii/S0034425717304856)
was specifically developed to
handle and correct for the direct sunlight reflected by the water surface and potentially reaching the sensor (i.e.,
sunglint signal) of Sentinel-2-like mission, that is nadir or near-nadir viewing sensor with SWIR bands. The GRS
processor consists of three main modules to correct for (i) gaseous absorption, (ii) diffuse light from sky and its
reflection by the air-water interface and (iii) the sunglint signal in order to retrieve the water-leaving signal at the
water surface level. 

First, the gaseous absorption (mainly CO2, H2O and O3) correction is performed based on parameterizations of the gas transmittances from full radiative transfer
computations using lidRadtran v2.0.4. Atmospheric pressure and gas concentrations are retrieved from bilinear
interpolation within the grid of the Copernicus Atmosphere Monitoring Service dataset (CAMS). Then, spectral radiances
are corrected for the diffuse sky light and its reflection on the air-water interface. For each pixel, the diffuse
radiance component is reconstructed for the given viewing geometry (i.e., sensor and Sun viewing angles and relative
azimuth) from pre-computed look-up tables (LUT). The Rayleigh optical thickness is rescaled based on the actual pressure
at the scene level to take into account the effects of the altitude on the scattering properties of the atmosphere.
Those LUTs were generated based on the radiative transfer model OSOAA (Chami et al., 2015) for a typical fine and coarse
mode aerosol models, encompassing weakly absorbing ones (Levy et al., 2009), and including the specific spectral
response of the sensor bands. The atmosphere plus surface diffuse signal $`L_{sky}`$ is obtained considering a bimodal aerosol
model (Wang & Gordon, 1994) as follows:

<img src="https://latex.codecogs.com/gif.latex?L_{sky}\left( {\lambda ,{\tau _a}} \right)
 = \gamma L_{sky}^{fine}\left( {\lambda ,{\tau _a}} \right) + \left( {1 - \gamma } \right)L_{sky}^{coarse}\left( {\lambda ,{\tau _a}} \right)"/>

where $`L_{sky}^{fine}`$ and $`L_{sky}^{coarse}`$are the radiances for the fine and coarse aerosol modes, respectively, 
for the aerosol optical thickness $`\tau _a`$; $`\gamma`$ is
the mixing coefficient corresponding to the relative amount of each mode in the atmosphere. Note that $`\tau _a`$ is obtained from
the CAMS dataset (Benedetti et al., 2008; Morcrette et al., 2009) and $`\gamma`$ is retrieved from non-linear fitting including the
LUT aerosol parameters with the spectral values of $`\tau _a`$ provided by CAMS. 

Regarding the sunglint correction, the main
principle is to estimate the bidirectional reflectance distribution function (BRDF) of the rough air-water interface
from the SWIR bands (i.e., ~1610 and ~2200 nm). The sunglint signal obtained in the SWIR is then extrapolated toward the
NIR and visible bands. Estimation of the sunglint radiance is based on the fact that water body is virtually totally
absorbing; water absorption coefficient in the SWIR is several orders of magnitude greater than that in the NIR. Once
corrected for atmosphere diffuse radiance, the remaining radiance in the SWIR is interpreted as the pure surface
component of the signal and then translated into BRDF. This BRDF in the SWIR is extrapolated to the other bands
considering the spectral variation of the refractive index of water and its important consequences onto the spectral
sunglint signal (see [Harmel et al., 2018](https://www.sciencedirect.com/science/article/pii/S0034425717304856) for details). The sunglint radiation is calculated for each pixel, for each
band, considering the estimated BRDF, atmosphere direct transmittance and the extraterrestrial sun radiance reaching the
atmosphere, and the water-leaving radiance is then corrected by removing this value. 

The water-leaving component at the
water surface level is eventually obtained after division by the total transmittance (i.e., diffuse + total
transmittances) calculated for the bimodal aerosol model from the LUT. The version used here accounts for the spectral
response of each band of Sentinel-2 A and B as well as Landsat-8 and it is based on the CAMS aerosol data for the
spectral value of $`\tau _a`$.


## Getting Started

> **Installing on CNES machines (TREX, HAL, PBS/SLURM clusters)?** See
> [README_CNES.md](README_CNES.md) for CNES-specific installation and run instructions.

## Installation
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Download the LUT files:
click  on [grsdata](https://drive.google.com/drive/folders/1N0-FtW-PTPblR4z-82fFrUTekMd8e3Vz?usp=sharing)
 to download and save in your desired path (your_GRSDATA_PATH) 

### please use conda environment
``` 
conda activate "name of your conda env"
```

Python >= 3.9 is recommended, example:
``` 
conda create python=3.10 -n grs_v2
conda activate grs_v2
```
Then, install python dependencies:
``` 
conda install -c conda-forge eoreader cdsapi netCDF4 docopt xmltodict numba
```
Set the `config.yml` file:
```
path:
  grsdata: your_GRSDATA_PATH
``` 

Finally, install grs with:
```commandline
pip install .
```

## Testing <a name="testing"></a>

After installation, you can type:
```commandline
grs -h
```

You should see something like:
```commandline
Executable to process Sentinel-2 L1C images for aquatic environment

Usage:
  grs <input_file> [--cams_file file] [-o <ofile>] [--odir <odir>] [--resolution res] [--scale_aot factor]   [--levname <lev>] [--no_clobber] [--allpixels] [--surfwater file] [--dem_file file] [--snap_compliant]
  grs -h | --help
  grs -v | --version

Options:
  -h --help        Show this screen.
  -v --version     Show version.

  <input_file>     Input file to be processed

  --cams_file file     Absolute path of the CAMS file to be used (mandatory)

  -o ofile         Full (absolute or relative) path to output L2 image.
  --odir odir      Ouput directory [default: ./]
  --levname lev    Level naming used for output product [default: L2Agrs]
  --no_clobber     Do not process <input_file> if <output_file> already exists.
  --resolution res  spatial resolution of the scene pixels
  --allpixels      force to process all pixels whatever they are masked (cloud, vegetation...) or not
  --surfwater file  Absolute path of the surfwater geotiff file to be used
  --dem_file file  Absolute path of the DEM geotiff file (already subset for the S2 tile)
  --scale_aot factor scaling factor applied to CAMS aod550 raster
                    [default: 1]
  --opac_model name  Force the aerosol model (OPAC) to be 'name'
                    (choice: ['ARCT_rh70', 'COAV_rh70', 'DESE_rh70',
                                'MACL_rh70', 'URBA_rh70'])
  --snap_compliant  Export output to netcdf aligned with "beam" for ESA SNAP software

  Example:
      grs /data/satellite/S2/L1C/S2B_MSIL1C_20220731T103629_N0400_R008_T31TFJ_20220731T124834.SAFE --cams_file /data/satellite/S2/cnes/CAMS/2022-07-31-cams-global-atmospheric-composition-forecasts.nc --resolution 60
  For CNES datalake:
      grs /work/datalake/S2-L1C/31TFJ/2023/06/16/S2B_MSIL1C_20230616T103629_N0509_R008_T31TFJ_20230616T111826.SAFE --cams_file /work/datalake/watcal/ECMWF/CAMS/2023/06/16/2023-06-16-cams-global-atmospheric-composition-forecasts.nc --odir /work/datalake/watcal/test --resolution 20 --dem_file /work/datalake/static_aux/MNT/COP-DEM_GLO-30-DGED_S2_tiles/COP-DEM_GLO-30-DGED_31TFJ.tif
```

> If you are on TREX CNES, see [README_CNES.md](README_CNES.md#running-grs-on-trex-with-a-slurm-interactive-job)
> for how to run the example above as a SLURM interactive job.

### To download CAMS data
[Register](https://apps.ecmwf.int/registration/) and [ask for a key](https://confluence.ecmwf.int/display/WEBAPI/Accessing+ECMWF+data+servers+in+batch#AccessingECMWFdataserversinbatch-key) to use ECMWF API

## Running the tests
From terminal:
```
grs test/data/S2B_MSIL1C_20180927T103019_N0206_R108_T31TGK_20180927T143835.SAFE --shape test/data/shape/SPO04.shp --odir test/results/ --aerosol cams_forecast --dem --resolution 20
```

You should get something like:

![image_output](images/example_snap_grs_image.png)

Another examples of output images before (1st column) and after  (2nd column) sunglint correction:

![image_output](images/Fig_valid_qualit_sea_scale.png)

> For how to launch GRS with Docker on CNES machines, see
> [README_CNES.md](README_CNES.md#running-grs-with-docker-on-cnes-deprecated).

## Deployment

See examples in [exe](exe).

## Compile Docker image locally
First and foremost, you must have the coresponding version of GRSdriver at the same level as grs2.
You should also make sure that the grsdata folder is full (it is a LTS).
Eventually, you must get Dockerfile out of grs2 folder to have a structure as diplayed below.

head_folder  
 ├grs2  
 ├GRSdriver  
 └Dockerfile

Note that anything in this folder tree will be added to the Docker build context, so make it light.
You might consider removing the notebooks and all useless files and directories from the grs2 & GRSdriver
folders to make the resulting image as light as possible (.git, illustration, notebook...).

Once all those requirements are met, you can compile the Docker image using the following command:
```
docker build -t grs2:<version_tag> *path_to_head_folder* -f *path_to_Dockerfile*
```

When the compilation has ended, you can access the image with the command:
```
docker images
```

To run the Docker image in a container on a S2 raster you can use the run_docker.sh script as follow:
```
./run_docker.sh <image_ID> <S2_raster_path> <CMAS_data_path> <desired_name_for_output> <desired_path_for_output> <desired_resolution> <surfwater_tif_path>
```
Example:
```
./run_docker.sh grs2:V2_CNES \
/DATA/S2_raster/S2B_MSIL1C_20220228T102849_N0400_R108_T31TFJ_20220228T123819.SAFE \
/DATA/CAMS/2022-02-28-cams-global-atmospheric-composition-forecasts.nc \
S2B_L2Agrs_20220228T102849_N0400_R108_T31TFJ_20220228T123819 \
/DATA/grs_outputs \
60 \
/DATA/Surfwater/SURFWATER_OPTICAL-SINGLE_T31TFJ_20220228T103850_20220228T103850_1-0-4_06/SURFWATER_OPTICAL-SINGLE_T31TFJ_20220228T103850_20220228T103850_1-0-4_06.tif
```

The docker containers will be called grs2, which mean that you cannot currently launch multiple ones simultaneously.
You can adapt the sh script to modify this behaviour.

## How the CI/CD Pipeline Works

GRSprocessor is developed on GitHub and mirrors part of its pipeline to an internal CNES
GitLab instance for HPC-based validation. Two independent pipelines run on every push.

### GitHub Actions (public, [`.github/workflows/`](.github/workflows))

**[`main.yml`](.github/workflows/main.yml) — Pipeline CI** (push to `feature/*`, `main`,
`develop`, tags; and pull requests to `main`/`develop`):

1. `python-tests` — installs GDAL and the Python dependencies, then runs the integration
   test suite (`grs/tests/integration_tests/`) with coverage. On a tag, it also regenerates
   and commits the coverage badge.
2. `pypi` (tags only, needs `python-tests`) — builds the sdist/wheel and publishes it to
   PyPI.
3. `podman-build` (needs `python-tests`) — builds the Docker image (`Dockerfile-github`);
   on a tag, it also pushes `guillaumeeb/grs:<tag>` and `:latest` to Docker Hub.
4. `lint` (pull requests and tags, needs `podman-build`) — runs `ruff` and `mypy` and
   uploads their reports as artifacts (non-blocking).
5. `sync` (needs `podman-build`) — strips Git-LFS pointer files from the whole history with
   `git-filter-repo`, then force-pushes the cleaned repository to the internal GitLab
   mirror (`gitlab.cnes.fr/waterquality/grs2.git`). This is what triggers the GitLab
   pipeline described below.

**[`doc.yml`](.github/workflows/doc.yml) — Pipeline documentation** (same push triggers):
builds the Sphinx documentation ([`docs/`](docs), furo theme, autoapi, mermaid diagrams);
on a tag it uploads the HTML as an artifact, and on `main` it publishes it to
[GitHub Pages](https://cnes.github.io/GRSprocessor/).

### GitLab CI (internal CNES, [`.gitlab-ci.yml`](.gitlab-ci.yml))

Runs on the mirrored repository at CNES. See
[README_CNES.md](README_CNES.md#gitlab-ci-internal-cnes) for the detailed stage-by-stage
breakdown.

In short: **GitHub Actions runs the public-facing pipeline** (tests, lint, PyPI, Docker
Hub, docs), while **the GitLab mirror runs CNES-internal validation** (HPC functional
tests against reference data, code quality/security gates) that cannot run outside the
CNES network.

## How to Contribute

Contributions are welcome! Quick start:

```bash
git clone https://github.com/CNES/GRSprocessor.git
cd GRSprocessor
conda create -n grs_dev python=3.11
conda activate grs_dev
pip install -r requirements.txt
pip install -e .[dev]
pytest grs/tests/integration_tests/
```

1. Create a branch off `develop` named `feature/<short-description>`.
2. Make your changes, with tests where relevant, and make sure they pass locally.
3. Open a pull request targeting `develop` — the `python-tests` job runs automatically on
   the PR (see [above](#how-the-cicd-pipeline-works)).
4. A maintainer will review your PR.

See [CONTRIBUTING.md](https://github.com/CNES/GRSprocessor/blob/main/CONTRIBUTING.md) for
the full guide: reporting issues, code style (`black`/`isort`, `ruff`/`mypy` in CI),
building the documentation, and the code of conduct.

## Authors

* **Tristan Harmel** - *Initial work* - [contact](mailto:tristan.harmel@ntymail.com)

See also the list of [contributors](https://github.com/CNES/GRSprocessor/graphs/contributors) who participated in this project.

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](https://github.com/CNES/GRSprocessor/blob/main/LICENSE) file for details

## Acknowledgments

* The [Step forum](http://forum.step.esa.int) and Marco Peters are acknowledged for their useful help to process Sentinel-2 data
with the snappy API.
* The authors are very grateful to Olivier Hagolle
for providing open source codes to perform gaseous absorption correction and massive Sentinel-2 data download.

