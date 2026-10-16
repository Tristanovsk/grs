# Changelog

All notable changes to GRS are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Dates are given for tagged releases only.

## [Unreleased] - 3.1.0

### Added
- New masking procedure including OmniCloudMask (building on the former s2cloudless development).
- Export to Zarr (Zarr v3) with multiscale pyramids and geo-information.
- AOT coarse resolution stored in output attributes (default 480 m, a multiple of the 10, 20, 30 and
  60 m pixel sizes).
- `--aot_resolution` option of the CLI (pixel size of the coarse grid of the aerosol retrieval, default
  480 m).
- Documentation of the kernel physics with equations; reworked Sphinx documentation for Read the Docs
  (light/dark logos, partner logos, reorganised example notebooks: processing, auxiliary data,
  visualization, Landsat).
- `CONTRIBUTING.md` and a separate `README_CNES.md`.

- `grs_config` command: prints the configuration in use and checks the LUTs (exit code 1 if one is
  missing); `--init [DATA_ROOT]` writes the user configuration file; `--download` downloads the
  LUTs from Zenodo (MD5 checked) into the LUT folder.
- Configuration page of the documentation.
- Automatic download of the missing LUTs from Zenodo before the first processing
  (`grs.config.ensure_luts`); switched off with `zenodo: auto_download: false` or `GRS_NO_DOWNLOAD=1`.
- PyPI packaging as `GRSprocessor` (package and commands still `grs`): version read from
  `grs.__version__`, dependencies matching the imports, wheel without tests nor LUTs; release recipe in
  `CONTRIBUTING.md` and GitHub workflow `pypi.yml` (trusted publishing on `v*` tags).

### Changed
- Configuration read at run time from the packaged defaults (`grs/config.yml`, now documented), the
  user file `~/.config/grs/config.yml` (or `GRS_CONFIG`) and `GRS_DATA_ROOT`, instead of editing
  `grs/config.yml` before installing (`grs.config` module). When `grsdata` is not set, the LUTs are
  read from the user data folder (`~/.local/share/grs/lut` on Linux). New dependency: `platformdirs`.
- `Process.execute` checks the LUT files before reading the image.
- Coarse raster now defined by spatial resolution instead of pixel factors.
- Optimised LUT interpolation; `acutils` module renamed to `utils`.
- Lazy package loading for faster imports; Sentinel-2 viewing geometry loaded lazily.
- AOT smoothing windows scale with the retrieval blocks (new `aot_smoothing_blocks` parameter of
  `Process.execute`, default 3).
- Python >= 3.11. The CLI reads the cloud cover of Sentinel-2 images from `MTD_MSIL1C.xml` directly:
  the GDAL Python bindings are no longer needed.
- `execute` and `write_output` use `processor: ncpu` dask threads (one per core before), the L1C bands
  are no longer held by the CAMS product, and the water masking is done in place: about 35 % less
  memory and 20 % less time for a Sentinel-2 tile at 20 m.
- Two glibc malloc arenas (`MALLOC_ARENA_MAX=2`) in the Docker images, and set at start-up by the `grs`
  command and `exe/launcher.py` (unless `MALLOC_ARENA_MAX` is set): the memory freed by the reading and
  processing threads is reused, about 2-3 GB less peak memory with the parallel band reading of GRSdriver.
- Blocks without water pixels are skipped in the final processing; the disabled sunglint slope
  probability is no longer computed.

### Fixed
- NetCDF export failed (`TypeError: illegal data type for attribute b'gas_absorption_correction'`):
  the boolean attributes are written as integers (0/1) in the NetCDF files; Zarr keeps the booleans.
- `Process.cams_dir` was set from the `trans_lut` setting instead of `path: cams`.
- NaN `Rrs` on the edges of the image (about half a pixel of the AOT grid), where the AOT map was not
  extrapolated.
- Vertical bands of NaN in Landsat products, and biased sunglint normalisation in the Landsat AOT
  retrieval (slope probability of the reference band taken from a row of the processing block).
- Weights of the AOT smoothing: the number of water pixels of each retrieval block is used as weight
  (it was overwritten, so that only the empty blocks were weighted); outliers and empty blocks get the
  median AOT with a null weight.
- `Kernel.filter2d`: window wrapping around the image edges and in-place computation from already
  smoothed values.
- `--no_clobber` behaviour for Zarr output and `--snap_compliant` help text in the CLI.
- Unresolved imports with lazy package loading.
- Autosummary class template path in the API reference.
- Compatibility with recent xarray versions (air-mass broadcasting in the kernel, `coarsen().mean()`
  in masking).
- Gas absorption correction of dask-backed bands (e.g. GRSdriver with `EOREADER_USE_DASK=1`), which
  failed on the `spatial_ref` coordinate of the CAMS grid.

## [3.0.0]

### Changed
- Merge of the 2.2.1 scientific developments with the CI/CD pipeline from 2.1.9.

### Fixed
- Bugs related to true (band-averaged) wavelengths.
- Support for the new ADS CAMS format.
- CLI: accept `--odir` as long alias for `-o`; pass a string path to `gdal.Open`.
- Docker: upgrade datashader to match the dask API.

## [2.2.1]

### Changed
- Refinement of the aerosol retrieval algorithm.

## [2.2.0]

### Changed
- Major update of rugosity (surface roughness) and aerosol retrievals.

## [2.1.9] - 2026-05-27

### Added
- L1C NetCDF accepted as input in the CLI.
- CI/CD: GitLab synchronisation, PyPI publication, Sphinx documentation build, coverage badge,
  CWL descriptor.

### Fixed
- Azimuth handling in the sunglint correction.

## [2.1.8] - 2025-03-19

### Changed
- Improved logging and CI.

## [2.1.7] - 2024-11-07

### Changed
- Account for absorption/scattering interaction in the transmittance of absorbing gases.

## [2.1.6] - 2024-09-24

### Changed
- Use `GRSdriver` instead of `s2driver`.

## [2.1.5]

### Added
- Landsat 8 and 9 imagery support.

### Changed
- Repackaging of some classes.

## [2.1.4]

### Changed
- `grs` process can be called as a simple function within a script (new input/output handling).
- Output variables now include bitmask flags.

## [2.1.3]

### Changed
- New automated CAMS loading and revised parameters.

## [2.1.2]

### Fixed
- Tiles straddling the Greenwich meridian.

## [2.1.1]

### Changed
- Package data.

## [2.1.0] - 2024-02-22

### Changed
- Major update of the radiative transfer look-up tables and aerosol models.

## [2.0.5]

### Changed
- New output format compliant with the SNAP "BEAM" format.

## [2.0.4]

### Changed
- Output format adapted for QGIS.

## [2.0.3]

### Added
- SurfWater integration.

## [2.0.2]

### Changed
- Memory optimisation.

## [2.0.1]

### Fixed
- I/O issues.

## [2.0.0]

### Removed
- snappy skeleton (SNAP dependency); simplification of previous options; Landsat handling dropped.

## [1.5.0]

### Changed
- Output parameters extended with `ndwi_nir` and `ndwi_swir`.

## [1.4.0]

### Changed
- Image processed by rectangular chunks.

## [1.3.4]

### Added
- Slope and shade from DEM.

## [1.3.3]

### Fixed
- Sunglint BRDF output.

## [1.3.2]

### Changed
- Major improvement of LUT interpolation and access.

## [1.3.1]

### Fixed
- CAMS data extraction with xarray and longitude conventions.

## [1.3.0]

### Changed
- CAMS data update (CDS version): new AOD wavelengths and spectral SSA from 2018 onwards.
- Adjustment for absorbing aerosols through SSA.

## [1.2.4]

### Added
- Option to load MAJA and WaterDetect masks and export them in the output file.
- Option to process WaterDetect water pixels only.

## [1.2.3]

### Changed
- Compliance with SNAP 8.
- Output written directly in compressed NetCDF4 at the end of the process.

### Added
- New utility `get_subset`.

## [1.2.2]

### Changed
- Nearest-neighbour interpolation to keep tile-edge pixels.
- `product.dispose` to minimise memory usage in the SNAP JVM.

## [1.2.1]

### Added
- High-latitude processing (DEM change).
- Option to process all pixels before masking non-water pixels.

## [1.2.0]

### Changed
- Full data matrix loaded from the image instead of line-by-line pixel extraction
  (preparation for multi-pixel retrieval).

## [1.1.4]

### Added
- Cloud mask; compressed NetCDF4 output format.

## [1.1.3]

### Changed
- Minor changes.

## [1.1.2]

### Changed
- Output in Rrs units.

## [1.1.1]

### Changed
- Minor changes.

## [1.1.0]

### Added
- Landsat 4, 5, 7 and 8 support.

## [1.0.0]

### Added
- GRS for Sentinel-2.
