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

### Changed
- Coarse raster now defined by spatial resolution instead of pixel factors.
- Optimised LUT interpolation; `acutils` module renamed to `utils`.
- Lazy package loading for faster imports; Sentinel-2 viewing geometry loaded lazily.
- AOT smoothing windows scale with the retrieval blocks (new `aot_smoothing_blocks` parameter of
  `Process.execute`, default 3).

### Fixed
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
