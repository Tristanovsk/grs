# CAMS auxiliary data

GRS takes the state of the atmosphere at the time of the image from the
[CAMS global atmospheric composition forecasts](https://ads.atmosphere.copernicus.eu/datasets/cams-global-atmospheric-composition-forecasts)
(Copernicus Atmosphere Monitoring Service), distributed by the Atmosphere Data Store (ADS): surface
pressure, total columns of the absorbing gases, spectral aerosol optical depth and 10 m wind. One file
covering the date and the footprint of the image is given to the processing (`cams_file`,
`--cams_file`).

The notebook shows:

- the access to the ADS (account, licence of the dataset, API key in `~/.cdsapirc`),
- the download of the data of one image with
  {meth}`CamsProduct.cams_download <grs.cams.CamsProduct.cams_download>`, and of global daily files,
  stored in a `YYYY/MM/DD` tree, for the processing of many images
  (`<cams_dir>/YYYY/MM/DD/YYYY-MM-DD-cams-global-atmospheric-composition-forecasts.nc`; older archive:
  one file per month),
- how to find the CAMS file of an image in this archive,
- what GRS reads from a CAMS file ({meth}`CamsProduct.load <grs.cams.CamsProduct.load>`): maps of the
  variables over a Sentinel-2 tile, spectral aerosol optical depth and Ångström exponent (used to
  choose the aerosol model), wind over the image.

The download cells are not run in this documentation (they need an ADS account); the rest of the
notebook is rendered with the outputs of its last run on a Sentinel-2 image of 2026-03-09 over Lake
Geneva (tile 31TGM), with the CAMS daily file of its date. The
[`ecmwf/`](https://github.com/CNES/GRSprocessor/tree/main/ecmwf) folder of the repository holds the
scripts used to download long periods on a computing cluster, and the daily update of the archive
(`download_cams_cds_day.py`).

```{toctree}
:maxdepth: 1

auxdata/cams_data
```
