# Processing a Sentinel-2 image

This notebook runs GRS on one Sentinel-2 L1C image from Python, from the input files to the checks of
the L2A product:

- inputs: the Level-1 image, the CAMS file of the month and an optional DEM; the look-up tables (see
  {doc}`../configuration`),
- processing with {class}`grs.Process <grs.grs_process.Process>`:
  {meth}`~grs.grs_process.Process.execute` sets up the processing lazily (dask), and
  {meth}`~grs.grs_process.Process.write_output` computes it block by block and writes the product,
- the L2A product as a Zarr store with overviews (or NetCDF, from the extension of the output path),
- maps of the remote-sensing reflectance $R_{rs}$, of the sunglint `BRDFg` and of the aerosol optical
  thickness retrieved from the image, compared with the CAMS forecast,
- the pixel classification (`flags`) and the `mask` of the valid water pixels,
- $R_{rs}$ spectra from the turbid Loire estuary to the open sea.

```python
import grs

process = grs.Process()
process.execute(l1c_file, cams_file=cams_file, resolution=20, aot_resolution=240,
                odir='S2B_MSIL2AGRS_20200411T110619_T30TWT.zarr')
process.write_output()
```

The image is the Sentinel-2B acquisition of 2020-04-11 over the Loire estuary and the Bay of Biscay
(tile 30TWT). The notebook is rendered with the outputs of its last run; the input images are not
distributed with GRS: set the paths of the settings cell to run it on your own data. The same
processing from the command line is described in {doc}`../usage`.

```{toctree}
:maxdepth: 1

processing/grs_s2_processing
```
