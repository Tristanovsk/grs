# Processing a Landsat 8/9 image

GRS processes the Level-1 products (collection 2, `L1TP` / `L1GT`) of the OLI sensor of Landsat 8 and
9 like the Sentinel-2 images, with a few differences:

- the bands of OLI from 443 to 2201 nm (865 nm for the NIR, 1609 and 2201 nm for the sunglint), at
  30 m;
- a single viewing geometry for all the bands;
- the clouds taken from the level-1 quality band `QA_PIXEL`, copied into the `flags` of the L2A
  product (bits `l1_*`, see {doc}`../api/mask`).

The notebook works on a subset of a Landsat 9 scene of 2022-11-20 over the Rhône delta and the Berre
lagoon (France): preview of the scene with [GRSdriver](https://github.com/Tristanovsk/GRSdriver),
reading of an area of interest only, processing of this subset (passed to
{meth}`Process.execute <grs.grs_process.Process.execute>` as an `xarray.Dataset`), maps of $R_{rs}$
and of the aerosol optical thickness, level-1 flags and spectra of the lagoons and of the Rhône
plume.

```python
import geopandas as gpd
from shapely.geometry import box
import GRSdriver
import grs

aoi = gpd.GeoSeries([box(4.35, 43.25, 5.3, 43.6)], crs=4326)
l1 = GRSdriver.LandsatDriver(l1_file, resolution=30, subset=aoi)
l1.load_mask()
l1.load_product()

process = grs.Process()
process.execute(l1.prod, cams_file=cams_file, resolution=30, odir='rhone_delta.zarr')
process.write_output()
```

The notebook is rendered with the outputs of its last run; the input images are not distributed with
GRS: set the paths of the settings cell to run it on your own data.

```{toctree}
:maxdepth: 1

landsat/grs_landsat_processing
```
