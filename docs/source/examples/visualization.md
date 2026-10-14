# Visualization

The L2A products of GRS are standard Zarr stores (or NetCDF files): they can be opened, plotted and
subset with xarray, rioxarray and matplotlib only. For the interactive exploration of the images, use
the GRS toolbox [grstbx](https://github.com/Tristanovsk/grstbx) (version 2.3.1 or later). It is not on
PyPI; install it from GitHub, with the Jupyter extras of the interactive viewers:

```bash
pip install "grstbx[notebook] @ git+https://github.com/Tristanovsk/grstbx.git"
```

Both notebooks are rendered with the outputs of their last run: the interactive maps are not
functional in this documentation (no Python kernel behind them), and the input images are not
distributed. Set the input paths to run them on your own products.

## Reading an L2A product with xarray

The first notebook opens the product of the {doc}`Sentinel-2 tutorial <processing>` with xarray:

- layout of a Zarr product: full resolution (group `0`), overviews (`1`, `2`...) and `ancillary`
  group, described by the root attribute `multiscales`,
- true-colour quick-looks from an overview level,
- decoding of the `flags` bit field from its attributes and `mask` of the valid water pixels,
- extraction of an area from longitude / latitude bounds with `rio.clip_box` (only the chunks of the
  area are read), spectra of a few pixels and export of $R_{rs}$ to GeoTIFF.

## Exploring an L2A image with grstbx

The second notebook (from the
[grstbx tutorials](https://github.com/Tristanovsk/grstbx/tree/main/docs/source/tutorials)) uses the
tools of grstbx:

`grstbx.L2grs`
: loads one or several L2A products (Zarr or NetCDF) into a datacube along `time`; for the Zarr
  products, `level` selects the pyramid level (0 = full resolution).

`grstbx.Masking`
: decodes the `flags` bit field (names, descriptions, proportion of raised pixels) and builds masks
  from the bits.

`grstbx.visual.ViewSpectral`
: interactive map of $R_{rs}$ over a free basemap (CARTO, OpenStreetMap, OpenTopoMap, Esri...).
  Only the visible pixels are read, at screen resolution, from the Zarr pyramid. It shows true-colour
  or single-band images with colormaps adapted to water reflectance, the full-resolution spectra of
  the drawn points and, since grstbx 2.3.1, the spectrum of the pixel under the mouse (dashed line,
  *Live spectrum* checkbox). The drawn points and polygons are retrieved with `get_points` and
  `get_geom`.

`grstbx.visual.ViewParam`
: same viewer for any 2-D parameter (e.g. water-quality products), with a colormap picker (`cmap=`
  for the initial colormap, *Reverse* toggle) and a log scale.

It also computes simple water-quality parameters (chlorophyll-a, CDOM absorption, suspended
particulate matter) from $R_{rs}$.

```{figure} ../_static/grstbx_visual_tool.gif
:alt: animated view of the interactive viewer of grstbx

Interactive exploration of an L2A image with grstbx.
```

```{toctree}
:maxdepth: 1

visualization/grs_l2a_xarray
visualization/grstbx_l2a_visu
```
