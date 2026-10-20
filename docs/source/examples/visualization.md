# Visualization

The L2A products of GRS are standard Zarr stores (or NetCDF files): they can be opened, plotted and
subset with xarray, rioxarray and matplotlib only. For the interactive exploration of the images, use
the GRS toolbox [grstbx](https://grstbx.readthedocs.io) (version 2.4.0 or later), installed from
[PyPI](https://pypi.org/project/grstbx/) with the Jupyter extras of the interactive viewers:

```bash
pip install "grstbx[notebook]"
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

The second notebook (adapted from the
[grstbx tutorials](https://grstbx.readthedocs.io/en/latest/tutorials/l2a_image.html)) uses the tools of
grstbx:

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
  the drawn points and the spectrum of the pixel under the mouse (dashed line, *Live spectrum*
  checkbox). The bands of the RGB composite (true colour, false colour, red edge, SWIR or any band),
  its gamma, the contrast and the brightness can be adjusted. The drawn points and polygons are
  retrieved with `get_points` and `get_geom`.

`grstbx.visual.ViewParam`
: same viewer for any 2-D parameter (e.g. water-quality products), with a colormap picker (`cmap=`
  for the initial colormap, *Reverse* toggle) and a log scale.

It then derives the water-quality parameters of the area with
[GRSl2bgen](https://grsl2bgen.readthedocs.io) (optical water types, chlorophyll-a, suspended
particulate matter and turbidity, CDOM, transparency, quality flags of the L2A reflectance) and browses
them with `ViewParam`.

The same viewers are available without notebook in the web dashboard of grstbx, which opens GRS L2A
and L2B products, time series of products or datacubes in the browser
(see the [dashboard tutorial](https://grstbx.readthedocs.io/en/latest/tutorials/dashboard.html)):

```bash
grstbx-viewer S2B_MSIL2AGRS_20200411T110619_T30TWT.zarr
```

```{figure} ../_static/grstbx_visual_tool.gif
:alt: animated view of the interactive viewer of grstbx

Interactive exploration of an L2A image with grstbx.
```

```{toctree}
:maxdepth: 1

visualization/grs_l2a_xarray
visualization/grstbx_l2a_visu
```
