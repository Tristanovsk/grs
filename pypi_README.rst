Usage
=====

WARNING: installation from `GitHub repository <https://github.com/Tristanovsk/grs>`__ is recommended
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Installation
~~~~~~~~~~~~~~~~~~~~~~~

These instructions will get you a copy of the project up and running on
your local machine for development and testing purposes. See deployment
for notes on how to deploy the project on a live system.

Download the LUT files
-----------------------

click on
`grsdata <https://drive.google.com/drive/folders/1N0-FtW-PTPblR4z-82fFrUTekMd8e3Vz?usp=sharing>`__
to download and save in your desired path (your_GRSDATA_PATH)

Installation with conda environment
------------------------------------

::

   conda activate "name of your conda env"

Python >= 3.9 is recommended, example:

::

   conda create python=3.10 -n grs_v2
   conda activate grs_v2

Then, install python dependencies:

::

   conda install -c conda-forge eoreader cdsapi netCDF4 docopt xmltodict numba

Set the ``config.yml`` file:

::

   path:
     grsdata: your_GRSDATA_PATH

Finally, install grs with:

.. code::

   pip install .

Testing 
~~~~~~~~

After installation, you can type:

.. code::

   grs -h

You should see something like:

.. code::

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
     --resolution=res  spatial resolution of the scene pixels
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
