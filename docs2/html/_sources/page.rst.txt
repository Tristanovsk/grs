Documentation for the GRS algorithm.
------------------------------------

GRS (Glint Removal for Sentinel-2-like sensors)
_______________________________________________

The GRS processing was designed to correct Sentinel-2-like satellite images for atmospheric effects and Sun reflection (sunglint) above water surfaces (e.g., ocean, coastal, estuary, river and lake waters).
Theoretical background can be found in `[Harmel et al., 2018] <https://www.sciencedirect.com/science/article/pii/S0034425717304856/>`_.

The overall structure of the algorithm can be summarized by the following flowchart:

.. figure:: ../../images/flowchart_sunglint_S2.png
  :scale: 60%
  :alt: check `README <../../README.md>`_

  Flowchart of the proposed GRS algorithm to remove the sunglint contribution from the MSI/Sentinel-2 images.


Getting Started
_______________

These instructions will get you a copy of the project up and running on
your local machine for development and testing purposes. See deployment
for notes on how to deploy the project on a live system.

Prerequisites
~~~~~~~~~~~~~

`Register <https://apps.ecmwf.int/registration/>`__ and `ask for a
key <https://confluence.ecmwf.int/display/WEBAPI/Accessing+ECMWF+data+servers+in+batch#AccessingECMWFdataserversinbatch-key>`__
to use ECMWF API

Download and install the `SNAP
software <http://step.esa.int/main/download/>`__. Configure the
`SNAP-python
interface <https://senbox.atlassian.net/wiki/spaces/SNAP/pages/50855941/Configure+Python+to+use+the+SNAP-Python+snappy+interface>`__
and link the obtained ``snappy`` folder to your python site-packages as
``esasnappy``. For example:

.. code:: bash

    ln -s /FULL_PATH/.snap/snap-python/snappy /PATH_TO_LIB_PYTHON/lib/python3.6/site-packages/esasnappy

Compilers such gcc and gfortran are needed to install the package.

Compile all C and fortran files into shared libraries:

::

    make


Generate the `config.py` file:

 * In the ./grs/grs folder, copy `config_local.py` to `config.py`.

 * Then, edit `config.py` according to your folders tree and path to your grs installation folder.

Installing
~~~~~~~~~~

To install the package:

::

    python setup.py install

or

::

    python setup.py install --user

If the installation is successful, you should have:

::

    $ grs
    Usage:
      grs <input_file> [--sensor <sensor>] [-o <ofile>] [--odir <odir>] [--shape <shp>] [--wkt <wktfile>]   [--longlat <longmax,longmin,latmax,latmin> ]    [--altitude=alt] [--dem] [--aerosol=DB] [--aeronet=<afile>]    [--aot550=aot] [--angstrom=ang] [--output param]   [--resolution=res] [--levname <lev>] [--no_clobber] [--memory_safe] [--unzip]
      grs -h | --help
      grs -v | --version


Running the tests
~~~~~~~~~~~~~~~~~

From terminal:

::

    grs test/data/S2B_MSIL1C_20180927T103019_N0206_R108_T31TGK_20180927T143835.SAFE --shape test/data/shape/SPO04.shp --odir test/results/ --aerosol cams_forecast --dem --resolution 20

You should get something like:

.. figure:: ../../images/example_snap_grs_image.png
   :alt: image\_output

   Example of an output image visualized in the `SNAP software <http://step.esa.int/main/download/>`__.


Another examples of output images before (1st column) and after (2nd
column) sunglint correction:

.. _fig3:
.. figure:: ../../images/Fig_valid_qualit_sea_scale.png
   :alt: image\_output

   RGB images obtained after subtraction of the atmospheric radiance from TOA signal but before (left column) and after (right column) removing the sunglint radiance. These images correspond to the areas surrounding the AERONET-OC sites of (a, b) Venice (July 18, 2016) and (c, d) WaveCIS (April 23, 2016). Note that the same color scale was used to generate the RGB images before and after removing the sunglint radiance.


Deployment
~~~~~~~~~~

See examples in :ref:`exe <deployment_modules.rst>`.


Package development
____________________

Contributing
~~~~~~~~~~~~

Please contact `authors <tristan.harmel@gmail.com>`_ for details on our code of conduct, and the process for submitting pull requests to us.

Authors
~~~~~~~~~~

* **Tristan Harmel** - *Initial work* - `contact <tristan.harmel@gmail.com>`_

See also the list of [contributors](...) who participated in this project.

License
~~~~~~~

This project is licensed under the MIT License - see the :download:`LICENSE.md <../../LICENSE.md>` file for details.

Acknowledgments
~~~~~~~~~~~~~~~

* The `Step forum <http://forum.step.esa.int>`_ and Marco Peters are acknowledged for their useful help to process Sentinel-2 data with the snappy API.

* The authors are very grateful to Olivier Hagolle for providing open source codes to perform gaseous absorption correction and massive Sentinel-2 data download.
