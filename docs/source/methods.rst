Methods
=======

This page describes what :py:meth:`Process.execute <grs.grs_process.Process.execute>` computes, step
by step. The algorithm is described in `Harmel et al. (2018)
<https://www.sciencedirect.com/science/article/pii/S0034425717304856>`__.

Notation
--------

.. list-table::
   :widths: 25 75

   * - :math:`\lambda`
     - central wavelength of the band
   * - :math:`\theta_s,\ \theta_v,\ \phi`
     - solar zenith, viewing zenith and relative azimuth angles
   * - :math:`R_{TOA}(\lambda)`
     - top-of-atmosphere reflectance of the Level-1 image
   * - :math:`\tau_{ref}`
     - aerosol optical thickness (AOT) at 550 nm
   * - :math:`\tau_a(\lambda),\ \tau_R(\lambda)`
     - aerosol and Rayleigh optical thicknesses
   * - :math:`P,\ P_0`
     - surface pressure and reference pressure (:math:`P_0` = 101500 Pa)
   * - :math:`T_g,\ T_g^{diff}`
     - gaseous transmittance of the direct and of the diffuse light
   * - :math:`U`
     - wind speed (m s\ :sup:`-1`)
   * - :math:`\varepsilon(\lambda)`
     - mean spectral shape of the sunglint reflectance

Processing steps
----------------

1. **Level-1 image.** Sentinel-2 ``.SAFE`` and Landsat 8/9 products are read by
   `GRSdriver <https://github.com/Tristanovsk/GRSdriver>`__ at the requested ``resolution``:
   :math:`R_{TOA}(\lambda)` and the angles :math:`\theta_s, \theta_v, \phi` (the viewing angles of
   Sentinel-2 depend on the band). A Level-1 product already converted to NetCDF can be given instead.
2. **Ancillary data.** The CAMS file gives the surface pressure, the gas columns, the spectral AOD and
   the 10 m wind (:py:class:`~grs.cams.CamsProduct`).
3. **Pixel classification.** :py:class:`~grs.mask.Masking` builds the ``flags`` bit field (no data,
   clouds with s2cloudless or the Landsat QA band, water indices, cirrus, high SWIR, optional
   SURFWATER mask; see the :doc:`bit table <api/mask>`).
4. **Gaseous absorption (O**\ :sub:`3`, **NO**\ :sub:`2`, **background gases).** The TOA reflectance is
   divided by the transmittance computed from the CAMS columns (`Gaseous transmittance`_).
5. **Water pixels.** With the band index

   .. math::

      NDWI_{SWIR} = \frac{R(490) - R(1610)}{R(490) + R(1610)}

   the processed pixels are those with :math:`NDWI_{SWIR} > 0` and :math:`R(2190) < 0.2` (unless
   ``allpixels=True``).
6. **Aerosol model.** Unless ``opac_model`` is given, the OPAC model :math:`m` whose spectral AOD is
   closest to the CAMS one is selected (the urban model is excluded):

   .. math::

      \hat m = \underset{m}{\operatorname{argmin}} \sum_\lambda
      \left| \frac{\overline{\tau_{CAMS}(\lambda)}}{\overline{\tau_{CAMS}(550)}} - \tau_m(\lambda; \tau_{ref} = 1) \right|

   where the bars are averages over the image.
7. **Aerosol optical thickness.** :math:`\tau_{ref}` is the AOD at 550 nm of CAMS, multiplied by
   ``scale_aot`` and interpolated to the pixels. It is not retrieved from the image in this version.
8. **Look-up tables.** The radiative transfer LUTs are reduced to the selected model, to the mean wind
   speed of the image, to the geometries of the image and to the range of :math:`\tau_{ref}`
   (`Atmospheric diffuse light`_).
9. **Surface pressure.** With a DEM of elevation :math:`z` (m), the pressure is computed from the CAMS
   mean sea level pressure :math:`P_{msl}`:

   .. math::

      P = P_{msl} \left(1 - \frac{0.0065\, z}{288.15}\right)^{5.255}

   otherwise the CAMS surface pressure is used.
10. **Atmospheric and sunglint correction** of every pixel, by blocks of 512 x 512 pixels processed in
    parallel by ``ncpu`` processes (`Remote-sensing reflectance`_).
11. **Flags and mask.** Bit 18 (``neg_rrs``) is set where :math:`R_{rs}(490)` or :math:`R_{rs}(565)` is
    below :math:`-5 \times 10^{-4}`. The binary ``mask`` is 0 (good) for the pixels with bit 3 (water
    index) set and none of the bits 0 (no data), 1 (s2cloudless cloud), 10 (SURFWATER cloud and shadow),
    13 (Landsat cloud), 14 (Landsat cloud shadow) and 18 (negative Rrs); it is 1 elsewhere
    (:py:meth:`Masking.create_mask <grs.mask.Masking.create_mask>`).
12. **Output.** :py:class:`~grs.output.L2aProduct` writes ``Rrs``, ``BRDFg``, ``aot550``, ``flags``,
    ``mask``, the angles, the DEM and the ancillary data as NetCDF (see :doc:`api/output`).

Gaseous transmittance
---------------------

For each gas, the band transmittance :math:`T_{gas}(\lambda; c)` is tabulated against the total
column :math:`c` and convolved with the spectral response of the bands
(:py:class:`~grs.acutils.GaseousTransmittance`). Step 4 corrects for the gases absorbing high in the
atmosphere, above the scattering layer:

.. math::

   R_{TOA}(\lambda) \leftarrow \frac{R_{TOA}(\lambda)}{T_{O_3}(\lambda)\, T_{NO_2}(\lambda)\, T_{bg}(\lambda)}

(:math:`T_{bg}`: CO, CO\ :sub:`2`, O\ :sub:`2` and O\ :sub:`4`, for the CAMS surface pressure). Water
vapour and methane, mostly in the lower layers, are taken into account in the radiative transfer
equations below through

.. math::

   T_g(\lambda) = T_{H_2O}(\lambda; C_{H_2O})\, T_{CH_4}(\lambda; C_{CH_4}), \qquad
   T_g^{diff}(\lambda) = T_{H_2O}(\lambda; 0.5\, C_{H_2O})\, T_{CH_4}(\lambda; 0.5\, C_{CH_4})

where the factor 0.5 accounts for the part of the diffuse light path that is above the gases.

Atmospheric diffuse light
-------------------------

The LUT of the selected aerosol model gives the reflectance :math:`R_{diff}` of the atmosphere and of
the sky light reflected by the air-water interface, as a function of the wavelength, of
:math:`\tau_{ref}` and of the geometry. It is removed from the TOA signal after scaling for the
pressure and the gaseous absorption:

.. math::

   R_{corr}(\lambda) = R_{TOA}(\lambda) - T_g^{diff}(\lambda)\, \frac{P}{P_0}\,
   R_{diff}(\lambda; \tau_{ref}, \theta_s, \theta_v, \phi)

Transmittances
--------------

With the two-way air mass :math:`M = 1/\cos\theta_s + 1/\cos\theta_v`:

.. math::

   T_{dir}(\lambda) = \exp\left[-\left(\tau_a(\lambda) + \tau_R(\lambda)\,\frac{P}{P_0}\right) M\right]

.. math::

   T_{tot}(\lambda) = T_{\downarrow}(\lambda; \theta_s)\; T_{\downarrow}(\lambda; \bar\theta_v)^{1.05}\; T_g(\lambda)

where :math:`T_{\downarrow}` is the total (direct + diffuse) transmittance of the downward irradiance
from the LUT and :math:`\bar\theta_v` the viewing angle averaged over the bands.

Sunglint
--------

The sea surface is described by the isotropic slope distribution of Cox and Munk (1954), whose mean
square slope depends on the wind speed:

.. math::

   \sigma^2 = \frac{U + 0.586}{195.3}

The water-leaving signal is negligible in the SWIR: the sunglint is measured in the SWIR bands
:math:`\lambda_s` (1610 and 2190 nm for Sentinel-2) and extrapolated to all the bands with its spectral
shape :math:`\varepsilon(\lambda)`, relative to the reference band :math:`\lambda_r` (the last SWIR band):

.. math::

   R_f = \min_{\lambda_s} \max\left(0,\ \frac{\varepsilon(\lambda_r)\, R_{corr}(\lambda_s)}
              {T_g(\lambda_s)\, T_{dir}(\lambda_s)\, \varepsilon(\lambda_s)}\right)

.. math::

   R_g(\lambda) = T_g(\lambda)\, T_{dir}(\lambda)\, \varepsilon(\lambda)\, \frac{R_f}{\varepsilon(\lambda_r)}

:math:`R_f` is exported as ``BRDFg``.

.. note::

   The probability of the wave slopes reflecting the sun towards the sensor,
   :py:meth:`Product.p_slope <grs.product.Product.p_slope>`, is currently set to 1 (fix of the azimuth
   handling of version 2.1.9): the spectral extrapolation then only depends on
   :math:`\varepsilon(\lambda)` and on the transmittances, and the wind speed is only used to select
   the look-up tables.

Remote-sensing reflectance
--------------------------

.. math::

   R_{rs}(\lambda) = \frac{R_{corr}(\lambda) - R_g(\lambda)}{\pi\, T_{tot}(\lambda)}

(the division by :math:`T_{tot}` is only applied where :math:`R_{corr} \ge R_g`). Every pixel uses its
own :math:`\tau_{ref}`, pressure and gaseous transmittances.

References
----------

* Harmel, T., Chami, M., Tormos, T., Reynaud, N., Danis, P.-A. (2018). Sunglint correction of the
  Multi-Spectral Instrument (MSI)-SENTINEL-2 imagery over inland and sea waters from SWIR bands.
  *Remote Sensing of Environment*, 204, 308-321.
  `doi:10.1016/j.rse.2017.10.022 <https://doi.org/10.1016/j.rse.2017.10.022>`__
* Cox, C., Munk, W. (1954). Measurement of the roughness of the sea surface from photographs of the
  sun's glitter. *Journal of the Optical Society of America*, 44(11), 838-850.
