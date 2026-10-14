Methods
=======

This page describes what :py:meth:`Process.execute <grs.grs_process.Process.execute>` computes, step
by step. The algorithm is described in `Harmel et al. (2018)
<https://www.sciencedirect.com/science/article/pii/S0034425717304856>`__; the implementation of the
equations is in :py:class:`~grs.grs_kernel.Kernel` (see also the :doc:`api/grs_kernel` overview).

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
     - CAMS surface pressure and reference pressure (:math:`P_0` = 101320 Pa)
   * - :math:`T_g,\ T_g^{diff}`
     - gaseous transmittance of the direct and of the diffuse light
   * - :math:`U`
     - wind speed (m s\ :sup:`-1`)
   * - :math:`\varepsilon(\lambda)`
     - mean spectral shape of the sunglint reflectance
   * - :math:`P_{slope}(\lambda)`
     - probability of the wave slopes reflecting the sun towards the sensor

Processing steps
----------------

1. **Level-1 image.** Sentinel-2 ``.SAFE`` and Landsat 8/9 products are read by
   `GRSdriver <https://github.com/Tristanovsk/GRSdriver>`__ at the requested ``resolution``:
   :math:`R_{TOA}(\lambda)` and the angles :math:`\theta_s, \theta_v, \phi` (the viewing angles of
   Sentinel-2 depend on the band). A Level-1 product already converted to NetCDF can be given instead.
2. **Ancillary data.** The CAMS file gives the surface pressure, the gas columns, the spectral AOD and
   the wind speed (:py:class:`~grs.cams.CamsProduct`).
3. **Pixel classification.** :py:class:`~grs.mask.Masking` builds the ``flags`` bit field (clouds with
   s2cloudless and OmniCloudMask, water indices, cirrus, high SWIR, optional SURFWATER mask; see the
   :doc:`bit table <api/mask>`).
4. **Gaseous absorption (O**\ :sub:`3`, **NO**\ :sub:`2`, **background gases).** The TOA reflectance is
   divided by the transmittance computed from the CAMS columns (`Gaseous transmittance`_).
5. **Water pixels.** With the band indices

   .. math::

      NDWI = \frac{R(490) - R(842)}{R(490) + R(842)}, \qquad
      NDWI_{SWIR} = \frac{R(490) - R(1610)}{R(490) + R(1610)}

   the processed pixels are those with :math:`NDWI_{SWIR} > 0` and :math:`R(2190) < 0.2` (unless
   ``allpixels=True``).
6. **Aerosol model.** Unless ``opac_model`` is given, the OPAC model :math:`m` whose spectral AOD is
   closest to the CAMS one is selected (the urban model is excluded):

   .. math::

      \hat m = \underset{m}{\operatorname{argmin}} \sum_\lambda
      \left| \frac{\overline{\tau_{CAMS}(\lambda)}}{\overline{\tau_{CAMS}(550)}} - \tau_m(\lambda; \tau_{ref} = 1) \right|

   where the bars are averages over the image.
7. **Look-up tables.** The radiative transfer LUTs are reduced to the geometries of the image and to
   the selected model (`Atmospheric diffuse light`_).
8. **Wind speed.** The wind speed of CAMS is used for the sunglint (the retrieval from the image,
   :py:meth:`~grs.grs_kernel.Kernel.rugosity_est_chunk`, is not used for the moment).
9. **Aerosol optical thickness** on a coarse grid (`Aerosol retrieval`_).
10. **Refined look-up tables** on the range of the retrieved AOT, with a 0.01 step.
11. **Atmospheric and sunglint correction** of every pixel, block by block
    (`Remote-sensing reflectance`_).
12. **Flags and mask.** Bit 18 (``neg_rrs``) is set where :math:`R_{rs}(490)` or :math:`R_{rs}(565)` is
    below :math:`-5 \times 10^{-4}`. The binary ``mask`` is 0 (good) for the pixels with bit 3 (water
    index) set and none of the bits 0 (no data), 1 (s2cloudless cloud), 10 (SURFWATER cloud and shadow),
    13 (OmniCloudMask shadow), 14 (Landsat dilated cloud) and 18 (negative Rrs); it is 1 elsewhere
    (:py:meth:`Masking.create_mask <grs.mask.Masking.create_mask>`).
13. **Output.** :py:class:`~grs.output.L2aProduct` writes ``Rrs``, ``BRDFg``, ``aot550``, ``wind``,
    ``flags``, ``mask`` and the ancillary data as Zarr or NetCDF (see :doc:`api/output`).

Gaseous transmittance
---------------------

For each gas, the band transmittance :math:`T_{gas}(\lambda; c)` is tabulated against the total
column :math:`c` and convolved with the spectral response of the bands. Step 4 corrects for the
gases absorbing high in the atmosphere, above the scattering layer:

.. math::

   R_{TOA}(\lambda) \leftarrow \frac{R_{TOA}(\lambda)}{T_{O_3}(\lambda)\, T_{NO_2}(\lambda)\, T_{bg}(\lambda)}

(:math:`T_{bg}`: CO, CO\ :sub:`2`, O\ :sub:`2` and O\ :sub:`4`, for the CAMS surface pressure). Water vapour and methane, mostly in the lower layers, are taken
into account in the radiative transfer equations below through

.. math::

   T_g(\lambda) = T_{H_2O}(\lambda; C_{H_2O})\, T_{CH_4}(\lambda; C_{CH_4}), \qquad
   T_g^{diff}(\lambda) = T_{H_2O}(\lambda; 0.5\, C_{H_2O})\, T_{CH_4}(\lambda; 0.65\, C_{CH_4})

where the factors 0.5 and 0.65 account for the part of the diffuse light path that is above the gases.

Atmospheric diffuse light
-------------------------

The LUT intensities :math:`I_m` of the selected aerosol model (weight :math:`w_m`) give the reflectance
of the atmosphere and of the sky light reflected by the air-water interface:

.. math::

   R_{diff}(\lambda; \tau_{ref}, \theta_s, \theta_v, \phi) =
   \frac{1}{\cos\theta_s} \sum_m w_m\, I_m(\lambda; \tau_{ref}, \theta_s, \theta_v, \phi)

which is removed from the TOA signal after scaling for the pressure and the gaseous absorption:

.. math::

   R_{corr}(\lambda) = R_{TOA}(\lambda) - T_g^{diff}(\lambda)\, \frac{P}{P_0}\, R_{diff}(\lambda; \tau_{ref})

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

   \sigma^2 = 0.003 + 0.00512\, U

The water-leaving signal is negligible in the SWIR: the sunglint is measured in a SWIR band
:math:`\lambda_s` and extrapolated to all the bands with its spectral shape :math:`\varepsilon(\lambda)`
and the slope probability of each band (which differs between the bands of Sentinel-2, because of
their different viewing angles):

.. math::

   R_f = \frac{\varepsilon(\lambda_r)\, P_{slope}(\lambda_r)\, R_{corr}(\lambda_s)}
              {T_g(\lambda_s)\, T_{dir}(\lambda_s)\, \varepsilon(\lambda_s)\, P_{slope}(\lambda_s)}

.. math::

   R_g(\lambda) = T_g(\lambda)\, T_{dir}(\lambda)\, \varepsilon(\lambda)\, P_{slope}(\lambda)\,
   \frac{R_f}{\varepsilon(\lambda_r)\, P_{slope}(\lambda_r)}

:math:`R_f` is exported as ``BRDFg``.

.. note::

   The slope probability :math:`P_{slope}` is currently set to 1 (see
   :py:meth:`~grs.grs_kernel.Kernel.p_slope`): the spectral extrapolation then only depends on
   :math:`\varepsilon(\lambda)` and on the transmittances.

Remote-sensing reflectance
--------------------------

.. math::

   R_{rs}(\lambda) = \frac{R_{corr}(\lambda) - R_g(\lambda)}{\pi\, T_{tot}(\lambda)}

(the division by :math:`T_{tot}` is only applied where :math:`R_{corr} \ge R_g`). Every pixel uses the
AOT interpolated from the coarse map of step 9, and its own pressure and gaseous transmittances.

.. _methods-aerosol:

Aerosol retrieval
-----------------

The AOT is retrieved on a coarse grid of pixel size ``aot_resolution`` (480 m by default): each coarse
pixel :math:`B` is the mean of its water pixels :math:`W_B`, and is kept if it has enough of them,

.. math::

   \bar R_B = \frac{1}{N_B} \sum_{\mathbf{x} \in W_B} R(\mathbf{x}), \qquad N_B > 0.5\, \max_{B'} N_{B'}

The coarse grid is cut into blocks of ``aot_chunk`` pixels. In each block, the AOT is the largest one
that does not make the reflectance of too many pixels negative. For a set of tested values,

.. math::

   N(\tau_{ref}) = \#\left\{ \mathbf{x} \in V : \min_\lambda R(\lambda, \mathbf{x}; \tau_{ref}) < 0 \right\},
   \qquad
   \hat\tau_{ref} = \max\left\{ \tau_{ref} : N(\tau_{ref}) < 0.05\, \max_\tau N(\tau) \right\}

where :math:`V` are the pixels that are not already negative for :math:`\tau_{ref} = 0` (e.g.
shadows). This is done in three passes:

1. **SWIR bound**, with :math:`R = R_{corr}` in the two SWIR bands
   (:py:meth:`~grs.grs_kernel.Kernel.aerosol_swir_chunk`);
2. the same, refined with 21 values between 0 and :math:`\max(\tau_{max}, 0.12)`;
3. **visible-NIR retrieval**, with :math:`R = (R_{corr} - R_g)/\pi` in the first seven bands, for 21
   values between 0 and :math:`1.02\, \tau_{max}` (:py:meth:`~grs.grs_kernel.Kernel.aerosol_visible_chunk`).
   Without a valid value, the mean CAMS AOD of the block is used.

The map of the block values is then cleaned and smoothed: values above 1.5 times the median are
replaced by the median, then a weighted mean filter (weights: number of water pixels of the blocks)
and a NaN-aware box filter are applied, and the map is interpolated back to the coarse grid.

The coarse resolution is stored in the ``aot_resolution`` attribute of the output product.

References
----------

* Harmel, T., Chami, M., Tormos, T., Reynaud, N., Danis, P.-A. (2018). Sunglint correction of the
  Multi-Spectral Instrument (MSI)-SENTINEL-2 imagery over inland and sea waters from SWIR bands.
  *Remote Sensing of Environment*, 204, 308-321.
  `doi:10.1016/j.rse.2017.10.022 <https://doi.org/10.1016/j.rse.2017.10.022>`__
* Cox, C., Munk, W. (1954). Measurement of the roughness of the sea surface from photographs of the
  sun's glitter. *Journal of the Optical Society of America*, 44(11), 838-850.
