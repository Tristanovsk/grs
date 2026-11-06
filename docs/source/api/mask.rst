grs.mask
========

.. automodule:: grs.mask
   :no-members:

:py:meth:`Masking.process <grs.mask.Masking.process>` builds the ``flags`` variable, a bit field where
each bit is one test. The names and descriptions of the bits (with the thresholds used) are stored in the
``flag_names`` and ``flag_descriptions`` attributes of ``flags``, and the fraction of the image raised
for each flag in the global attributes ``flag_<name>``.

.. list-table::
   :header-rows: 1
   :widths: 8 27 65

   * - Bit
     - Name
     - Description
   * - 0
     - ``nodata``
     - no data in the input image
   * - 1
     - ``cloud_p06``
     - low confidence cloud mask from s2cloudless (probability > 0.6, dilation 12), Sentinel-2
   * - 2
     - ``cloud_p08``
     - high confidence cloud mask from s2cloudless (probability > 0.8, dilation 5), Sentinel-2
   * - 3
     - ``water_swir_visible_index``
     - water from the normalized index between the 490 and 1610 nm bands (> 0.2)
   * - 4
     - ``water_red_visible_index``
     - water from the normalized index between the 490 and 842 nm bands (> 0); can fail for turbid
       waters
   * - 5
     - ``thin_cirrus``
     - cirrus band reflectance > 0.003
   * - 6
     - ``opac_cirrus``
     - cirrus band reflectance > 0.009
   * - 7
     - ``high_swir``
     - reflectance at 2190 nm > 0.2 (bright clouds, too bright reflection...)
   * - 8
     - ``surfwater_land``
     - land from the SURFWATER file
   * - 9
     - ``surfwater_water``
     - water from the SURFWATER file
   * - 10
     - ``surfwater_cloud_and_shadow``
     - cloud and shadow from the SURFWATER file
   * - 11 - 17
     - ``l1_dilated_cloud``, ``l1_cirrus``, ``l1_cloud``, ``l1_cloud_shadow``, ``l1_snow``,
       ``l1_clear``, ``l1_water``
     - flags copied from the QA band of Landsat 8/9
   * - 18
     - ``neg_rrs``
     - negative :math:`R_{rs}` (< -0.0005) at 490 or 565 nm, set after the correction

The binary ``mask`` of the output is made by
:py:meth:`Masking.create_mask <grs.mask.Masking.create_mask>` from the bits to keep (3) and the bits
to mask (0, 1, 10, 13, 14, 18): 0 for the good water pixels, 1 elsewhere.

To select a pixel on one flag:

.. code-block:: python

   cirrus = (l2a.flags & (1 << 5)) > 0
