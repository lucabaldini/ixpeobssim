.. _irf:

Response functions
==================

The instrument response functions (IRFs) are a critical part of ixpeobssim, and
they are used (in identical form) both in the event simulation and in the
analysis of the data products from the simulation itself.

All the response functions are stored in FITS files in the OGIP format defined
in `CAL/GEN/92-002 <https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html>`_ (and modifications
`CAL/GEN/92-002a <https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002a/cal_gen_92_002a.html>`_) and are intended to be fully compatible
with the spectral analysis tools provided by
`XSPEC <https://heasarc.gsfc.nasa.gov/xanadu/xspec>`_ (see
`this page <https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/fits/fitsfiles.html>`_
for more details).

We identify six different types of response functions:

* effective area (.arf), see :class:`ixpeobssim.irf.arf.xEffectiveArea`;
* vignetting, see :class:`ixpeobssim.irf.vign.xVignetting`;
* energy dispersion (.rmf), see :class:`ixpeobssim.irf.rmf.xEnergyDispersion`;
* point-spread function, see :class:`ixpeobssim.irf.psf.xPointSpreadFunction`;
* modulation factor, see :class:`ixpeobssim.irf.modf.xModulationFactor`;
* modulation response function, see :class:`ixpeobssim.irf.mrf.xModulationResponse`.

If you are familiar with basic spectral analysis in XSPEC, the .arf and
.rmf files have exactly the meaning that you would expect, and can be in fact
used in XSPEC; the .mrf files have the exact same format as the .arf files.

ixpeobssim provides facilities for generating, reading displaying and using
IRFs, as illustrated below.


IRF summary
-----------

This is a short, top level description of the version v11 of the default IXPE
response functions used by ixpeobssim. If you want to dig into all the gory details
of the process of generating the response files, the section about :ref:`irfgen`
is the place you want to start from.

All the response functions are defined between 1 and 12 keV in steps of 40 eV,
and are detector-unit based (the pulse-invariant channel space is defined
between 0 and 15 keV in steps of 40 eV). As such, there are differences between
the three modules in terms of effective area, modulation factor and
point-spread function.

Since the eight iteration of the response functions, all response files come
in two distinct flavors---un-weighted and weighted. Accordingly, all the analysis
tools support weights.

.. warning::

   At this point in test the analysis of IXPE data with weights has not been
   thoroughly tested. While weight do allow for a measurable sensitivity
   improvement, cross-checking your analysis using the plain, old, un-weighted
   strategy is a sensible way to make sure everything is in order.


Effective area
~~~~~~~~~~~~~~

Below is the nominal IXPE on-axis effective area, in both the un-weighted and
weighted flavors. The latter is roughly 15% smaller at the peak value but, as
we shall see in a second, the increase in the modulation factor drives the
overall sensitivity slightly up.

.. _figure-ixpe_effective_area:
.. figure:: figures/irf/ixpe_effective_area_v11.*
   :width: 80%

   On-axis effective area as a function of the energy. The solid lines represent
   the total effective area for the sum of three telescopes, in the un-weighted and
   weighted version, while the un-labeled dashed lines represent the curve for each
   of the three single telescopes.

The effective-area curves for the three telescopes are within a few % from each
other, the small differences being due to the slightly different mirror effective
areas measured during the MMA calibration, as well as the different asymptotic
pressure values for the three GPD at the focal plane.

The effective-area calculation in ixpeobssim includes all the relevant
contributions, namely:

* the mirror effective area;
* the transparency of mirror-module-assembly thermal shield;
* the transparency of the detector-unit UV filter;
* the transparency of the GPD window;
* the efficiency of the GPD gaseous active medium;
* the efficiency of the event weighting (in the weighted flavor).

The plots below show the principal ingredients that go into the calculation.

.. _figure-mma_effective_area:
.. figure:: figures/irf/mma_effective_area_v11.*
   :width: 80%

   On-axis effective area as a function of the energy for three Mirror-Module
   Assemblies (MMA) and for a single module.

.. _figure-uv_filter_transparency:
.. figure:: figures/irf/uv_filter_transparency_v11.*
   :width: 80%

   Transparency of the UV filter as a function of the photon energy.

.. _figure-gpd_quantum_efficiency:
.. figure:: figures/irf/gpd_quantum_efficiency_v11.*
   :width: 80%

   Quantum efficiency of the GPD as a function of the energy, broken up in its
   two main components---the Be window transparency, and the gas cell absorbing
   efficiency.


The vignetting function shown below comes from a preliminary study by MSFC
based upon ray-trace simulations for a perfect mirror module assembly, and is
relevant for the simulation of extended sources.

.. _figure-mma_vignetting:
.. figure:: figures/irf/mma_vignetting_v11.*
   :width: 80%

   Preliminary estimation of the vignetting of the optics as a function of
   energy and off-axis angle.

The vignetting, along with the relative orientation of the three IXPE detector
units, defines the relative exposure across the field of view of the instrument,
as shown in the following two plots. (Note that above 6 keV the drop of the
effective area at the edge of the field of view is relatively more important.)

.. _figure-field_of_view_at_3_kev:
.. figure:: figures/irf/field_of_view_at_3_kev_v11.*
   :width: 80%

   Relative exposure at 3 keV across the field of view for the set of three
   telescopes clocked in the IXPE configuration.

.. _figure-field_of_view_at_8_kev:
.. figure:: figures/irf/field_of_view_at_8_kev_v11.*
   :width: 80%

   Relative exposure at 8 keV across the field of view for the set of three
   telescopes clocked in the IXPE configuration.

Of course, unless you specifically decide to disable this functionality,
ixpeobssim handles all of this behind the scene, so you don't have to worry
about it---but keep it in mind when you do back-of-the-envelope
calculations.



Energy dispersion
~~~~~~~~~~~~~~~~~

The energy dispersion (a.k.a. the response matrix) comes from a series of line
Monte Carlo simulations performed with the IXPE GPD Geant 4 simulation framework.
Below is a color representation of the energy dispersion as a function of the energy,
which is essentially the content of the binary table in the ``MATRIX`` extension of
the rmf file.

.. _figure-energy_dispersion:
.. figure:: figures/irf/energy_dispersion_v11.*
   :width: 80%

   Representation of the GPD response matrix.

For illustration purposes, here are the corresponding one-dimensional pdfs
at a few fixed true energies (i.e., these are just vertical slices of the
color plot above).

.. _figure-energy_resolution:
.. figure:: figures/irf/energy_resolution_v11.*
   :width: 80%

   Energy dispersion (one-dimensional probability density function) at a set
   of discrete energies. The FWHM energy resolution is indicated for
   completeness.


Point-spread function
~~~~~~~~~~~~~~~~~~~~~

The PSF model is derived from one of the early point-source observations, as
described in `issue #158 <https://github.com/lucabaldini/ixpeobssim/issues/158>`_.

.. note::
  Starting with version 6 of the instrument response function each DU comes
  with a different PSF scaling factor to account for the differences
  measured during the mirror calibration. As shown in
  `issue #387 <https://github.com/lucabaldini/ixpeobssim/issues/387>`_, MMA 1
  has a significantly better PSF (less than 20 arcsec HPD) than MMAs 2 and 3
  (running at more than 25 arcsec HPD).


.. _figure-psf_eef:
.. figure:: figures/irf/psf_eef_v11.*
   :width: 80%

   Encircled energy fraction (EEF) for the PSF of the three IXPE telescopes.

For completeness, the current set of response functions do not include the
mirror aberration, which is nonetheless much smaller than the PSF half-power
diameter across the entire field of view and is therefore, to first order,
negligible.


Modulation factor
~~~~~~~~~~~~~~~~~

Our parametrization of the modulation factor comes from a series of line-type
Monte Carlo simulations, informed by the ground calibrations of the three
detector units.

.. _figure-ixpe_modulation_factor:
.. figure:: figures/irf/ixpe_modulation_factor_v11.*
   :width: 80%

   Modulation factor as a function of the photon energy for the IXPE detectors.
   The solid line represents the average for the three GPD, in the un-weighted
   and weighted version, while the un-labeled dashed lines (admittedly, barely
   visible) represent the curve for each of the three detectors.

.. warning::
   The noticeable edge around 9 keV is due to the K-edge of the copper,
   above which the extraction of photoelectrons from X-rays absorbed in the GEM
   becomes significantly more likely. This causes an increase of effective area,
   accompanied by a dilution of the modulation. While we provide a tabulation of
   all the IRFs in the standard grid between 1 and 12 keV, significant more
   work is needed to validate the response of the detector above the Cu
   K-edge, and simulations outside the 2--8 keV standard range should be
   interpreted with caution.


Minimum detectable polarization
-------------------------------

The effective area curve (for the sum of three mirror modules) and the
modulation factor are enough for a crude estimation of the minimum detectable
polarization for a point source, and for reference we produce the basic
performance plot below for each iteration of the response functions using
`xppimms` for definite sets of spectral indices and rescaling for the
source flux and the observing time.

.. _figure-ixpe_mdp:
.. figure:: figures/irf/ixpe_mdp_v11.*
   :width: 80%

   IXPE Minimum Detectable Polarization (MDP) as a function of the source flux
   for several different exposure times and spectral indices.


Reading and visualizing IRFs
----------------------------

In a nutshell, the recommended way to load the default set of response
functions (whatever that means at any point in time) is

.. code-block:: python

    from ixpeobssim.irf import load_irf_set

    # Load all the default response functions.
    irf_set = load_irf_set(du_id=1)

    # Access the actual response functions.
    aeff = irf_set.aeff
    vign = aeff.vignetting
    edisp = irf_set.edisp
    psf = irf_set.psf
    modf = irf.modf

    # Print the effective area and modulation factor at 5 keV.
    print(aeff(5.))
    print(modf(5.))

The reader is referred to the documentation and the source code of the relevant
classes for a full description of the interfaces that ixpeobssim provides.

For completeness, ixpeobssim makes available `xpirfview.py` as a single
visualization interface to all the response file. Just type

.. code-block:: shell

   xpirfview.py path/to/the/response/file.fits

and you should get back some sensible visualization of the thing.


Pseudo-CALDB
------------

For convenience, at this point in time, ``ixpeobssim`` is effectively implementing
its own, self-contained CALDB---that we sometimes refer to as the ``ixpeobssim``
`pseudo-CALDB`. The plans for interfacing ``ixpeobssim`` with the actual IXPE CALDB
are not yet defined as there are definitely peculiarities on both sides
(simulation and real data) that make having a drop-in replacement structure less
than trivial.

The actual FITS files with the response data live (provisionally) in the
`ixpeobssim/caldb <https://github.com/lucabaldini/ixpeobssim/tree/main/ixpeobssim/caldb/ixpe>`_
folder and the basic logic determining the naming and the file location is defined
`here <https://github.com/lucabaldini/ixpeobssim/blob/main/ixpeobssim/irf/__init__.py>`_.

.. note::

  .. versionadded:: 21.0.0

  Starting from version 21.0.0 the structure of the pseudo-CALDB has been
  drastically changed to match as closely as possible that of the actual
  CALDB submitted to HEASARC.

  The ``ixpeobssim`` internal rules for the IRF-name designation have also been
  modified changing the delimiter, in order to have a better match with the
  CALDB file names and provide support for weights in a more straightforward
  fashion.

  ``ixpeobssim`` will, at least provisionally, maintain a separate version
  numbering with respect to the official CALDB.

The convention we use to name response file is
``[base][unit][calibtype][intent][ver]``, where:

* [base] is the base name for the set of response functions, e.g., `ixpe_`;
* [unit] indicates the telescope unit (`du1`, `du2`, `du3`).
* [calibtype] provides an identifier for the calibration data (e.g., `vign` or
  `psf`), with the exception of the `arf`, `rmf` and `mrf` files, where the
  data type is indicated by the file extension;
* [intent] is the intent of a particular set of response functions, e.g.,
  `_obssim_`;
* [version] is the CALDB version number for any given file.

Additionally, each coherent set of response functions is identified within
ixpeobssim by a name, in the form ``[base]:[intent]:[ver]``. ixpeobssim
is able to parse a string formed according to this rule and resolve all the
relevant paths to the actual response files.

.. warning::

   The pseudo-CALDB contains a number of response files that are not shipped
   with the real CALDB, including pre-launch estimates that we retain for
   bookkeeping purposes, but should not be used to analyze flight data.

   In a nutshell: all the set of response files named as ``*_legacy_*``
   should be considered of solely historical interest and should never be
   used in conjunction with flight data samples.


Response file versioning
~~~~~~~~~~~~~~~~~~~~~~~~

This is a short description of the main features of different sets of response
files that are useful for simulation and science analysis:

* ``ixpe:obssim:v11``: identical to ``ixpe:obssim:v10``, except that the PSF
  parametrization (used on the simulation side of things) has been improved to
  match the on-orbit radial dependence measured with point sources.
* ``ixpe:obssim:v10``: this is the first iteration of the response files
  matching the structure of the actual CALDB, and the first that can be used
  with real data.


Mapping to the real CALDB
~~~~~~~~~~~~~~~~~~~~~~~~~

Although the pseudo-CALDB and the real-CALDB are fundamentally different in
some key aspects, most of the relevant files (e.g., those containing the
effective area, the response matrix and the modulation response function)
have a definite, one-to-one correspondence between the two databases---meaning
that they are `identical`, modulo a few header keywords.

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - Pseudo CALDB
     - Real CALDB
   * - gpd/cpf/arf/ixpe_d?_obssim_v010/011.arf
     - gpd/cpf/arf/ixpe_d?_20170101_01.arf
   * - gpd/cpf/arf/ixpe_d?_obssim_alpha075_v010/011.arf
     - gpd/cpf/arf/ixpe_d?_20170101_alpha075_01.arf
   * - gpd/cpf/rmf/ixpe_d?_obssim_v010/011.rmf
     - gpd/cpf/rmf/ixpe_d?_20170101_01.rmf
   * - gpd/cpf/rmf/ixpe_d?_obssim_alpha075_v010/011.rmf
     - gpd/cpf/rmf/ixpe_d?_20170101_alpha075_01.rmf
   * - gpd/cpf/mrf/ixpe_d?_obssim_v010/011.mrf
     - gpd/cpf/mrf/ixpe_d?_20170101_01.mrf
   * - gpd/cpf/mrf/ixpe_d?_obssim_alpha075_v010/011.mrf
     - gpd/cpf/mrf/ixpe_d?_20170101_alpha075_01.mrf
   * - gpd/cpf/modfact/ixpe_d?_obssim_mfact_v010/011.fits
     - gpd/cpf/modfact/ixpe_d?_20170101_mfact_01.fits
   * - gpd/cpf/modfact/ixpe_d?_obssim_alpha075_mfact_v010/011.fits
     - gpd/cpf/modfact/ixpe_d?_20170101_mfact_alpha075_01.fits



Differences with the real CALDB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The structure of the ``ixpeobssim`` pseudo-CALDB is designed to match as closely
as possible that of the real IXPE CALDB.

.. code-block::

  ixpe
  |----gpd
      |----bcf
           |----chrgparams
      |---- cpf
           |----arf
           |----modfact
           |----mrf
           |----rmf
  |----mma
      |----bcf
          |----psf
          |----vign


The naming conventions for the FITS files have been aligned to the IXPE
CALDB starting from version 10, and the file names for the previous iterations
have been changed accordingly.

There are two subtle but noticeable differences between the IXPE CALDB and the
``ixpeobssim`` pseudo-CALDB, namely:

* the pseudo-CALDB has no concept of validity date, and that is reflected both
  in the file names and in the corresponding header keywords;
* the versioning scheme has a different meaning in the two contexts, and
  version numbers are physically assigned by different people at different times
  (more specifically, ``ixpeobssim`` uses a unique, sequential identifier that
  is tied to the basic ingredients going into the response functions, while
  the actual IXPE CALDB can have the same identifiers for files with a different
  validity epoch).

The latter difference is also reflected in the format string for the version
identifier in the file name, which is ``'%02d'`` for the IXPE CALDB and
``'v%03d'`` in the pseudo-CALDB.

If you are really careful you will also notice that, for weighted response files,
the position of the weight identifier in the file name `for the modulation factor`
is different in the pseudo-CALDB, compared with the real one, i.e., the mapping
is in this case

.. code-block::

   ixpe/gpd/cpf/modfact/ixpe_d1_obssim_alpha075_mfact_v010.fits
   ixpe/gpd/cpf/modfact/ixpe_d1_20170101_mfact_alpha075_01.fits

This is just an historical accident that, at this point, is not worth correcting.
(And you probably will never need the modulation factor, anyway.)


Historical notes
----------------

The release process and the differences with respect to the previous iterations
are summarized on our issue tracker at:

* `issue #580 <https://github.com/lucabaldini/ixpeobssim/issues/580>`_
  (release of version 11);
* `issue #496 <https://github.com/lucabaldini/ixpeobssim/issues/496>`_
  (release of version 10);
* `issue #460 <https://github.com/lucabaldini/ixpeobssim/issues/460>`_
  (release of version 9);
* `issue #402 <https://github.com/lucabaldini/ixpeobssim/issues/402>`_
  (release of version 7);
* `issue #333 <https://github.com/lucabaldini/ixpeobssim/issues/333>`_
  (release of version 6);
* `issue #344 <https://github.com/lucabaldini/ixpeobssim/issues/344>`_
  (release of version 5);
* `issue #294 <https://github.com/lucabaldini/ixpeobssim/issues/294>`_
  (release of version 4);
* `issue #258 <https://github.com/lucabaldini/ixpeobssim/issues/258>`_
  (release of version 3);
* `issue #161 <https://github.com/lucabaldini/ixpeobssim/issues/161>`_
  (release of version 2 and differences with respect to version 1).

In iterations of the response functions up to v3, ixpeobssim used
to ship combined versions of the effective area and modulation factor, that
were useful for back-of-the envelope sensitivity calculations. From
version 4 onward this is no more the case, and all the relevant applications
have been modified to make the appropriate loop over the three detector units
where the combined response functions were used before.

Likewise, all the non-standard versions of the response files (e.g., without
the standard cuts or with the MMA alone) have been dropped altogether
starting from version 4.

Version 6 of the response function is the first iteration taking full
advantage of the flight DU calibration and the telescope end-to-end
calibration. This is also the last iteration of the response function
using the now infamous 80% cut---the next ones will hopefully support
ensemble-weighted analyses.

Starting with version 6 of the instrument response function the PI runs
from 0 to 374 (included), corresponding to a physical-space binning
spanning the 0--15 keV energy range in steps of 40 eV.
(In previous iterations the PI spanned the very same energy interval used
to define the response functions, i.e., 1--12 keV.)

Version 7 of the response files features the first non-diagonal response matrix.

Version 8 of the response files is the first supporting XSPEC spectro-polarimetric
analysis with weights. Version 9 is fairly similar, with a refined parametrization
of the MMA effective area. Version 10 features a few new header keywords, and
is the one on which the first version of the CALDB submitted to HEASARC is
based.
