.. _background:

Backgrounds
===========

Although for relatively bright point sources background is typically not an issue
for IXPE, it might become important for dimmer point sources and/or extended
sources, especially at high energy.

We distinguish three different types of backgrounds: galactic, extra-galactic, and
instrumental---as it turns out, the latter is typically the most important.
ixpeobssim provides ways to simulate background components and handle backgrounds
in high-level analysis. In the remaining of this section we shall briefly cover
the simulation part.

.. note::

   All the background components are assumed to be completely unpolarized, which
   is consistent with our best understanding of the detector response.


.. seealso::

   The :mod:`ixpeobssim.srcmodel.bkg` module contains all the background-related
   modeling facilities.



Galactic background
-------------------

This is, admittedly, the part that is less extensively covered in ixpeobssim.
We do have a :class:`ixpeobssim.srcmodel.bkg.xGalacticBkg` class that is designed
to interface to the the
`HEASARC background tool <https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xraybg/xraybg.pl>`_
and, particularly, use the ROSAT count rate in the 1--2 keV energy band (R7).
The model is overly simple (uniform in sky coordinates over the field of view
and with a fixed spectral index), with the overall intensity being the only
parameter that can be changed by the user.

.. note::

   Tuning the galactic X-ray background to a specific observation is still tricky
   in the current ixpeobssim setup. Any help in making this more user-friendly
   (and, possibly, developing an iterface to the e-ROSITA data) is more than
   welcome.


Extra-galactic background
-------------------------

The extra-galactic background, although not really the relevant contribution in
any practical situation, is relatively straightforward to handle. We provide a
:class:`ixpeobssim.srcmodel.bkg.xExtragalacticBkg`, with the basic parametrization
taken from `Gruber et al., 1999 <https://iopscience.iop.org/article/10.1086/307450/pdf>`_.



Instrumental background
-----------------------

This is most often the primary source of background, its main characteristics
being that the instrumental background is modeled in instrument coordinates and
not convolved with the instrument response functions---not with the effective area,
nor with the vignetting, and typically, not with the energy dispersion.

In terms of simulation facilities, the two most important data structures are

* :class:`ixpeobssim.srcmodel.bkg.xPowerLawInstrumentalBkg`: a background
  component with a power-law energy spectrum with adjustable normalization and
  spectral index;
* :class:`ixpeobssim.srcmodel.bkg.xTemplateInstrumentalBkg`: a background
  component with an energy spectrum derived from an actual observation---typically
  by using a dim point source and removing the source counts at the center of the
  field of view, and then properly rescaling the the relevant area.

By deafault the instrumental background is generated uniformly in instrument
coordinates, although, based on actual observations, we do support a generic
linear radial dependence, as briefly summarized in the next section.


Radial dependence
~~~~~~~~~~~~~~~~~

To be filled in.


Creating a template
~~~~~~~~~~~~~~~~~~~

To be filled in.
