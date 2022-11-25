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

IXPE celestial observations show that the distribution of the instrumental
background in detector coordinates is not rigorously flat---in fact it is typically
higher close to the walls of the detector---and can be reasonably modeled with a
linear function of the radial distance from the center.

The radial slope alpha represents the fractional half-excursion of the variation
across the size h of the fiducial rectangle. For alpha = 0 the detector position
are distributed uniformly over the fiducial rectangle. For alpha = 2 the radial
dependence is maximal, and the density of events is zero at the center of the
detector.

The values of alpha found experimentally vary from observation to observation
(at a level which is consistent with expectations based on the typical variations
of the radiation environment in low-Earth orbit), and are typically between
0 and 0.2.


Creating a template
~~~~~~~~~~~~~~~~~~~

:ref:`reference-xpbkgtemplate` provides a tool to generate an instrumental background starting
from an input pha1 file. The tool is designed specifically to work with dark fields
extracted from ixpe observation, and as such it requires the pha1 file to have
a backscal keyword defined in its header to account for the size of the extraction
region.

The rate of the background is then extracted from each file that is provided,
multiplied by its livetime and normalized with respect to the total detector area.
Everything is then averaged and rescaled back to physical units to create an
average spectrum which is smoothed with a spline and saved as an ascii file that
can be used as an agument to generate an :class:`ixpeobssim.srcmodel.bkg.xTemplateInstrumentalBkg`
object into a config file to be used by |xpobssim|.
Note that while this is designed to create the average spectrum from a single
observation, the algorithm will run just fine on multiple observation and treat
each file in the same way.

The relevant piece of code for the generation of the templateresides in
:mod:`ixpeobssim.bkg.instr` in which the function create_background_template
is defined and invoked by |xpbkgtemplate|.
