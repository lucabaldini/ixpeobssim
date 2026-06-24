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

   All the background components treated by our software are assumed to be completely
   unpolarized. Some instrumental background components may be polarized, but over the
   years no way has been found to predict its polarization angle and degree so the only
   way to asses its polarization is to measure it directly from the data, and then
   subtract it from the total signal.

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
   (and, possibly, developing an interface to the e-ROSITA data) is more than
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
* The current implementation of the X-ray instrumental background is based on the
  average value of counts in the full IXPE band obtained from years of
  data taken in occultation (e.g.: when the source is blocked by the Earth).
  In addition an in-eclipse (e.g..: when the sun is not shining on the detector)
  filter and a particle cut (Di Marco, et al. AJ 165, 143 (2023)) are applied.
* To our knowledge, this background component is stable over
  time and shows no significant polarization, although a mild trend of counting
  rate decrease over multiple years may be expected. It is always a good idea to
  check that the static background component does not overshoot the data.

By default the instrumental background is generated uniformly in instrument
coordinates, although we do support a generic linear radial dependence,
as briefly summarized in the next section.


Radial dependence
~~~~~~~~~~~~~~~~~

IXPE celestial observations show that the distribution of the instrumental
background in detector coordinates is not necessarily rigorously flat. In
this case, we provide the means to model it with a linear function of the
radial distance from the center.

The radial slope alpha represents the fractional half-excursion of the variation
across the size h of the fiducial rectangle. For alpha = 0 the detector position
are distributed uniformly over the fiducial rectangle. For alpha = 2 the radial
dependence is maximal, and the density of events is zero at the center of the
detector.

The values of alpha found experimentally vary from observation to observation
(at a level which is consistent with expectations based on the typical variations
of the radiation environment in low-Earth orbit), and are typically lower than
0.2.

Flaring background
~~~~~~~~~~~~~~~~~

In addition to the static background, IXPE observations are affected by flaring in
time coincidence with solar flares when the sun is shining on the detector. The time
scale of this effect is of the order of minutes to a few hours and can be characterized
by the difference between the spectrum obtained when the sun is shining on the detector
and the spectrum obtained when the sun is not shining on the detector.
There are currently two ways of getting rid of the flaring background

* The flares can then be estimated in the full detector area and subtracted from the data.
  This relies on the assumption that the flaring activity is uniform across the detector,
  and recovers a relatively high statistics of the flare by integrating over it. It can then
  be subtracted from any spatial selection (or map binning) so that both the spectrum and the
  polarization properties of the flare are properly accounted for. This is handled by
  :ref:`reference-xpbin` and :ref:`reference-xpsun`.
* :ref:`reference-xpsun` splits the observation dataset in two parts, the "insun", where
  the sun is shining on the detector, and the "ineclipse", where the sun is not shining on the
  detector. The two datasets can be used to create a "deflared" dataset by subtraction
* :ref:`reference-xpbin` launched with the level 2 files with the appropriate spatial (not time!)
  cuts as positional argument, the full-detector insun and ineclipse with the --insun and
  --ineclipse options will provide a _deflared.fits file which has been already flare-subtracted.
* The flares can be removed with an ad-hoc cut in the time intervals, removing the time
  intervals in which the count rate increases over some arbitrary threshold. A standard
  way of doing this is the so-called "sigma-clipping" method with a 3 sigma threshold, or
  multiple pass sigma clipping. Unfortunately, the solar activity has been proven to be
  often sub-threshold and carries a small number of counts that sometimes are however
  persistently strongly polarized (Bucciantini et al.  A&A 699, A33 (2025)). The software
  does not currently support this method.


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
can be used as an argument to generate an :class:`ixpeobssim.srcmodel.bkg.xTemplateInstrumentalBkg`
object into a config file to be used by |xpobssim|.
Note that while this is designed to create the average spectrum from a single
observation, the algorithm will run just fine on multiple observation and treat
each file in the same way.

The relevant piece of code for the generation of the template resides in
:mod:`ixpeobssim.bkg.instr` in which the function create_background_template
is defined and invoked by |xpbkgtemplate|.
