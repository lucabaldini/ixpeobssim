.. _quickstart:

Quick start
===========


The main purpose of this simulation package is to simulate an observation of a
given source model, based on suitable detector response functions.

xppimms
-------

As a first test case we can evaluate the MDP for a single point source using the
|xppimms| tool. Its capabilities and options are described in the
:doc:`reference` section of the documentation.

In order to reproduce the results reported
`here <http://bigfoot.iaps.inaf.it:8080/xwiki/wiki/ixpeglobal/view/Main/IXPEsensitivityFiles/?srid=TfFqiP8o>`_
for a standard source (power law spectrum with index 2 and normalization 10 at
1 keV, no absorption), in 100 ks and in the 2-8 keV energy range, we run the
command

.. code-block:: shell

    xppimms --duration 10000

Which prints out the following MDP table:

.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xppimms.py --duration 10000
   :ellipsis: 0,-5
   :shell:

This program calculates the MDP by direct numerical integration of the power-law
input spectrum. As a consequence, part of the richness of the detector response
(most notably, the energy dispersion and the effective area vignetting) is not
captured here. For use cases beyond simple stationary point sources, the use
of |xpobssim| and |xpbin| in PCUBE mode are recommended, as this approach offers
the maximum flexibility (see next section for an example).


xpobssim
--------

The main Monte Carlo simulation application is |xpobssim|. For a quick reference
on this tool see :doc:`reference` section. Assuming that the current working
directory is the ixpeobssim root folder, the command

.. code-block:: shell

		xpobssim --configfile ixpeobssim/config/toy_point_source.py --duration 100000

Will produce three event (FITS) files (one for each Detector Unit) with a 100 ks
simulation of a stationary source with a power-law spectrum (with an index of 2
and normalization of 10) with energy- and time-independent polarization degree
(0.5) and angle (30°), correctly folded with all the instrument response
functions: effective area, modulation factor, energy dispersion and point-spread
function.

In order to bin the event files in energy (from 2 keV to 8 keV in one single bin)
we run the command

.. code-block:: shell

		xpbin $IXPEOBSSIM_DATA$/toy_point_source_du?.fits --algorithm PCUBE --emin 2. --emax 8. --ebins 1

that will produce three new FITS files (called modulation cubes), containing the
MDP value (at 99%) and the histogram of the azimuthal distribution of
photoelectrons emission for every bin.

The binned output files can be easily visualized using the |xpbinview| tool:

.. code-block:: shell

		xpbinview $IXPEOBSSIM_DATA$/toy_point_source_du?_pcube.fits

We are already fully equipped for a basic spectral analysis with XSPEC.
The first step is to bin again the event files by running the |xpbin| tool with
the PHA1 algorithm.

.. code-block:: shell

		xpbin $IXPEOBSSIM_DATA$/toy_point_source_du?.fits --algorithm PHA1


We can feed the binned files (along with the corresponding .arf and .rmf
response functions) into XSPEC and recover the input parameters of our source.

.. code-block:: shell

		xpxspec $IXPEOBSSIM_DATA$/toy_point_source_du?_pha1.fits
