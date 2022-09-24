.. _ixpesim:

ixpesim interface
=================

If you have never heard about it, ``ixpesim`` (a C++ application living in the
`gdpsw <https://bitbucket.org/ixpesw/gpdsw>`_ repository) is the IXPE full detector
simulation based on the `Geant4 <https://geant4.web.cern.ch>`_ framework.
The process described in :ref:`irfgen` owns much of the external, auxiliary
information on (long) ``ixpesim`` simulations.

.. seealso::

   :ref:`irfgen`

``ixpesim`` is fundamentally different from ``ixpeobssim``, in that it follows the
impinging photons one at a time, simulating the interaction with the detector
at a microscopic level. As such, it provides a much richer information (e.g., real
track images, with all the associated topological information), well beyond
the ultimate capabilities of ``ixpeobssim``, due to the fact that the latter
only operates on a top-level description of the instrument (the response functions).
It goes without saying, this does not come at no cost, and generating events with
``ixpesim`` is about 1000 times slower than the fast-simulation counterpart---while
``ixpeobssim`` can throw a million events in a few seconds, it takes an hour or
so to ``ixpesim`` to achieve the same thing.

.. note::

   .. versionadded:: 16.13.0

      Over the years the inter-operation between ``ixpeobssim`` and ``ixpesim`` has
      improved at the point that it is now possible to simulate a simple celestial
      point source in ``ixpesim`` and feed the output back in to XSPEC for a
      fully-fledged spectro-polarimetric fitting.

   .. versionadded:: 19.0.0

      This is the first version with partial support of the concept of photon
      lists, i.e., list of photon energies, times and coordinates at the top of the
      Be window that can be fed into ``ixpesim`` to simulate arbitrary source
      models. Note this new functionality requires ``ixpesim`` version 13.6.0 or
      later.

   .. versionadded:: 21.3.0

      This version saw a complete refactoring of the photon list mechanism, and
      it is the first in which the latter can be considered fully functional,
      including the projection of the photon direction in the sky, with
      vignetting and dithering fully taken into account. This version
      requires ``ixpesim`` version 13.10.0 or later.


The basic work-flow for the inter-operation of ``ixpeobssim`` and ``ixpesim``
reads basically as follows:

* create a custom spectrum or, even better, a photon list to be fed into ``ixpesim``;
* run a simulation within ``ixpesim``;
* reconstruct the simulation with ``ixperecon``;
* process the ``ixpesim`` output with :ref:`reference-xpsimfmt` to make it
  compatible with XSPEC;
* feed all the ingredients in XSPEC for the spectro-polarimetic fit.

The various steps are further described in the remaining of this section.


Custom spectra
--------------

``ixpesim`` was mainly developed to support the characterization, test, and
calibration of the IXPE detector units, and, as such, is limited in scope to
the simulation of the gas pixel detectors at the focal plane of the
instrument---the simulation flow starts `right before the Berillium entrance window`.

Interestingly enough, ``ixpesim`` supports the use of custom input spectra in the
form of simple, two-column text files (one line per energy, the overall normalization
being irrelevant). This allows to incorporate the effect of all the elements
upstream the GPD (i.e., the mirrors and the UV filter) in such custom spectra.

.. warning::

   Admittedly, the flexibility of ``ixpesim`` in simulating celestial sources
   is limited---there are no real handles for the source morphology beyond a
   Gaussian beam, and the program is limited to energy-independent polarization
   degree and angle. This is, nonetheless, enough to perform a spectro-polarimetric
   fit within XSPEC. For more advanced applications, the use of photon lists
   is the obvious way to go.

``ixpeobssim`` provides the capability of creating arbitrary custom spectra to be
fed into ``ixpesim`` for simulating power-law sources, including the effect
of the galactic absorption with an arbitrary H column density, through the
:ref:`reference-xpsimspec` application.

:ref:`reference-xpsimspec` creates the convolution of:

* the intrinsic source spectrum;
* the galactic absorption;
* the mirror effective area;
* the transparency of the UV filter;
* the effect of the contaminants in the Be windows---which are accounted for
  in the response functions shipped with ``ixpeobssim``, but are not
  simluated in ``ixpesim``.

The output is a text file that can be used directly as a custom spectrum in
``ixpesim``.


Photon lists
------------

Since version 19.0.0 (but you are strongly encouraged to upgrade to version
21.3.0 or later, if you are playing this game) ``ixpeobssim`` provides the
:ref:`reference-xpphotonlist` application, that allows to generate a list of photon
definitions (in the form of a FITS file) that can then be fed into ``ixpesim``
to be propagated through the detector and turned into a list of actual photo-electron
tracks for the events triggering the detector.

.. warning::

   The photon list mechanism is a relatively new feature that should be
   considered in beta testing mode. Since it has the potential for becoming
   our final weapon when it comes to detailed simulations of celestial sources,
   it is very important to proceed with all the necessary tests. You are
   welcome to try it and you are kindly encouraged to report as soon as possible
   any issue you might run into.

The photon lists contain all the relevant information about the photons
impinging on the Be window of the GPD, i.e., time, energy, position and
polarization degree and angle. The full format specification is as follows.

.. include:: _phlistfmt.rst
   :start-line: 3

The format is consistent with that expected with the full rewrite of the particle
source class made specifically for this purpose in ``ixpesim``, and requires
``gpdsw`` version 13.10.0 or later.

.. note::

   Since not all the photons impinging on the Be window make it to the active
   gas volume and are absorbed, generating a photoelectron track with enough
   energy to trigger the detector, you should expect the event list in output
   from ``ixpesim`` to be significantly shorter than the photon list that
   is given as an input. The exact ratio depends on the energy spectrum of the
   source. As a rule of thumb, for a typical celestial source you should
   expect about 10% of the photons to generate a track.

   As of ``ixpeobssim`` version 21.3.0 the photon list FITS files include
   by default a copy of the ``SC_DATA`` binary table that is used downstream
   to project the photons back in the sky taking into account the dithering of
   the observatory.


Fitting ixpesim data sets
-------------------------

The ``ixpesim`` output files are not readily usable in XSPEC. Beyond a few
missing header keywords, the main problem is that the measured energy is
expressed in (un-calibrated) PHA counts, and there is no concept of pulse
invariant.

In real life level-1 files are processed by a series of tools, correcting
for a number of effects (gain non-uniformity, secular gain variations,
temperature effects, and charging). For simulated data we can use the same
model that is used at the generation stage of the response functions, as
described in the section about :ref:`irfgen`.

.. seealso:: :ref:`irfgen`

:ref:`reference-xpsimfmt` performs all the steps that are necessary to close the
loop and have the ``ixpesim`` output inter-operate correctly with all the
analysis tools shipped with ``ixpeobssim``---namely:

* add some relevant keywords in the proper headers;
* add the ``PI``, ``Q``, ``U``, ``RA``, ``DEC``, ``X``, ``Y``, and ``W_MOM``
  columns in the ``EVENTS`` extension.

.. warning::

   For the process to round-trip correctly, one should make sure that all the
   relevant adjustable simulation parameters (e.g., the GEM gain) are strictly
   the same as those used in the data sets upon which the response functions are
   based.

   Care should be taken when simulating events with an ``ixpesim`` version
   different from that used for the generation of the response functions.


An exemplary use case
---------------------

Below is some sample code that you might use to generate an ``ixpesim`` simulation
of an arbitrarily complex source model, suitable for a spectro-polarimetric fit
in XSPEC.


.. code::

   # Create the photon list to be fed into ixpesim. Be careful to the version of
   # the IRF you are using, because you will need a non-diagonal response matrix.
   # If you don't know what to do, use the default ``ixpe_mc_v9``
   xpphotonlist.py --config path/to/srcmodel.py \
                   --duration 5000\
                   --irfname "ixpe:obssim:v10"

   # Run the ixpesim simulation. A few notes:
   # - if you want to simulate the three DUs, you will have to run them separately
   # - to process the entire photon list, use a very large number for the -n switch
   # - make sure you have the pressure right for the DU you are targeting.
   #
   # For version 9 of the IRFs, the asymtotic pressure values read:
   # GPD_DME_ASYMPTOTIC_PRESSURE_DICT = {1: 645., 2: 631., 3: 638.}
   ixpesim --src-photon-list ~/ixpeobssimdata/srcmodel_du1_photon_list.fits \
           --output-file ~/ixpeobssimdata/srcmodel_ixpesim_du1.fits \
           -n 1000000
           --dme-pressure 645.

   # Run the event reconstruction.
   ixperecon ~/ixpeobssimdata/srcmodel_ixpesim_du1.fits

   # add the missing columns for the use in XSPEC
   xpsimfmt.py ~/ixpeobssimdata/srcmodel_ixpesim_du1_recon.fits

   # And if you want a reference ``ixpeobssim`` simulation to compare with...
   xpobssim.py --config path/to/srcmodel.py \
               --duration 5000 \
               --deadtime 0. \
               --irfname ixpe_mc_v9

From this point on you should be able to bin the data and fit them in
XSPEC as usual.
