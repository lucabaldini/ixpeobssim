.. _datafmt:

Data format
===========

Below is an (automatically-generated) dump of the ``ixpeobssim`` output data
format.

.. note::

   While we have tried since the very beginning of the code development
   to keep the data format of the output event list synchronized with that
   of the actual files that we were planning to distribute to the community,
   ``ixpeobssim`` does not in fact produce plain level-2 files.

IXPE Level-2 data format
------------------------

IXPE filtered level-2 files are fairly minimal, in that contain two
extensions:

* ``EVENTS``: contains the event data;
* ``GTI``: contains the good time intervals for the observation;

The ``EVENTS`` extension comes in the form of a binary table with 10 columns:

* ``TRG_ID`` (1J): the trigger identifier. The use of trigger ID, instead of event ID, is to emphasize
  that the DAQ can discard triggers based on the ROI size.
* ``TIME`` (1D): sum of ``SEC`` and ``MICROSEC`` [s].
* ``STATUS`` (16X): 16-bits of processing status/error flags.
* ``STATUS2`` (16X): 16-bits of processing status/error flags.
* ``PI`` (1J): pixel-equalization and gain-corrected event signal.
* ``W_MOM`` (1E): statistical weight of this event (from moments analysis).
* ``X`` (E): calculated position, projected onto the J2000 tangent plane axis parallel
  to celestial equator using the preliminary aspect correction [pixel].
* ``Y`` (E): calculated position, projected onto the J2000 tangent plane axis parallel
  to celestial equator using the preliminary aspect correction [pixel].
* ``Q`` (D): value of Stokes parameter q in J2000 tangent plane axis.
* ``U`` (D): value of Stokes parameter u in J2000 tangent plane axis.

Note the primary header does not include the WCS information to
map ``X`` and ``Y`` in the sky.
Note a few differences remain between 'EVENTS' columns format,
in data files and in ixpeobssim simulated files.

``ixpeobssim`` event lists
--------------------------

Event lists consists of a primary header and a number of distinct binary tables.
The two main extensions (``EVENTS`` and ``GTI``) mimic the expected content of
real celestial observations, although the former contains many more columns,
partly for historical reasons, and partly as a mean to provide useful debug
information specific to simulated data.

.. note::

   The main interface to event lists, :class:`ixpeobssim.evt.event.xEventFile`
   is designed with the level-2 files in mind, and we strive to use as much as
   possible the information in there, in such a way that the analysis tools
   can inter-operate transparently with both real and simulated data.


Additional extensions, specific to simulated data, include:

* ``MONTE_CARLO``: contains the ground truth (e.g., true energy and sky direction)
  for all the events, which is useful for development and debugging purposes
  (it goes without saying, this will not be available for flight data);
* ``ROITABLE``: contains the mapping from the values in the ``SRC_ID`` column in
  the ``MONTE_CARLO`` extension and the name of the corresponding source in
  region of interest, which is useful to map the events into the corresponding
  model components;
* ``TIMELINE``: contains the timeline of the observation in the form of a series
  of homogeneous epochs with a given SAA and Earth occultation flags (the table
  is optional and its generation is controlled by the ``--timelinedata``
  command-line switch in :ref:`reference-xpobssim`);
* ``SC_DATA``: contains the basic spacecraft data (e.g., position and pointing)
  on a regular time grid (the table is optional and its generation is controlled
  by the ``--scdata`` command-line switch in :ref:`reference-xpobssim`, while
  the step of the grid is controlled by the ``--scdatainterval`` switch);
* ``OCTI``: contains the on-orbit calibration time interval for a given DU
  (the table is only generated if the ``--onorbitcalib`` :ref:`reference-xpobssim`
  command-line switch is set);
* ``CHRG_MAP``: contains the charging map at the end of the observation, in a
  form that can be used as an input for subsequent simulations
  (the table is only generated if the ``--charging`` :ref:`reference-xpobssim`
  command-line switch is set);

While the content on the ``TIMELINE`` and ``SC_DATA`` contain partially overlapping
information (the SAA and Earth occultation flags) they are fundamentally different
and serve two different purpose: the former is the most succinct summary of
the different states a given observation is traversing, while the latter is a
mere sampling of a few interesting quantities on a regularly spaced time grid.

We also emphasize that the GTI and OCTI do not necessarily fill the entirety
of the time intervals that can be in principle allocated for science data taking
and calibration, the exact behavior of the observatory being controlled by the
``--gtiminduration``, ``--gtistartpad``, ``--gtistoppad``, ``--onorbitcaldemult``,
``--onorbitcalminduration``, ``--onorbitcalstartpad``, ``--onorbitcalstoppad``,
and ``--onorbitcalrate`` flags of :ref:`reference-xpobssim`.

.. note::

  ixpeobssim provide a small application, :ref:`reference-xpobsview` to take a
  quick-look to the timeline of an observation.


.. include:: _datafmt.rst


Removing Monte Carlo information
--------------------------------

.. warning::
  .. versionadded:: 16.6.0

  Support for stripping the ground truth from the photon lists simulated
  with ``ixpeobssim`` was added in version 16.6.0 to support the first IXPE
  data challenge. This new feature required a series of modification for the
  inter-operability with the analysis tools shipped with ``ixpeobssim``, and
  should be provisionally considered experimental---please report any issue
  you should encounter with Monte-Carlo-less data sets!


While the Monte Carlo information in the simulated photon lists is extremely
useful for developing analysis tools, having the ability to generate files
without any information that would not be ordinary available for real data can
be useful as well, e.g., to test the inter-operability of any given software tool.


xpstripmc
---------

|xpstripmc| is a simple script processing a series of simulating photon lists
and removing the ground truth information---more specifically:

* the ``MONTE_CARLO`` extension;
* the ``ROITABLE`` extension.

Although this would seem like an operation that should be essentially transparent
to the end user, one should realize that the these two extensions do contain valid
information in the spirit of streamlining the simulation and analysis process.
More specifically, the ``MONTE_CARLO`` extension header contains the ``IRFNAME``
keyword, encapsulating the version of the response functions used for the
simulation; in normal condition this is automatically picked up by all the
tools downstream to ensure that the simulation and the analysis are performed
self-consistently---which is usually a sensible thing to do (in real life the
IRF version will be supplied by the user).

All the usual analysis tools have been properly updated to support photon lists
with no Monte Carlo information, provided that the user supplies the extra
information that is needed. (You should get a sensible error message if that is
not the case.)


Using xpbin
-----------

When using |xpbin| with stripped-down data files, and depending on the particular
binning algorithm being used, you should make sure you provide the relevant
information, where appropriate, through the proper command line options.
More precisely:

* for all the ``PHA1`` related binned files (i.e., un-normalized and normalized
  Stokes spectra), you should provide the name of the response functions to be
  used (via the ``--irfname`` command-line switch) so that this is correctly
  picked up by ``XSPEC`` downstream;
* for the ``CMAP``, ``MDPMAP``, ``MDPMAPCUBE``, ``PCUBE``, ``PMAP``, and ``PMAPCUBE``,
  algorithms, you should also provide the name of the response functions to be used
  (via the ``--irfname`` command-line switch) so that the appropriate effective area
  and modulation factor can be correctly loaded.

It goes without saying: with no Monte Carlo information available you will not
be able to run with the ``--mc`` switch.


Using xpselect
--------------

|xpselect| should work exactly as before.
(The previous remark about the ``--mc`` switch still holds.)
