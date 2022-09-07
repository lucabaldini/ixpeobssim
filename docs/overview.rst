.. _overview:

Overview
========

The basic simulation and analysis facilities implemented in ixpeobssim are made
available through the following Python script:

* |xpobssim|: produce a photon list for a specific observation
  time, given a :ref:`source model <source-models>` and a set of
  :ref:`instrument response functions <irf>`;
* |xpphase|: calculate the phase of a periodic source based
  on ephemeris and generates a new FITS file appending the phase column.
* |xpophase|: calculate the orbital phase for binary systems;
* |xpselect|: select sub-samples of photons in a given event file, based
  on the event energy, direction, time or phase, producing a new event file;
* |xpbin|: allows to bin the data using several different algorithms,
  producing as output files :ref:`binned events lists <binning>`;
* |xpbinview|: provide the visualization interface to the binned data given
  as input files;
* |xppimms|: calculate the MDP for a a specific observation
  time, given a :ref:`source model <source-models>` and a set of
  :ref:`instrument response functions <irf>`;
* |xpmdp|: calculate the MDP for a a specific observation
  time, assuming a power-law spectrum with adjustable parameters;
* |xpxspec|: perform a spectro-polarimetric fit in XSPEC.

Where applicable, the data formats are consistent with the common display and
analysis tools used by the community, e.g., the binned count spectra can be
fed into XSPEC, along with the corresponding response functions, for doing
standard spectral analysis (note that the response files used are the *same*
for the simulation and the analysis tasks.)

All the ixpeobssim simulation and analysis tools are fully configurable via
command-line and the corresponding signatures are detailed here. In addition,
ixpeobssim provides a pipeline facility that allow to script in Python all the
aforementioned functionalities (e.g., for time-resolved polarimetry this would
mean: create an observation simulation for the system under study, run
|xpselect| to split the photon list in a series of phase bins, run |xpbin| to
create suitable modulation cubes for each of the data subselections and
analyze the corresponding binned output files).
