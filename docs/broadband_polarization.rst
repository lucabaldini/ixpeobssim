.. _broadband-polarization:

Broadband polarization
======================

Estimating the broadband polarization degree and angle, or Stokes parameters,
from a (real or simulated) observation is a fairly common task---possibly the
most common one. This section is devoted to a succinct discussion of the various
means ixpeobssim is offering for the job and the associated merits and pitfalls.

.. warning::
   The section is largely focused on point sources, at the moment. Updates for
   extended sources are on their way.


(For completeness, by the term *broadband*, in this context, we mean *averaged
over an energy range much larger than the intrinsic binning of the underlying
response functions*.)


.. note::
   This is a slippery business, and mileage might vary considerably depending on
   your precise analysis setup, so take the whole discussion in this section
   *cum grano salis* and make an effort to gauge its validity to the problem at
   hand.

   To make a long story short, going through the effort of a full
   spectro-polarimetric fit with a proper model is likely to give you the right
   answer, but there are simpler and faster ways around that might bring
   you close enough in specific circumstances.


The ixpeobssim toolbox
----------------------

ixpeobssim offers several different tools designed to estimate the broadband
polarization, given a simulated photon list. We briefly recap the laundry list,
here; more details can be found in the chapters about :ref:`binning` and
:ref:`xspec`.


Fitting with XSPEC
~~~~~~~~~~~~~~~~~~

Whenever practically possible, one should fit the I, Q and U Stokes spectra with
a proper spectro-polarimetric model (e.g., in XSPEC). This is properly
taking into account both the effective area and the energy dispersion, and is
supposed to provide the ultimate statistical precision---at least for point
sources.

Fitting with a *constant* polarimetric multiplicative model (i.e., ``polconst``
in the ixpeobssim language) will often yield sensible results in terms of the
average broadband polarization degree, but it is not guaranteed to be unbiased.

If a suitable spectral model, for any reason, is not available, fitting the
normalized Stokes parameters Q/I and U/I with an additive polarimetric model is
a possible alternative, although it is not guaranteed to be unbiased, either---the
reason being that, strictly speaking, the detector response matrix is not
directly applicable to the normalized Stokes spectra.

.. seealso::
   :ref:`xspec`.


Polarization cubes
~~~~~~~~~~~~~~~~~~

Polarization cubes are a legitimate, model-independent, alternative to a fully
fledged forward-folding fit. The automatically incorporate the detector
acceptance correction, but, by their very nature, cannot handle the energy
dispersion. Whether this is an important effect or not depends very much on the
setup under study, but the good news is that one can gauge the magnitude of any
possible bias by comparing polarization cubes binned in measured and true energy.
*Polarization cubes binned in true energy are supposed to provide the right answer*.

(It goes without saying that the true energy is not available in real life, but
still this is a perfectly valid thing to do for sensitivity studies and/or
debugging purposes.)

.. seealso::
   :ref:`binning`.


Modulation cubes
~~~~~~~~~~~~~~~~

In a nutshell: *do not use them*. Modulation cubes perform the relevant sums
and/or averages over the count spectrum (as opposed to a proper proxy of the
input source spectrum), and do not handle neither the effective area nor the
energy dispersion of the detector. Binning in true energy will save you from the
latter, but not from the former, which is typically a bigger effect.

Modulation cubes are retained in ixpeobssim, mainly for diagnostic purposes,
but their use for science analysis is deprecated.


A toy case study
----------------

In the spirit of substantiate the somewhat generic statements in the previous
section, the following plot summarizes the performance of the various methods
for retrieving the broadband polarization (between 2 and 8 keV) for a toy
setup in which the polarization degree increases linearly with energy and the
polarization angle is constant. (For completeness, the spectrum was a simple
power law.)

We generated 1000 independent, high-statistics (but unrealistic) realizations of
the very same model and, for each single one we applied all the methods described
in the previous section.

For completeness, the model source file is at :repourl:`ixpeobssim/config/toy_pollin.py`
and the full benchmark code at :repourl:`ixpeobssim/benchmarks/toy_pollin.py`.

.. _figure-toy_pollin_polarization_degree_precision:
.. figure:: figures/benchmarks/toy_pollin_polarization_degree_precision.*
   :width: 80%

   Average (over 1000 independent realizations of the same model) of the
   broadband 2--8 keV polarization, estimated by various means. The orange
   dashed line represents the expectation from the input model.

A few comments are in order:

* The full XSPEC fit (with a proper model) to the I, Q and U spectra provides the
  right answer; fitting the normalized Q/I and U/I Stokes parameters to the
  polarimetric part of the model features a measurable negative bias, as does
  a fit with a constant polarization model;
* a polarization cube with a single energy bin between 2 and 8 keV features a
  small, but measurable, negative bias---which is completely recovered
  when binning in true energy;
* the corresponds single-bin modulation cube is way off in both flavors.
