.. _magnetar_models:

Magnetar models
===============

Roberto and Roberto kindly provided the magnetar models described in
`Taverna et al. 2020 <https://ui.adsabs.harvard.edu/link_gateway/2020MNRAS.492.5057T/EPRINT_PDF>`_
in a tabular form suitable for fitting with XSPEC. The data files are available on
our `auxiliary file repository <https://bitbucket.org/ixpesw/ixpeobssim_auxfiles>`_
and, as of version 14.0.0, ixpeobssim provides the necessary infrastructure
to read and use them.

.. seealso::

  :ref:`large_files`.


Model description
-----------------

The models come as three separate FITS files for the I, Q, and U Stokes parameters,
the basic structure of each one being of the form

.. code-block::

  No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU      15   ()
  1  PARAMETERS    1 BinTableHDU     47   6R x 10C   [12A, J, E, E, E, E, E, E, J, PE(13)]
  2  ENERGIES      1 BinTableHDU     21   49R x 2C   [D, D]
  3  SPECTRA       1 BinTableHDU     21   314496R x 2C   [6E, 49E]

(The actual dimensions of the cards might change, but the very fundamental
structure of the files is fixed by the OGIP standard.)

The ``PARAMETERS`` binary table encapsulate the following model parameters:

* `Model`: the underlying Physical model, i.e., ``ATMOSPHERE``, ``BLACKBODY``,
  ``SOLID_SURFACE_FREE_IONS`` and ``SOLID_SURFACE_FIXED_IONS``.
* `chi`: angle between the line of sight and the rotation angle of the star.
* `xi`: angle between the magnetic axis ans the rotation angle of the star.
* `delta_phi`: twist angle.
* `beta`: electron velocity in units of c.
* `phase`: the phase bins.

In the initial implementation we have 4 physical models x 13 chi x 8 xi x 12
delta_phi x 7 beta x 9 phase bins, that is, 314496 combinations.
Spelling things out explicitly:

.. code-block::

  Model     -> [1, 2, 3, 4]
  chi       -> [0.1, 15, 30, 45, 60, 75, 89.9, 105, 120, 135, 150, 165, 179.9]
  xi        -> [0.1, 15, 30, 45, 60, 75, 85, 89.9]
  delta_phi -> [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]
  beta      -> [0.2, 0.3, 0.34, 0.4, 0.5, 0.6, 0.7]
  phase     -> [0.350066, 1.0482, 1.74633, 2.44446, 3.14259,
                3.84072, 4.53886, 5.23699, 5.93512]

The ``ENERGY`` binary extensions encapsulates the energy binning, and in the
initial implementations has 49 bins between about 0.5 keV and 10 keV.

.. warning::

  Some of the models can be noisy at high energy, and the fact that tables are
  limited to less than 10 keV might result in unstable behaviors if you actually try and
  extrapolate to higher energies. While this should not be an issue for a
  standard in the nominal 2--8 keV energy band, be aware of this potential
  issue.


The ``SPECTRA`` extension contains the spectra for all the possible configurations.
It comes with two columns, the first of which contains the 6 model parameters and
the second the 49 spectral values.
Note the ``SPECTRA`` extension is filled by looping over the model parameters in
reverse order, i.e., from `phase` to `Model`.


Parametrization details
^^^^^^^^^^^^^^^^^^^^^^^

.. note::

  This section contains some technical details about how we parse, store, and
  manipulate the data in the input FITS tabular model. You can probably skip
  at a first reading, but if you do run into issues, this might turn out to be
  useful for debugging.

All the underlying model quantities are binned in 49 bins in energy between 0.5
and 10 keV and 9 bins in phase between 0 and 2 pi, with the corresponding values
in the energy grids where the Stokes spectra are sampled representing the centers
of the bins. As it turns out, the first and the last energy bins act as
underflow/overflow accumulators in the underlying Monte Carlo simulation, and
therefore the corresponding spectral and polarimetric quantities are un-phyisical
and need to be discarded. This leaves us with I, Q and U tabulated on the
following grids.

.. code-block::

  # Energy grid in keV.
  [0.549539  0.5841565 0.6209545 0.6600705 0.701651  0.74585   0.7928335
   0.842777  0.8958665 0.9523005 1.0122895 1.076055  1.14384   1.215895
   1.29249   1.37391   1.460455  1.552455  1.65025   1.754205  1.864705
   1.98217   2.10704   2.23977   2.38086   2.530835  2.69026   2.85973
   3.03987   3.231365  3.434925  3.6513    3.88131   4.125805  4.3857
   4.661975  4.95565   5.267825  5.59966   5.9524    6.32737   6.72595
   7.14964   7.60002   8.07877   8.587685  9.12865  ]

  # Phase grid divided by 2pi.
  [0.05571473  0.1668262   0.27793705  0.38904786 0.50015873 0.61126953
   0.72238195  0.83349282  0.94460368]

The first slight complication is that if we create an interpolated bi-variate
spline on this grid, since technically we are `extrapolating` at phases 0 and 1,
there is no guarantee that the two values will be identical, or even close---i.e.,
in general we will have discontinuities when the phase roll over.
In order to fix this we add two values to both sides of the energy grid,
namely

.. math::
   \max(\phi) - 1 \quad \text{and} \quad 0

on the left side, and

.. math::
   1 \quad \text{and} \quad 1 + \min(\phi).

In addition, we impose boundary conditions

.. math::
   I(\max(\phi) - 1) = I(\max(\phi)) \quad
   I(1 + \min(\phi)) = I(\min(\phi)) \quad
   I(0) = I(1) = \frac{1}{2}\left( I(\min(\phi)) + I(\max(\phi)) \right)

(and the equivalent for Q and U) such that, by using a spline of order 2 for the
interpolation, we are guaranteed that there are no discontinuities when the
phase rolls over.

The other issue is that the energy grid is effectively limited to just above
9 keV on the right end, while ixpeobssim, by default, simulates photons between
1 and 12 keV in energy. While this is not a major concern, as photons between
9 and 12 keV are relatively sparse and, even when the energy dispersion is taken
into account, they hardly have any effect on any analysis in the 2--8 keV
energy range, we need a sensible way to extrapolate the spectra. We accomplish
this internally at the stage where the model table is interpolated in the parameter
space to create 2-dimensional splines in the phase-energy space to be used for
the simulation. More precisely, this is a two-step process where:

* I is fitted with a power law between 6 and 9 keV (adjustable) and the fit model
  is used to extrapolate at high energy;
* the reduced Stokes parameters Q/I and U/I are fitted with a straight line
  in the same interval, and the linear model is used for the extrapolation.




Simulation interface
--------------------

ixpeobssim provides the :class:`ixpeobssim.srcmodel.magnetar.xMagnetarModelsT2020`
class as the basic interface the magnetar table models described above.
Instantiating the model table takes one Python line.

.. code-block:: python

  from ixpeobssim.srcmodel.magnetar import xMagnetarTableModelT2020

  # Load the table model. If you need to download any of the large auxiliary
  # FITS files, you should be prompted with detailed instructions (and the program
  # will exit).
  table = xMagnetarTableModelT2020()

This will load all the input FITS files, make some integrity check on the
underlying binary tables, and cache the necessary arrays for later use.
At this point you have the full I, Q and U spectra for the entire range of the
model parameters and the of the magnetar rotational phase.

At simulation time, you typically want to calculate the actual spectro-polarimetric
model for a given set of input parameters, in a form that can be fed
directly into the simulation code. ixpeobssim provides two main interfaces,
allowing to find the *nearest* model to a set of input parameters, or to
*interpolate* the underlying model table to the aforementioned parameters.
Typically you want to do the second, and the following snippet should be fairly self-explanatory.


.. code-block:: python

  from ixpeobssim.srcmodel.magnetar import xMagnetarModelsT2020, xMagnetarTableModelT2020


  # Source name and position.
  source_name = 'AXP 1RXS J1708'
  ra, dec = 257.2042, -40.1528

  # Source period.
  nu0 = 0.09089812328

  # Magnetar model parameters.
  model = xMagnetarModelsT2020.BLACKBODY
  chi = 89.9
  xi = 60.
  delta_phi = 0.5
  beta = 0.34

  # Phase-averaged integral flux between 2 and 10 keV, in erg/cm2/s.
  flux = 3.5e-11

  # Load the magnetar model table.
  table = xMagnetarTableModelT2020()

  # Interpolate the model table to the target parameters.
  spec, pol_deg, pol_ang = table.interpolate(model, chi, xi, delta_phi, beta, flux)

  # Define the actual ROI model.
  ROI_MODEL = xROIModel(ra, dec)
  ephem = xEphemeris(0., nu0)
  src = xPeriodicPointSource(source_name, ra, dec, spec, pol_deg, pol_ang, ephem)
  ROI_MODEL.add_source(src)

.. note::
   For completeness, a full model table with the QED effects switched off
   is available for debugging and illustration purposes.
   The usage is exactly the same, with the class named ``xMagnetarTableModelT2020DedOff``.
