.. _dc1:

Data Challenge 1
================

The photon lists for the first IXPE Data Challence (DC1) were released on
June 1, 2021.
(The `data repository <https://drive.google.com/drive/folders/1cI1lOzS6E3SkYIfxJQm4nSCkfbYCNx65>`_
contains exactly three files per observation---one for each of the three telescopes.)

All the input models have been merged into the main ixpeobssim repository as
of version 16.16.0, and this page is acting as the documentation for them.

All the model files live in the usual ``ixpeobssim/config`` directory.


Pulsar ephemeris
----------------

Based on the consideration that, under most circumstances, the ephemeris for
known periodic sources will be available when IXPE observes the target, we provide
below those used for the simulation for the two pulsars involved in the data
challenge.


PSR J1708-4008
~~~~~~~~~~~~~~

.. code-block::

  * Source name: PSR J1708-4008
  * R. A.      : 257.204167
  * Dec.       : -40.152778
  * Epoch MJD  : 59823.000000 (MET 178761600.000 s)
  * nu0        : 0.090798237 Hz
  * nudot0     : -1.598360e-13 Hz s^{-1}
  * nuddot     : 0.000000e+00 Hz s^{-2}


RX B1509-58
~~~~~~~~~~~

.. code-block::

  * Source name: RX B1509-58
  * R. A.      : 228.481333
  * Dec.       : -59.135778
  * Epoch MJD  : 59611.000000 (MET 160444800.000 s)
  * nu0        : 6.572649741 Hz
  * nudot0     : -6.582043e-11 Hz s^{-1}
  * nuddot     : 1.918559e-21 Hz s^{-2}


Multi-wavelength context
------------------------

For Mrk 421 we do have quasi-simultaneous optical polarization information available
for the six IXPE observations, as listed in the following table.

=========  ==============  ===============  ================
MJD        R [mag]         Pol. Degree [%]  Pol. Angle [deg]
=========  ==============  ===============  ================
59692.904  12.53 +/- 0.05  2.60 +/- 0.20    39.3 +/- 3.0
59706.556  12.51 +/- 0.05  2.88 +/- 0.20    47.6 +/- 3.0
59720.754  11.97 +/- 0.05  6.01 +/- 0.20    -47.1 +/- 3.0
59728.883  12.44 +/- 0.05  3.14 +/- 0.20    41.0 +/- 3.0
59734.496  12.44 +/- 0.05  2.82 +/- 0.20    42.2 +/- 3.0
59748.431  12.49 +/- 0.05  3.40 +/- 0.20    41.5 +/- 3.0
=========  ==============  ===============  ================


.. _figure-mkr421_optical_polarization:
.. figure:: figures/dc1/dc1_mrk421_optical_polarization.*
   :width: 70%

   Quasi-simultaneous optical polarization for Markarian 421.



Observing Plan
--------------

Cygnus X-1
~~~~~~~~~~

* Source: Cygnus X-1 in the hard state
* Observation pattern: single ~300 ks observation
* Observation start: 2022-04-01T00:00:00.000000

.. _figure-cygx1_visibility:
.. figure:: figures/dc1/cygx1_visibility.*
   :width: 70%

   Visibility for Cygnus X-1 in 2022.


.. automodule:: ixpeobssim.config.dc1_cygx1



Markarian 421
~~~~~~~~~~~~~

* Source: Markarian 421
* Observation pattern: six ~50 ks observations spaced by ~two weeks
* Observation start:

  * 2022-04-23T00:00:00.000000 (obs_id 1)
  * 2022-05-07T00:00:00.000000 (obs_id 2)
  * 2022-05-21T00:00:00.000000 (obs_id 3)
  * 2022-05-28T00:00:00.000000 (obs_id 4)
  * 2022-06-04T00:00:00.000000 (obs_id 5)
  * 2022-06-19T00:00:00.000000 (obs_id 6)

.. _figure-mkr421_visibility:
.. figure:: figures/dc1/mrk421_visibility.*
   :width: 70%

   Visibility for Markarian 421 in 2022.


.. automodule:: ixpeobssim.config.dc1_mrk421



MSH 15-52
~~~~~~~~~

* Source: MSH 15-52
* Observation pattern: single 1.5 Ms observation
* Observation start: 2022-02-01T00:00:00.000000


.. _figure-msh1552_visibility:
.. figure:: figures/dc1/msh1552_visibility.*
   :width: 70%

   Visibility for MSH 15-52 in 2022.



.. automodule:: ixpeobssim.config.dc1_msh1552



RX J1708-4008
~~~~~~~~~~~~~

* Source: RX J1708-4008
* Observation pattern: single 1 Ms observation
* Observation start: 2022-09-01T00:00:00.000000


.. _figure-1rxs_j1708_visibility:
.. figure:: figures/dc1/1rxs_j1708_visibility.*
   :width: 70%

   Visibility for PSR J1708-4008 in 2022.


.. automodule:: ixpeobssim.config.dc1_1rxs_j1708
