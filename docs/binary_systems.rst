.. _binary_systems:

Binary systems
==============

The class :class:`ixpeobssim.srcmodel.roi.xBinarySource` allows to simulate a
point-like periodic source in a binary system. The basic ingredients to describe
the system are the same as in the simple periodic source configuration but here
the input can be in the more generic form of a parameter file including the system
ephemeris. The spin and orbital ephemeris are handled by the
:class:`ixpeobssim.srcmodel.roi.xOrbitalEphemeris` class.


The ``toy_binary`` configuration file distributed with ``ixpeobssim`` describes
a pulsar in a binary system with a power-law photon spectrum with normalization
that varies as a function of the pulse phase. To illustrate the flexibility of
the code, the polarization degree also varies as a function of both energy and
pulsar phase.


.. literalinclude:: ../ixpeobssim/config/toy_binary.py
   :start-after: from __future__ import print_function, division
