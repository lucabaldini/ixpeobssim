.. _interoperability:

Inter-operability
=================

This section covers the inter-operability between the ``ixpeobssim`` modules
and data products and the standard tools of the X-ray community.

.. warning::

   This section is a stub and needs to be expanded.


Spectro-polarimetry
-------------------

XSPEC
~~~~~

`XSPEC <https://heasarc.gsfc.nasa.gov/xanadu/xspec/>`_ is the primary benchmark
for the ``ixpeobssim`` tools and data product, and the relevant interfaces are
covered in the section about :ref:`xspec`.


Sherpa
~~~~~~

`Sherpa <https://cxc.cfa.harvard.edu/sherpa/>`_ is the
`CIAO <https://cxc.cfa.harvard.edu/ciao/>`_ modeling and fitting application.
Among other things, it supports the XSPEC spectral models and provides a more
modern (and Pythonic) interface to the underlying facilities.

The inter-operation between ``ixpeobssim`` and Sherpa has been demonstrated at
the prototypical level, and the related activity is tracked at
`issue #291 <https://github.com/lucabaldini/ixpeobssim/issues/291>`_.


3ML
~~~

The `3ML <https://threeml.readthedocs.io/en/stable/>`_ Multi-Mission Maximum
Likelihood framework provides a common high-level interface and model definition,
interfacing under the hood to the official software of the various missions.

3ML provides experimental support to IXPE, which is illustrated, e.g.,
`here <http://xwiki.ssdc.asi.it/xwiki/bin/download/Science%20Analysis%20and%20Simulations/SASWG%20meeting%20agendae/WebHome/210715_threeML_DC1.pdf>`_.


Timing
------

HENDRICS
~~~~~~~~

`HENDRICS <https://hendrics.stingray.science/en/latest/>`_, the
High-Energy Data Reduction Interface from the Command Shell, is sophisticated
timing analysis package for X-ray astronomy.

HENDRICS has been tested against ``ixpeobssim`` simulations, see, e.g.,
`here <http://xwiki.ssdc.asi.it/xwiki/bin/download/Science%20Analysis%20and%20Simulations/SASWG%20meeting%20agendae/WebHome/Quicklook%20timing%20analysis%20of%20simulated%20IXPE%20data.pdf>`_.
