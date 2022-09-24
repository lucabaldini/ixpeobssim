.. ixpeobssim documentation master file, created by
   sphinx-quickstart on Wed Dec  9 11:52:55 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ixpeobssim
==========

``ixpeobssim`` is a simulation and analysis framework specifically developed for the
`Imaging X-ray Polarimetry Explorer <https://ixpe.msfc.nasa.gov/>`_ (IXPE).
Given a source model and the response functions of the telescopes, it is designed
to produce realistic simulated observations, in the form of event lists in FITS
format, containing a strict superset of the information included in the publicly
released IXPE data products.

The core simulation capabilities are complemented by
a suite of post-processing applications which support the spatial, spectral, and
temporal models needed for analysis of typical polarized X-ray sources, allowing
for the implementation of complex, polarization-aware analysis pipelines, and
facilitating the inter-operation with the standard visualization and analysis
tools traditionally in use by the X-ray community.

.. important::

   ``ixpeobssim`` is made available and maintained on a best effort basis under
   the `GNU General Public License v3 <https://www.gnu.org/licenses/gpl-3.0.html>`_,
   in the hope that it will be useful, but it is not officially maintained by
   either `NASA <https://www.nasa.gov/>`_ or `ASI <https://www.asi.it/>`_.

   The `IXPE page on the HEASARC <https://heasarc.gsfc.nasa.gov/docs/ixpe/>`_
   should be considered the ultimate source of information for everything concerning
   the observatory and the analysis of the science data it provides. Should there
   be any discrepancy between the information included in the ``ixpeobssim`` code
   and documentation and that provided by the HEASARC, the latter should be assumed
   to be correct. (And we kindly ask users to report any such instance on our
   `issue tracker <https://github.com/lucabaldini/ixpeobssim/issues>`_, so that we
   can fix it.)


The ``ixpeobssim`` development takes place on a
`git repository <https://github.com/lucabaldini/ixpeobssim>`_ hosted on on
github.
The :ref:`installation`, :ref:`overview` and :ref:`quickstart` pages should be
enough to get you up and running. If that is not the case, or should you encounter
any issue at all, we encourage you to use our
`issue tracker <https://github.com/lucabaldini/ixpeobssim/issues>`_ to file for
bug reports and feature requests.

For convenience, a static version of the latest documentation in pdf format is available
at `this link <https://ixpeobssim.readthedocs.io/_/downloads/en/latest/pdf/>`_.


.. admonition:: Latest changes

   .. include:: release_notes.rst
       :start-line: 5
       :end-line: 20

   ...

   See the :ref:`release_notes` for a full change log.


If you found ``ixpeobssim`` useful, feel free to include an acknowledgment in
papers and/or presentations---`here
<https://www.sciencedirect.com/science/article/pii/S2352711022001169>`_
is the most up-to-date reference.


.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Introduction

   overview
   showcase
   installation
   quickstart
   data_format

.. toctree::
   :hidden:
   :caption: Simulating observations

   tutorial
   reference
   source_models
   toy_models
   dc1
   large_files
   irf
   orbit

.. toctree::
   :hidden:
   :caption: Analysis tools

   filtering
   binning
   xspec
   pipeline
   benchmarks
   broadband_polarization
   interoperability

.. toctree::
   :hidden:
   :caption: Specific applications

   magnetar_models
   extended_sources
   binary_systems
   sonification

.. toctree::
   :hidden:
   :caption: Advanced topics

   irfgen
   ixpesim

.. toctree::
   :hidden:
   :caption: Developer zone

   implementation
   api
   code_development
   release_notes
   about
