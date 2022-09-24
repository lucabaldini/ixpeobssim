.. _reference:

Application reference
=====================

Here is the synopsis of all the applications in the IXPE observation simulation
framework, along with the complete, up-to-date summary of the corresponding
command-line switches.



Simulation facilities
---------------------

.. _reference-xpobssim:

xpobssim
~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpobssim.py --help
   :ellipsis: 0,11
   :shell:

.. _reference-xpphotonlist:

xpphotonlist
~~~~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpphotonlist.py --help
   :ellipsis: 0,11
   :shell:


.. warning::

   ``xpphotonlist`` is designed to generate list of photons at the top of
   the focal plane, to be fed into the full simulation of the focal-plane
   detectors---which we are not planning to release to the public. This
   application is, therefore, of limited use for the Community; yet, it
   can be used as an inspiration for a useful feature that could be adapted
   to different applications.


.. _reference-xpcalib:

xpcalib
~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpcalib.py --help
   :ellipsis: 0,11
   :shell:



Sensitivity estimation
----------------------

.. _reference-xpmdp:

xpmdp
~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpmdp.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xppimms:

xppimms
~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xppimms.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpvisibility:

xpvisibility
~~~~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpvisibility.py --help
   :ellipsis: 0,11
   :shell:



Analysis facilities
-------------------

.. _reference-xpbin:

xpbin
~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpbin.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpselect:

xpeselect
~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpselect.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpxspec:

xpxspec
~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpxspec.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpphase:

xpphase
~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpphase.py --help
   :ellipsis: 0,11
   :shell:

.. warning::

  When using ``xpphase`` from command line, since the derivatives of the
  frequence are typically (small) negative numbers, it is customary to bump
  into an odd corner of the Python
  `argparse <https://docs.python.org/3/library/argparse.html>`_ module, where
  the "e" character of the exponent specifier, in conjunction with the leading
  minus sign, tricks Python into thinking that the value for the ``nudot0``
  and/or the ``nuddot`` command line arguments are actually a separate option.
  The deal, here, is to use, e.g., the ``nudot0=-1.e13`` form of the options
  specification, `with the equal sign`.

  See `this issue <https://github.com/lucabaldini/ixpeobssim/issues/440>`_
  for more details.


.. _reference-xpophase:

xpophase
~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpophase.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpstokesalign:

xpstokesalign
~~~~~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpstokesalign.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpstokesrandom:

xpstokesrandom
~~~~~~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpstokesrandom.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpstokesshuffle:

xpstokesshuffle
~~~~~~~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpstokesshuffle.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpstokessmear:

xpstokessmear
~~~~~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpstokessmear.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpexposure:

xpexposure
~~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpexposure.py --help
   :ellipsis: 0,11
   :shell:



Visualization facilities
------------------------

.. _reference-xpbinview:

xpbinview
~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpbinview.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpirfview:

xpirfview
~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpirfview.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpobsview:

xpobsview
~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpobsview.py --help
   :ellipsis: 0,11
   :shell:



Miscellanea
-----------

.. _reference-xpsonify:

xpsonify
~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpsonify.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpevtstat:

xpevtstat
~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpevtstat.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpgrppha:

xpgrppha
~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpgrppha.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpsrccoords:

xpsrccoords
~~~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpsrccoords.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpancrkey:

xpancrkey
~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpancrkey.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpchrgmap:

xpchrgmap
~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpchrgmap.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xppicorr:

xppicorr
~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xppicorr.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpsimfmt:

xpsimfmt
~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpsimfmt.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpsimspec:

xpsimspec
~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpsimspec.py --help
   :ellipsis: 0,11
   :shell:


.. _reference-xpstripmc:

xpstripmc
~~~~~~~~~
.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpstripmc.py --help
   :ellipsis: 0,11
   :shell:
