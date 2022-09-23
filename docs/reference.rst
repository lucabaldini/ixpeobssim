.. _reference:

Application reference
=====================

Here is the synopsis of all the applications in the IXPE observation simulation
framework, along with the complete, up-to-date summary of the corresponding
command-line switches.

Simulation facilities
---------------------

.. _reference-xpmdp:
xpmdp
~~~~~
.. command-output:: python ../ixpeobssim/bin/xpmdp.py --help
   :shell:

.. _reference-xpobssim:
xpobssim
~~~~~~~~
.. command-output:: xpobssim --help










.. _reference-xpancrkey:

xpancrkey
---------
.. command-output:: xpancrkey.py --help



.. _reference-xpbin:

xpbin
-----
.. command-output:: xpbin.py --help



.. _reference-xpbinview:

xpbinview
---------
.. command-output:: xpbinview.py --help



.. _reference-xpcalib:

xpcalib
-------
.. command-output:: xpcalib.py --help



.. _reference-xpchrgmap:

xpchrgmap
---------
.. command-output:: xpchrgmap.py --help



.. _reference-xpevtstat:

xpevtstat
---------
.. command-output:: xpevtstat.py --help



.. _reference-xpexposure:

xpexposure
----------
.. command-output:: xpexposure.py --help



.. _reference-xpgrppha:

xpgrppha
--------
.. command-output:: xpgrppha.py --help



.. _reference-xpirfview:

xpirfview
-----------
.. command-output:: xpirfview.py --help







.. _reference-xpobsview:

xpobsview
---------
.. command-output:: xpobsview.py --help


.. _reference-xpophase:

xpophase
--------
.. command-output:: xpophase.py --help



.. _reference-xpphase:

xpphase
-------
.. command-output:: xpphase.py --help


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

   See `this issue <https://bitbucket.org/ixpesw/ixpeobssim/issues/440>`_
   for more details.


.. _reference-xpphotonlist:

xpphotonlist
------------
.. command-output:: xpphotonlist.py --help



.. _reference-xppimms:

xppimms
-------
.. command-output:: xppimms.py --help



.. _reference-xppiscale:

xppicorr
--------
.. command-output:: xppicorr.py --help



.. _reference-xpselect:

xpeselect
---------
.. command-output:: xpselect.py --help



.. _reference-xpsimfmt:

xpsimfmt
--------
.. command-output:: xpsimfmt.py --help



.. _reference-xpsimspec:

xpsimspec
---------
.. command-output:: xpsimspec.py --help



.. _reference-xpsonify:

xpsonify
--------
.. command-output:: xpsonify.py --help



.. _reference-xpsrccoords:

xpsrccoords
-----------
.. command-output:: xpsrccoords.py --help


.. _reference-xpstokesalign:

xpstokesalign
-------------
.. command-output:: xpstokesalign.py --help



.. _reference-xpstokesrandom:

xpstokesrandom
--------------
.. command-output:: xpstokesrandom.py --help



.. _reference-xpstokesshuffle:

xpstokesshuffle
---------------
.. command-output:: xpstokesshuffle.py --help



.. _reference-xpstokessmear:

xpstokessmear
-------------
.. command-output:: xpstokessmear.py --help



.. _reference-xpstripmc:

xpstripmc
---------
.. command-output:: xpstripmc.py --help




.. _reference-xpvisibility:

xpvisibility
------------
.. command-output:: xpvisibility.py --help



.. _reference-xpxspec:

xpxspec
-------
.. command-output:: xpxspec.py --help
