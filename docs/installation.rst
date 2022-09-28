.. _installation:

Installation
============

If you are in a rush, and to make a long story short, as of version
29.2.0 ixpeobbsim is hosted on `PyPI <https://pypi.org/project/ixpeobssim/>`_
and you can install it via pip

.. code-block:: shell

   pip install ixpeobssim --user

This should install all the necessary dependences and work out of the box
(we recommend using the ``--user`` option not to pollute everybody else's
environment on the host machine, but you're definitely free to opt for
a system-wide installation by simply omitting this specific option).

.. warning::

   ixpeobssim includes `NumPy <http://www.numpy.org/>`_ as a binary dependence,
   and it is generally not recommended to install Numpy from pip. Should you
   encounter problems with the Numpy installation, you might consider to use a
   `science-ready version of Python <https://scipy.org/install/#scientific-python-distributions>`_
   that comes with NumPy (and the associated Python scientific ecosystem) built-in,
   like the `Anaconda <https://www.anaconda.com/products/distribution>`_ distribution.

That said, most of the remaining of this section is left to the reader for
historical reasons, and in the hope it will be useful to debug installation
problems, should any arise.


Prerequisites
-------------

The package is based on the `Python <https://www.python.org/>`_ scripting
language and the `SciPy <http://www.scipy.org/>`_ Python-based ecosystem.
You will need a working Python installation including several different
third-party packages, most notably:

* `NumPy <http://www.numpy.org/>`_: the fundamental package for scientific
  computing with Python.
* `SciPy <http://www.scipy.org/>`_: a Python-based ecosystem of open-source
  scientific software.
* `matplotlib <http://matplotlib.org/>`_: a Python 2D plotting library.
* `Astropy <http://www.astropy.org/>`_: a Python astronomy package (including
  tools to manipulate FITS files).
* `regions <https://astropy-regions.readthedocs.io/en/latest//>`_: an Astropy
  affiliated package for handling region (including ds9 region files).
* `skyfield <https://rhodesmill.org/skyfield/>`_: a package computing
  positions for stars, planets and satellites.
* `PyXspec <https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/>`_:
  the Python binding for XSPEC (optional).

Loosely speaking you should be able to open the Python terminal and execute the
following ``import`` statements with no errors.

>>> import numpy
>>> import scipy
>>> import matplotlib
>>> import astropy
>>> import regions
>>> import skyfield

If any of the required packages fails to import, take a deep breath and fix the
issue before you move on.


Dependencies in depth
~~~~~~~~~~~~~~~~~~~~~

.. warning::
   ixpeobbsim as a project started in late 2015 when requiring users to
   install Python 3 was more a political decision that anything else. Over the
   years, it has evolved along with the Python ecosystem and the Python language
   itself.

   Although we have not deliberately taken any action to break backward
   compatibility with older versions of the dependencies (and in fact, you might
   very well find that ixpeobbsim is still working with Python 2), you are
   strongly encouraged to make sure that your setup matches the indications in
   this section, as new features will not be tested against older environments.

We recommend installing Python 3.6 or later. In case you missed it, the Python
2x series reached the `end of life <https://www.python.org/doc/sunset-python-2/>`_
on January 1, 2020, which means that no new bug reports, fixes, or changes will
be made to Python 2, and Python 2 is no longer supported. If you are still using
it, you should definitely consider upgrading to Python 3 for your projects.

As far as the dependencies are concerned, our requirement specification reads:

.. literalinclude:: ../requirements.txt
   :end-before: pytest-cov

You can easily find out which version of any specific package you are using
through the interactive Python interpreter, e.g.:

.. code-block:: shell

   [lbaldini@nblbaldini ixpeobssim]$ python
   Python 3.8.3 (default, May 15 2020, 00:00:00)
   [GCC 10.1.1 20200507 (Red Hat 10.1.1-1)] on linux
   Type "help", "copyright", "credits" or "license" for more information.
   >>> import numpy
   >>> numpy.__version__
   '1.18.4'
   >>>

In addition, there is a unit test in the test folder (``test/test_environment.py``)
that you can run to diagnose your environment.


Installing Python through Anaconda
----------------------------------

.. tip::
   Making a long story short: we encourage people to use Python 3 and
   install Python and the associated ecosystem through
   `Anaconda <https://www.anaconda.com/download/>`_. If you know what you
   are doing you can find your way through the framework through alternative
   ways (e.g., the Python installation provided by your GNU/Linux distribution)
   feel free to do so, but we do recommend Anaconda as a decent,
   platform-independent Python experience.

You can download the installer for your os from the `download <https://www.anaconda.com/download/>`_
page (make sure you pick the Python 3 branch) and follow the
`installation instructions <https://docs.anaconda.com/anaconda/install/>`_.

In a nutshell, Anaconda will create on your hard drive a self-contained
directory structure with the Python executable and associated goodies, and all
you have to do to use it is to prepend the folder containing the Python
executable itself to your ``PATH`` environmental variable. The installer is
pretty simple in that it comes in the form of a single bash script (under
GNU/Linux and Mac) or a self-extracting installer (under Windows) that you
execute and does everything for you---just follow the instructions.

It is your responsibility to have the proper environment set up for Anaconda
Python. As mentioned in the previous paragraph, this simply boils down to have
the ``PATH`` environmental variable pointing to the folder containing the
Python executable, which is

* the ``bin`` folder in the Anaconda installation folder under GNU/Linux and
  Mac;
* the top-level Anaconda installation folder under Windows.

Once you have correctly done that, you should be able to launch Python from
the shell or the command prompt and see something in response along these
lines:

.. code-block:: shell

   [lbaldini@nblbaldini doc]$ python
   Python 3.6.4 |Anaconda, Inc.| (default, Jan 16 2018, 18:10:19)
   [GCC 7.2.0] on linux
   Type "help", "copyright", "credits" or "license" for more information.
   >>>


Installing the Python dependencies
----------------------------------

The default Anaconda installer comes with most of the Python scientific
ecosystem that you need (i.e., numpy, scipy, matplotlib and astropy among
many other packages).

For the missing bits, Anaconda comes with the ``pip`` utility that allows you
to install them in a matter of seconds (with the exception of the Python
bindings for Xspec, see the following section). If you're starting from a fresh
installation you will typically need to do

.. code-block:: bash

    pip install regions --user
    pip install skyfield --user

and you should be all set.

.. warning::

   If you have an error when trying and install regions, you might want to
   try and do

   .. code-block:: shell

      pip install wheel --user

   (There are more extensive resources on the web about this very issue,
   including `this <https://stackoverflow.com/questions/34819221>`_)

   Under Windows, depending on your installation, you might come across
   the infamous

   .. code-block:: shell

      error: Unable to find vcvarsall.bat

   while installing regions or skyfield. This is essentially because pip is trying
   to compile part of the library, and you can find all the gory details of the
   issue `here <https://blogs.msdn.microsoft.com/pythonengineering/2016/04/11/unable-to-find-vcvarsall-bat/>`_. In a nutshell, assuming that you are using
   Python 3.5 or later, you will need to install the `Visual C++ Build Tools 2015 <http://go.microsoft.com/fwlink/?LinkId=691126/>`_.


Installing the Xspec Python bindings
------------------------------------

This is an interesting corner of the pre-requisite installation, as PyXspec is
used in some of the ixpeobssim analysis-related modules, and it is
not readily available through Anaconda.

.. note::

   Getting PyXspec up and running might prove less than trivial depending on
   your OS and the details of your installation. It is important to point
   out that you can get away through the vast majority of the ixpeobssim
   functionalities without PyXspec, and if you have Xspec installed you might
   use that netively after all.

   If, on the other end, you are interested in making full use
   of the ixpeobbsim analysis pipelines, detailed instructions about how
   get the Xspec Python bindings up and running follow.

   At any time you can check whether you have PyXspec up and running on
   your ixpeobssim installation by doing something along the lines of

   .. code-block:: python

      >>> from ixpeobssim.utils.environment import PYXSPEC_INSTALLED
      >>> if PYXSPEC_INSTALLED:
      >>>     import ixpeobssim.evt.xspec_ as xspec_

   Please keep in mind to protect the PyXspec-related sections of your
   code with a similar guard in order not to break everything for the
   users who do not have PyXspec installed.

The PyXspec Python bindings for Xspec are fully integrated into the general
`HEASOFT build procedure <http://heasarc.gsfc.nasa.gov/lheasoft/install.html>`_
and they will be built and installed automatically with the rest of
XSPEC/HEASOFT. Note that the source code distribution of XSPEC is required for
using PyXspec. What you want to do is to go to the HEASOFT
`download page <https://heasarc.gsfc.nasa.gov/lheasoft/download.html>`_ and:

* select the source code distribution software type (STEP 1) for your
  target architecture;
* select the desired packages (XANADU/Xspec at the very minimum) at STEP 2;
* click the submit button to get the installation tarball tailored to your
  choices.

At this point you will have an archive, say ``heasoft-<version>src.tar.gz``,
that you should decompress into your target installation folder. There are
detailed installation instructions for all the major platforms,
including

* `Red Hat-based GNU/Linux <https://heasarc.gsfc.nasa.gov/lheasoft/fedora.html>`_
* `Debian-based GNU/Linux <https://heasarc.gsfc.nasa.gov/lheasoft/ubuntu.html>`_
* `Mac OS X <https://heasarc.gsfc.nasa.gov/lheasoft/macos.html>`_
* `Windows with cygwin <https://heasarc.gsfc.nasa.gov/lheasoft/cygwin.html>`_

You should look carefully at the list of prerequisites for your target OS
before you start compiling the code.


Compiling PyXspec under GNU/Linux
---------------------------------

Compiling is pretty easy, the only real thing to pay attention to being the fact
that by default Xspec will try to create the bindings for Python 2 by default,
which is not what you want if you have been following the instructions upstream.
In order to use Python 3 the correct sequence of events is to configure
the build process

.. code-block:: bash

   cd BUILD_DIR/
   ./configure

and, before you start the compilation, tweak the Xspec configuration file
``heasoft-<version>/Xspec/BUILD_DIR/hmakerc`` and explicitly set the
``PYTHON_INC`` and ``PYTHON_LIB`` variables along the lines of

.. code-block:: bash

   PYTHON_INC="-I/path/to/your/anaconda_py3/include/python3.6m"
   PYTHON_LIB="-L/path/to/your/anaconda_py3/lib -lpython3.6m"

At this point you are ready to trigger the actual build process, i.e.,

.. code-block:: bash

   make
   make install

Refer to the PyXspec `installation page <https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/buildinstall.html>`_ for all the gory details and the
troubleshooting information. Once everything is compiled and installed you will
have to source the proper setup file, e.g.

.. code-block:: bash

   export HEADAS=/path/to/your/installed/heasoft-<version>/<platform>
   . $HEADAS/headas-init.sh


See the section about :ref:`xspec` for a more in-depth discussion about
the interactions between ixpeobssim and XSPEC.


Installing ixpeobssim
---------------------

As mentioned at the beginning, the easiest way to install ixpeobssim is via pip.
If you plan on actively contributing to the software development (as opposed
to just using it) you will need to clone the github repository, as explained
in the :ref:`code-development` page.


Installing ixpeobssim as a user
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Installing ixpeobbsim as a user should be as simple as doing

.. code-block:: shell

   pip install ixpeobssim --user

Loosely speaking, if at this point you can open a Python prompt and do

>>> import ixpeobssim

without getting back an error message like this

>>> import ixpeobssim
    Traceback (most recent call last):
      File "<stdin>", line 1, in <module>
    ImportError: No module named ixpeobssim
>>>

you should be all set. You should be able to call all the executables from
their entry points (without the .py extension) as well, e.g.

.. code-block:: bash

   xpobssim --help

In case, for any reason, you do need to uninstall the package, just type

.. code-block:: bash

   pip uninstall ixpeobssim

If you installed PyXspec and you want to use it, make sure also to define the
HEADAS path to HEASOFT package and run the corresponding setup script, e.g.

.. code-block:: bash

   export HEADAS=/path/to/your/installed/heasoft-<version>/<platform>
   . $HEADAS/headas-init.sh


.. warning::

   When you install ``ixpeobbsim`` in user mode you want to make sure you
   don't have any leftover from previous installation, as that is bound to
   lead to cryptic errors that would be hard to debug. Make sure you uninstall
   any old version and you don't leave any .pyc file behind.

   Likewise, if you ever decide to switch from user to developer mode (see next
   section), make sure you clean up any previous installation first.


Up and running in development mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are actively developing code, having it installed (possibly within
system folders) is unlikely to be what you want. Once you have cloned the
github repository all you really need to do in this case is to make sure that
the root folder of the repository is included in the ``PYTHONPATH``
environmental variable.

ixpeobssim comes with its own minimal shell (:repourl:`setup.sh`,
for GNU/Linux or Mac) and batch (:repourl:`setup.bat`, for Windows) setup
scripts  that you're welcome to take advantage of---you should be up and
running by just doing

.. code-block:: bash

   source setup.sh

under GNU/Linux or Mac, and

.. code-block:: shell

   setup.bat

under Windows. (Mind this will also add the path to the folder with the
ixpeobssim executable scripts to the ``PATH`` evironmental variable.)

.. warning::
   If you're operating in development mode you won't have the entry points
   for the ixpeobssim applications installed, e.g., you won't be able to
   call the executables without the .py extensions.

   Most of the documentation is written from the user standpoint, so keep this
   in mind: whenever you read, e.g.,

   .. code-block:: bash

      xpobssim --help

   if you're working in developer mode you will actually need to type

   .. code-block:: bash

      xpobssim.py --help



Changing the default output directory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, ixpeobssim will create an ``ixpeobssimdata`` folder in your
``HOME`` directory and stuff all the output files (e.g., FITS event lists)
in there---unless you specify the path to the output files from the command
line.

If you want to change the default behavior and select your favourite directory
for storing all the output files, you just need to define the
``IXPEOBSSIM_DATA`` environmental variable pointing to the directory itself.

.. code-block:: bash

   export IXPEOBSSIM_DATA=/path/to/the/new/output/folder
