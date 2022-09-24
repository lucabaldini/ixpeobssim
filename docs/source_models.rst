.. _source-models:

Creating source models
======================

.. warning::
   While this is possibly one of the most important pieces of documentation,
   it is still work in progress. Do not hesitate to complain if something
   is not clear incomplete, or simply wrong.

Source models are specified in ixpeobssim through Python configuration
files---effectively simple Python modules. As we shall see in a second, this is
the key of the flexibility in terms of source specification provided by the
simulation framework. All the source model configuration files shipped with
ixpeobssim are included in the :repourl:`ixpeobssim/config` folder, and that is
a good starting point to get up and running with creating your own
observation-simulation setup.

In a nutshell, the one thing that a source configuration file `must` contain
is an instance of the :class:`ixpeobssim.srcmodel.roi.xROIModel` class
called ``ROI_MODEL`` (mind that the name and case of the variable are important,
as when you run a simulation the relevant application will try and import that
exact name from the Python configuration module---and fail miserably if that is
not possible). In the rest of this section we shall walk through the relevant
classes for the definition of a region of interest (ROI) to be simulated.



Spectra and polarization patterns
---------------------------------

Understanding how photon spectra and polarization patterns are specified in
ixpeobssim is key to being able to create a functional simulation setup.
The two basic rules, here, are:

#. the photon spectrum can be an arbitrary function of energy and time,
   expressed in units of :math:`\text{cm}^{-2}~\text{s}^{-1}~\text{keV}^{-1}`;
#. the polarization degree and angle can be arbitrary functions of energy, time,
   `and` position in the sky.

From the implementation standpoint, the photon spectrum and the polarization
degree and angle are simply Python functions with the proper signature:

.. code-block:: python

   def spec(E, t):
       """Definition of the photon spectrum (cm^-2 s^-1 keV^-1).
       """
       # Function definition here.

   def pol_deg(E, t, ra, dec):
       """Definition of the polarization degree (between 0 and 1).
       """
       # Function definition here.

   def pol_ang(E, t, ra, dec):
       """Definition of the polarization angle (in radians).
       """
       # Function definition here.

These functions are then used internally by the simulation code to generate a
proper photon list, and we shall see in the next section how they are passed
to the relevant classes. We emphasize that this approach is fairly flexible, as
the user can put in the function bodies literally anything: a constant, an
analitic formula, an interpolator constructed from a numeric table, etc.
`The only important thing is to make sure that at runtime the function can
be called with the proper signature` (and, by the way, can operate on numpy
arrays---so always make sure that you use the numpy version of the mathematical
functions that you need in your function body).

.. note::
   At this point you might wonder why the rules are comparatively restrictive
   for the photon spectra which, in contrast to the polarization patterns,
   cannot be dependent on the position on the sky. Well, the reason is that
   a fully arbitrary source model doesn't really fit with the current
   implementation of the basic simulation strategy. But before you panic, a
   couple of things:

   #. for extended sources, as we shall see in the next section, the intensity
      of the emission is determined by the spatial template and therefore
      it depends on the position---in other words it is really the spectral
      shape that is bound to be constant across the source extension;
   #. there are ways to break an extended source in an arbitrary number
      of smaller components, each with its own independent spectral
      shape---which effectively allows to overcome the limitation, albeit in
      a possibly inefficient way;
   #. the Chandra-to-IXPE conversion tool shipped with ixpeobssim provides
      an alternative observation-simluation strategy fully preserving the
      correlations between the sky position and the spectrak shape.

Say that you want, e.g., simulate a stationary source with a power-law energy
spectrum; all you have to do is something along the lines of the following
snippet.

.. code-block:: python

   def spec(E, t):
       """Definition of the photon spectrum (cm^-2 s^-1 keV^-1).

       Change the power-law normalization and index to the values that suit
       you. Also, note that t is not used, i.e., the source is stationary, but
       you are free to put in the function body whatever complicated expression
       of energy `and` time, e.g., a prefactor and an index that vary with time.
       """
       return 10. * (E**-2.)

On a related note, a polarization degree constant in energy, time and
position in the sky is readily specified.

.. code-block:: python

   def pol_deg(E, t, ra, dec):
       """Definition of the polarization degree (between 0 and 1).
       """
       return 0.5

Note that in all cases the signature of the function must be the correct
one, no matter whether any specific dynamical variable is effectively used in
the function body or not.

.. tip::
   Since we don't necessarily want you to reinvent the wheel every time,
   ixpeobssim provides convenience wrappers to power laws and constants, among
   other things. The first and last snippets above would be more succintly
   coded as

   .. code-block:: python

      from ixpeobssim.srcmodel.spectrum import power_law
      from ixpeobssim.srcmodel.polarization import constant

      spec = power_law(10., 2.)
      pol_deg = constant(0.5)

   And you should take a look at the documentation of the two functions for
   more details (note that the power-law wrapper support time-varying
   prefactor and index).



Using XSPEC models
~~~~~~~~~~~~~~~~~~

Arbitrary XSPEC models can be fed into xpobssim, providing that PyXSPEC is
up and running on your system. (If you are in a hurry, the section about
:ref:`toy_models` provides a working example that you can take inspiration from.)

In order to effectively use an underlying XSPEC model, you have to create a
:class:`ixpeobssim.srcmodel.spectrum.xXspecModel` instance. This is essentially
an univariate, interpolated spline (coming with all the standard spline facilities),
with the only exception that it can be called with a phony, second argument
(in addition to the energy) so that it can be readily used to define a source
component---in this case the second argument is acting as the time and has no
effect.

The easiest way to create an XSPEC model is through a generic expression and a
list specifying all the model parameters.

.. code-block:: python

   from ixpeobssim.srcmodel.spectrum import xXspecModel

   spec = xXspecModel('powerlaw', [2., 0.1])

The :class:`ixpeobssim.srcmodel.spectrum.xXspecModel` constructor mimics the
underlying XSPEC bindings and allows passing the parameters in the form of a
dictionary indexed by sequential identifiers, which is handy when one wants to
only set a subset of them and use the XSPEC defaults for the others:

.. code-block:: python

   from ixpeobssim.srcmodel.spectrum import xXspecModel

   spec = xXspecModel('powerlaw', {1: 2., 2: 0.1})
   spec = xXspecModel('powerlaw', {2: 0.1})

Alternatively, the model can be read from a text file with a suitable syntax,
see :meth:`ixpeobssim.srcmodel.spectrum.xXspecModel.from_file()`

.. code-block:: python

   from ixpeobssim.srcmodel.spectrum import xXspecModel

   spec = xXspecModel.from_file('path/to/your/file')


In addition, XSPEC models provide a convenient facility to dump the data points
to an ASCII file, see :meth:`ixpeobssim.srcmodel.spectrum.xXspecModel.save_ascii()`.


.. _model-components:

Model components
----------------

The concept of a `model component` is central in the definition of an
observation simulation setup, and is used throughout ixpeobssim to encapsulate
the relevant properties of a few slightly different, but related, things:

* a simple celestial source that cannot be further decomposed, such as that
  defined in :repourl:`ixpeobssim/config/toy_point_source.py`;
* different physical components of the same source---e.g., the thermal and non
  thermal emission in Cas A as modeled in :repourl:`ixpeobssim/config/toy_casa.py`;
* different sources in the same ROI---e.g., the Crab pulsar and nebula in
  :repourl:`ixpeobssim/config/crab_complex.py`.

From a technical standpoint, a model component is basically an instance of
the :class:`ixpeobssim.srcmodel.roi.xCelestialModelComponentBase` class (or, to
be more precise, an instance of a subclass of the base class) and encapsulates,
among other things, the source name, photon spectrum, and polarization degree
and angle, the three latter quantities being specified at construction time
as explained in the previous section. (Take a look at the signature of the
constructors in the module documentation for more details.)

The subclasses of :class:`ixpeobssim.srcmodel.roi.xModelComponentBase`
deal for the most part (with the exception of that describing periodic sources)
with the source morphology, re-implementing the method for extracting vectors
of positions in the sky (ra, dec). For completeness, the actual subclasses that
can be instantiated and used in simulations are

======================== ====================================================
Source type              Specialized subclass
======================== ====================================================
Point source             :class:`ixpeobssim.srcmodel.roi.xPointSource`
Periodic source          :class:`ixpeobssim.srcmodel.roi.xPeriodicPointSource`
Binary system            :class:`ixpeobssim.srcmodel.roi.xBinarySource`
Uniform disk             :class:`ixpeobssim.srcmodel.roi.xUniformDisk`
Gaussian disk            :class:`ixpeobssim.srcmodel.roi.xGaussianDisk`
Uniform annulus          :class:`ixpeobssim.srcmodel.roi.xUniformAnnulus`
Extended source          :class:`ixpeobssim.srcmodel.roi.xExtendedSource`
Chandra observation      :class:`ixpeobssim.srcmodel.roi.xChandraObservation`
Instrumental background  :class:`ixpeobssim.srcmodel.bkg.xInstrumentalBkg`
Instrumental background  :class:`ixpeobssim.srcmodel.bkg.xPowerLawInstrumentalBkg`
Extragalactic background :class:`ixpeobssim.srcmodel.bkg.xExtragalacticBkg`
Galactic background      :class:`ixpeobssim.srcmodel.bkg.xGalacticBkg`
======================== ====================================================

and provide a fairly broad and diverse selection of morphological
characteristics that should be sufficient for most of the applications.

Now, the simplest possible model component (which in this case is a
fully-fledged source) can be defined with only a few lines of code.

.. literalinclude:: snippets/point_source_definition.py

And, for completeness, the output of the last ``print`` function should
be something along the lines of the following snippet. We're almost there.

.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python snippets/point_source_definition.py
   :shell:
   :ellipsis: 0,1


.. note::

  `ixpeobssim` has a concept of Calibration sources and associated region of
  interest. The reader is referred to the `toy_flat_field.py` and `fcw_calC.py`
  configuration files for usage examples.

  (Note that calibration configuration files are supposed to be run through
  :ref:`reference-xpcalib` rather than :ref:`reference-xpobssim`.)



Regions of interest
-------------------

A complete model of a region of interest is basically a collection of an
arbitrary number of model components---and it is the basic simulation unit in
the context of ixpeobssim.

From an implementation standpoint, it is an instance of the
:class:`ixpeobssim.srcmodel.roi.xROIModel` class, and is specified passing the
coordinates of its center; components are added after the fact. A small tweak
to the last snippet in the previous section

.. literalinclude:: snippets/roi_definition.py

make it for a valid ixpeobsim configuration file that can be used to run an
actual simulation. The ``print`` output will look like this.

.. program-output:: export PYTHONPATH=../:$PYTHONPATH;python snippets/roi_definition.py
   :shell:
   :ellipsis: 0,1


.. _chandra-observation:

Chandra region of interest
--------------------------
The region of interest for a Chandra observation follows the same concept
as described above, with the difference that the spectral and
morphological information of the source is directly taken from the provided
Chandra event list and does not need to be specified. Then the xpobssim tool will
take care to convert the Chandra observation into a correponding IXPE event
list, opportunely down-sampling and smearing the events with the instrument
response functions. This simulation technique has the advantage to preserve the
full correlation between the source morphology and the spectrum, especially
important for extended sources.

From tecnical standpoint, it is an instance of the
:class:`ixpeobssim.srcmodel.roi.xChandraROIModel` class, derived from the
:class:`ixpeobssim.srcmodel.roi.xROIModel` base class, and is initialized by
passing the Chandra event file and choosing the ACIS detector used for the
observation (I or S). Sub-regions of the ROI can be defined (with individual
polarization models) as instances of the
:class:`ixpeobssim.srcmodel.roi.xChandraObservation` class by specifying the id,
polarization degree and angle and the region definition (via a  ds9 region file).

.. literalinclude:: snippets/chandra_roi_definition.py

It is also possible to remove sources in the region of interest by
passing the exclude option in the instance of
the :class:`ixpeobssim.srcmodel.roi.xChandraObservation` class and to add
sources using the standard model components described above. These features are
particularly useful for treating the pileup phenomenon occurring in
Chandra CCD detectors for bright sources, whenever two or more photons are
detected as a single event. In these cases, the simplest strategy is to
remove the bright source suffering of pileup from the conversion ROI and replace
it with a standard model component with an appropriate source model description.
The converter will then take care of it, simulating the source in the standard
way and merging the photon lists at the end of the process. For reference, here
are the additional lines that need to be included to the source model example
above to remove a region and add a point source to the chandra region of
interest:

.. literalinclude:: snippets/chandra_roi_pointsource_definition.py


Source model etiquette
----------------------

Although the ``ROI_MODEL`` top-level object is all ixpeobssim is looking for
in a configuration file (and everything will run happily as long as
you provide that in a valid form) there is undoubtedly value in trying and
keep a consistent style when creating source models. People will get
acquainted to the conventions and will understand new models more easily.
This section is essentially a list of suggestions on how to go about
creating your own models.

We note, in passing, that the section about :ref:`toy_models` is a good read,
along with the material in this section.

Give your model file a sensible name (which, in most cases, should be obvious).
It is conceivable that we might have multiple, alternative, configuration
files for a given source; in this case, use your best judgement to try and
convey with the file name the `intent` of any of these alternative models.
(Having crab.py, crab2.py, crab_new.py, crab_newer.py stops being funny
after a while.) Even more important, try and factor out as much as possible
of the code that is shared across multiple models---resist the temptation to
cut and paste code around!

Try and document your model in a reasonable fashion. This is typically
achieved by putting arbitrary `reStructuredText <http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_ that can be then parsed in
`sphinx <http://www.sphinx-doc.org/en/master/index.html>`_ at the beginning
of the Python configuration module, right after the license notice and before
the ``__future__`` Python imports.

.. warning::
   For the documentation to be parsed by sphinx, the reStructuredText
   `has` to be before any actual Python code. Keep this in mind.

When you run a simulation and you want to compare the ixpeobssim output
with the input model, you will need to have all the input parameters at
hand. Keep this in mind when writing configuration files:

.. code-block:: python

    spec = power_law(1., 2.)

is more succinct than

.. code-block:: python

    pl_norm = 1.
    pl_index = 2.
    spec = power_law(pl_norm, pl_index)

but the second formulation makes it easier for you to retrieve the input
spectral index when you fit your binned ``PHA1`` spectrum in xspec and want to
compare the best-fit parameters with the models.

Take a second to read our :ref:`coding-guidelines` and try and stick to them.
In the long run people will be grateful to you for doing it.

Create a ``display()`` method in each configuration file that, when called,
will plot anything that might be interesting in you model. (Each source model is
different and it is hard to achieve this through a common, completely
automated interface.) This method will not partecipate in the observation
simulation when the model is parsed by ixpeobssim, but is something that,
when protected by a ``if __name__ == '__main__'`` can be useful as a
quick-look interface and for creating figures for the documentation.
To this end, we have implemented a small bootstrapping function that you are
welcome to use by doing something along the lines of

.. code-block:: python

    if __name__ == '__main__':
	from ixpeobssim.config import bootstrap_display
        bootstrap_display()

This will create a custom option parse for you with a limited set of
options (e.g., to save the plots in a specific folder) and will fire up
your ``display()`` method right away. Try at any time

.. code-block:: bash

    python path/to/config/file.py --help

to have an up-to date help output.

Being consistent in terms of axis labels and units, conventions for file names
and all that is cool. Really. Take a seconds to acquaint yourself with the
following small utility functions and classes:

* :meth:`ixpeobssim.utils.matplotlib_.save_gcf`: save the current figure
* :meth:`ixpeobssim.utils.matplotlib_.setup_gca`: setup the current plot
* :class:`ixpeobssim.utils.fmtaxis.fmtaxis`: list of axis labels and units
