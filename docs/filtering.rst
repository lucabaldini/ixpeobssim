.. _filtering:

Filtering event lists
=====================

``ixpeobssim`` provides facilities to filter event lists in different flavors,
selecting in time, phase, energy and position in the sky. Or work-horse
application in this respect is :ref:`reference-xpselect`.


.. program-output:: export PYTHONPATH=../:$PYTHONPATH; python ../ixpeobssim/bin/xpselect.py --help
   :ellipsis: 0,11
   :shell:

.. note::

   .. versionadded:: 22.0.0

   :ref:`reference-xpselect` was largely refactored in ixpeobssim version
   22.0.0, mainly in response to
   `issue #169 <https://github.com/lucabaldini/ixpeobssim/issues/169>`_.
   If you are using ``xpselect`` and lagging behind the ``ixpeobssim`` version,
   you should definitely consider upgrading.

   At the time of writing, the `XSELECT <https://heasarc.gsfc.nasa.gov/ftools/xselect/>`_
   FTOOL is being updated to support IXPE, with a new
   `HEASOFT <https://heasarc.gsfc.nasa.gov/docs/software/heasoft/>`_ release
   expected in February, 2022, and it is conceivable that :ref:`reference-xpselect`
   will be modified in the future to avoid overlap in functionality with the
   tools officially supported by `HEASARC <https://heasarc.gsfc.nasa.gov/>`_.
   In a nutshell: we shall try and strip off xpselect anything that will be
   provided out off the shelf by XSELECT.

The basic idea is that `filtered event lists should be functionally identical to
their parents level-2 files` and they should, at least in principle, inter-operate
with all the analysis tools in exactly the same fashion. If that is not the case,
it is likely a bug :-)


Selecting in time or phase
--------------------------

Selections in time or phase are achieved by means of the ``--tmin/max`` and
``--phasemin/max`` command-line switches. For the sake of making it easier to
keep track of the livetime and time-related header keywords when selecting in
time and phase, ``xpselect`` performs some minimal validation on the input
values of the command line switches, and will exit with an error message if
any of the conditions below is not met:

* ``TSTART <= TMIN <= TSTOP``
* ``TSTART <= TMAX <= TSTOP``
* ``TMAX > TMIN``
* ``0. <= PHASEMIN <= 1.``
* ``0. <= PHASEMAX <= 1.``

.. warning::

   :ref:`reference-xpselect` does not currently support concurrent selections
   on time `and` phase, and there is no specific plan to add this functionality,
   unless proved to be important in practice. In other words, you can select
   an arbitrary time `or` phase interval, but not at the same time.

   Chaining multiple :ref:`reference-xpselect` calls in time or phase is also
   discouraged as, depending on the exact setup, this makes keeping track of the
   livetime extremely difficult. This should not be a limitation for any real
   analysis.

It goes without saying that, in order to be able to select in phase, you must
run the :ref:`reference-xpphase` tool first, in order to have the ``PHASE``
column added to the event list.


Propagating the livetime
~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

   .. versionadded:: 22.0.0

   As of version :ref:`reference-xpselect` provides optional command-line switches to
   propagate the relevant time-related header keywords to the filtered event lists.
   We are fully aware that the treatment is approximate, at the limit of
   naiveness, and we fully anticipate it will undergo modifications as we
   gain experience with flight data.

.. warning::

   The new infrastructure for the livetime correction is only suitable to experiment
   with simulated data, and does not yet inter-operate with filtered level-2
   flight data. This is the reason why the correction is disabled by default.


:ref:`reference-xpselect` tries to do a decent job in keeping track of the
livetime when selecting in time or phase. Unfortunately the matter is not completely
trivial, and a short discussion is in order, here.

IXPE event lists provide a ``LIVETIME`` column encapsulating
`the livetime from the previous event in microseconds`, or the total time
in which the trigger was enabled from the end of the readout of the previous
event from the trigger time of the current event. Under normal conditions the
sum of the ``LIVETIME`` column should be equal to the values of the ``LIVETIME``
and ``EXPOSURE`` header keywords up to numerical rounding errors.

The default mean for :ref:`reference-xpselect` to re-calculate the livetime, when
time or phase selections are applied is to simply sum the ``LIVETIME`` column
`over the events surviving the selections`. This is what happens when the
``--ltimealg`` command-line switch is set to ``LTSUM``, and the thing works as
advertised when the selected events are contiguous, e.g., when a simple time
selections is applied.

When selecting in phase the thing gets more complicated, as in some circumstances
(i.e., when the average time distance between subsequent events is larger than
the period of the ephemeris we use for folding) the ``LIVETIME`` value for any
particular event will in general be referred to an event in a `different phase bin`.
In order to cope (approximately) with this situation, :ref:`reference-xpselect`
provide an alternative mean of propagating the livetime, that is enabled by
using the ``--ltimealg LTSCALE`` command-line switch. In a nutshell, we first
calculate the average dead time per event for the original event list by simply
dividing the difference between the ``ONTIME`` and ``LIVETIME`` header keywords
by the total number of events. (This is typically in the ms ball-park.) This
average dead time per event can in turn be used to calculate a total `scaled livetime`
within the selection, by simply multiplying it by the number of events within
the selection.

.. warning::

   The livetime scaled we just described is only valid in the limit that the
   average dead time per event in the original photon list is representative
   of the average dead time per event within the selection. Since the dead time
   per event depends on the size of the region of interest of the event, which
   in turns depends on the energy, this assumption might not be accurate in the
   presence of strong spectral variations as a function of the pulse phase.
   (In practice, all the dependences are mild, and this should not be an issue
   in most practical situations.)

In many conditions, and particularly for selections in time, the default
``LTSUM`` mechanism will work fine. As a rule of thumb, you can use ``LTSUM``
when selecting in time and ``LTSCALE`` when selecting in phase, but be advised
that there might be cases when one need to examine more closely the situation and
opt for an ad-hoc solution.

For completeness, when we select in time, the header keywords in the releavant
extensions are update as follows

* ``TSTART -> max(TSTART, TMIN)``
* ``TSTOP -> min(TSTOP, TMAX)``
* ``ONTIME -> TSTOP - TSTART``

and when we select in phase we have instead.

* ``TSTART -> TSTART (unchanged)``
* ``TSTOP -> TSTOP (unchanged)``
* ``ONTIME -> `ONTIME * selected_phase_fraction``

In all cases, the livetime-related columns are updated according to:

* ``EXPOSURE -> modified livetime (LTSUM or LTSCALE)``
* ``LIVETIME -> same as EXPOSURE``
* ``DEADC -> EXPOSURE / ONTIME``
