.. _sonification:

Sonification
============

Surprisingly enough, ``ixpeobssim`` comes with built-in sonification facilities.
The :ref:`reference-xpsonify` application allows to tranform a given event list
(in the form of a level-2 FITS file) into a midi and/or wave file.

More precisely, the event list is transformed into a `MIDI <https://en.wikipedia.org/wiki/MIDI>`_
file, which can in turn be played with any MIDI sequencer, or tranformed into
an audio file using a virtual (or real, for what it's worth) synthetizer.

.. versionadded:: 19.5.0

   First version of the sonification facilities.

.. versionadded:: 19.6.0

   Minimal support for audio-video animations added.



Manifesto: listening to the X-ray Universe in stereo
----------------------------------------------------

The ideas behind the IXPE-specific sonification facilities are hardly new.
It is customary to transform particles in notes by translating the energy into
a frequency---for photons this is really straightforward, as there is a perfect
analogy between electromagnetic and acoustic waves.

When examined from the standpoint of its spectral response, IXPE is peculiar in
that the standard band-pass for scientific analysis (2--8 keV) corresponds to
exactly 2 octaves, comparable to the typical vocal range of a singer. The
`extended` band in which we define the instrument response (1--12 keV) is in
slight excess of 3.5 octaves, which is pretty much the same as the range of an
acoustic guitar. Phrased in a slightly different way, the IXPE `dynamics` is
comparable with that of typical musical instruments, which makes the conversion
between energy and pitch very natural.

When translating X-ray photons into musical notes, pitch is not the only thing
that we care about. Notes comes with at least three additional properties that
need to be defined:

* velocity: a measure of how vigorously the player initially presses a key
  on the keyboard, to use a piano analogy---the note is louder when the key is
  struck more forcefully (and there are other subtle tonal differences, as well);
* duration: the time difference between the key press and release;
* pan: the location of the sound within the stereo field of the instrument
  (e.g., left, right, center).

The pan is where IXPE really shines. The IXPE polarization-sensitive detectors
are unique in that they not only provide a time, energy and direction
in the sky for each photon---they also provide a proxy of the photon polarization
in the form of the azimuthal angle of emission of the photoelectron. Being an
angle in the plane, this can be mapped in a straightforward way into a point in
a stereo field, i.e., a position in the left-center-right space. This is something
that has essentially no parallel in the vast majority of previously flown X-ray
detectors.

In a sense, IXPE allows us to listen to the Universe in stereo for the first time.


Installation
------------

Working with audio files requires a couple more Python packages on top of the
standard ``ixpeobssim`` dependencies. Since this feature is not of general
interest, these packages are not covered in the installations instructions, and
are not included in the requirement file. If you try and run the sonification
code without the relevant dependencies, you should get sensible error messages
pointing you in the right direction, but essentially everything boils down to
two things:

* `mido <https://mido.readthedocs.io/en/latest/>`_, a MIDI manipulation Python
  library;
* `fluidsynth <https://www.fluidsynth.org/>`_, a cross-platform software synthesizer
  based on the `SoundFont 2 <https://en.wikipedia.org/wiki/SoundFont>`_ specifications.

.. note::

   If you're mot familiar with music production, you can picture a MIDI file
   as a LaTeX source file, and fluidsynth like a typesetting engine producing
   output pdf document ready for visualization. The fonts in this process are
   called _soundfonts_.

Being ``mido`` a pure-Python package, installing it should be as easy as

.. code-block::

   pip install --user mido

For ``fluidsynth`` you should refer to the
`online documentation <https://www.fluidsynth.org/download/>`_. If, e.g.,
you're running a GNU Linux Fedora distribution you are all set:

.. code-block::

   sudo dnf install fluidsynth

In addition, you will also need some soundfonts---after all, fluidsynth is the
rough equivalent of a typesetting system like LaTeX, and you cant't do much
without fonts. If you get an error message sounding like

.. code-block::

   fluidsynth: error: fluid_is_soundfont(): fopen() failed: 'File does not exist.'
   ...
   fluidsynth: error: Failed to load SoundFont "/usr/share/soundfonts/default.sf2"
   Rendering audio to file '/home/lbaldini/ixpeobssimdata/toy_point_source_du1.wav'..
   fluidsynth: warning: No preset found on channel 0 [bank=0 prog=98]

this is almost certainly your problem. You distro might include some soundfont
bank out of the box

.. code-block::

   sudo dnf install fluid-soundfont*

Alternatively, there's plenty of place out there where you can find
General MIDI soundfonts banks---just start from
`here <https://github.com/FluidSynth/fluidsynth/wiki/SoundFont>`_.

The animation facilities are based on the matplotlib
`animation <https://matplotlib.org/stable/api/animation_api.html>`_ module, and you
will need `ffmpeg <http://ffmpeg.org/>`_, a cross-platform, open-source
video recording and editing software. Again, if you're running a GNU Linux
Fedora distribution you are all set:

.. code-block::

   sudo dnf install ffmpeg

(Otherwise refer to the `documentation <http://ffmpeg.org/download.html>`_ online.)


Sonification in details
-----------------------

.. warning::

   This section is a stub.


:ref:`reference-xpsonify` provides a large variety of settings to translate the
properties of the X-rays impinging on the gas pixel detector into notes.

The velocity is determined by the fractional amount of energy released
in the active volume---the photoelectron track is not always fully collected,
and, at any given energy, the amount of energy `collected` by the detector is
a sensible analogous to `how hard the piano player presses the keys`.

The note duration can be related to the track size in the detector or, more
precisely, to the number of pixels in the region of interest for the readout:
the larger the readout window, the longer it takes to the readout electronics
for processing, formatting and transferring out the track.

The pan is determined by the photoelectron azimuthal angle, as explained in the
beginning of the section.

.. note::

   There are others expressive means that the MIDI protocol offers (such as the
   `aftertouch`), and could be incorporated into the sonification process in the
   future, assuming they are not too subtle to be effective.



API
---

Most of the relevant API dealing with sonification are include in the
:mod:`ixpeobssim.evt.animate` and :mod:`ixpeobssim.evt.sonification`
Python modules.
