#!/usr/bin/env python
#
# Copyright (C) 2021, the ixpeobssim team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""Sonification utilities.
"""

from __future__ import print_function, division

from functools import total_ordering
import subprocess
import wave

import numpy

from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.matplotlib_ import plt

try:
    from mido import Message, MetaMessage, MidiFile, MidiTrack, second2tick, bpm2tempo
    MIDO_INSTALLED = True
except ImportError as e:
    logger.error(e)
    logger.info('mido is a library to manipulate MIDI files.')
    logger.info('See https://mido.readthedocs.io/en/latest/ for more details.')
    logger.info('Type `pip install --user mido` to install it')
    MIDO_INSTALLED = False
    # This is needed for the xMidiFile class below not to chocke on the inheritance tree.
    MidiFile = object


CENTRAL_A_FREQ = 440.
NOTES_PER_OCTAVE = 12

# pylint: disable=invalid-name, too-few-public-methods, too-many-locals

class xMidiNote:

    """Small class encapsulating a MIDI note.

    There's an infinite number of place on the web where one can find conversion
    tables between MIDI notes and physical characteristics, see, e.g.,
    https://www.inspiredacoustics.com/en/MIDI_note_numbers_and_center_frequencies

    A MIDI note is univiquely identified by a note number, ranging from 0 to 127.
    Here we shall restrict ourselves to the piano keys, i.e., note numbers from
    21 to 127 (the lower notes are too low to be useful, anyway), with the
    understanding that:

    * note 21 is A0, at 27.500 Hz;
    * note 127 is G9, at 12543.854 Hz;
    * note 69 is A4, at 440.000 Hz.

    A MIDI note is initialized by its MIDI note number, and encapsulates all the
    facilities to calculate the frequency, note name and alike.
    """

    NOTE_NAMES = ['A', 'A#', 'B', 'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#']
    NOTE_OFFSET = 21
    CENTRAL_A_NUMBER = 69

    def __init__(self, note_number):
        """Constructor.
        """
        assert self.NOTE_OFFSET <= note_number < 128
        self.note_number = note_number

    def frequency(self):
        """Return the note frequency.
        """
        return CENTRAL_A_FREQ * 2.**((self.note_number - self.CENTRAL_A_NUMBER) / 12.)

    def note_name(self):
        """Return the note name.
        """
        return self.NOTE_NAMES[(self.note_number - self.NOTE_OFFSET) % NOTES_PER_OCTAVE]

    def octave(self):
        """Return the note octave.
        """
        return (self.note_number - 12) // NOTES_PER_OCTAVE

    def name(self):
        """Return the note name and octave.
        """
        return '%s%d' % (self.note_name(), self.octave())

    def __str__(self):
        """String formatting.
        """
        return '%3d - %3s (%.3f Hz)' % (self.note_number, self.name(), self.frequency())



class xMusicalScale:

    """Simple class describing a musical scale.

    A scale is essentilly given by a series (in the form of a numpy array) of
    note number covering the piano dynamic range.
    """

    SCALE_DICT = {
        'Chromatic': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        'Ionian': [1, 3, 5, 8, 10],
        'Dorian': [1, 3, 6, 8, 11],
        'Phrygian': [1, 4, 6, 9, 11],
        'Myxolydian': [1, 3, 6, 8, 10],
        'Aeolian': [1, 4, 6, 8, 11],
    }

    OFFSET_DICT = {'C': 11}

    def __init__(self, mode='Ionian', key='C'):
        """Constructor.
        """
        self.mode = mode
        self.key = key
        base = numpy.array(self.SCALE_DICT[mode])
        offset = self.OFFSET_DICT[key]
        self.notes = numpy.array([], dtype=int)
        # Compose the entire scale.
        for octave in range(11):
            self.notes = numpy.append(self.notes, base + offset + octave * NOTES_PER_OCTAVE)
        # Trim the notes outside the piano dynamic range.
        mask = numpy.logical_and(self.notes >= xMidiNote.NOTE_OFFSET, self.notes < 128)
        self.notes = self.notes[mask]

    def snap(self, values):
        """Snap a series of notes to the scale.
        """
        idx = numpy.searchsorted(self.notes, values, side='left')
        mask = values - self.notes[idx - 1] <= self.notes[idx] - values
        idx = (idx - 1) * mask + idx * numpy.logical_not(mask)
        return self.notes[idx]

    def __str__(self):
        """String formatting.
        """
        text = '%s scale, key = %s\n' % (self.mode, self.key)
        text += '\n'.join([str(xMidiNote(num)) for num in self.notes])
        return text



@total_ordering
class xMidiEvent:

    """Small classe describing a midi event.

    This is intended as a temporary storage for actual MIDI messages expressed in
    absolute time and not necessarily time-ordered.
    """

    def __init__(self, type_, timestamp, **kwargs):
        """Constructor.
        """
        self.type_ = type_
        self.timestamp = timestamp
        self.kwargs = kwargs

    def __eq__(self, other):
        """Overloaded operator.
        """
        return self.timestamp == other.timestamp

    def __lt__(self, other):
        """Overloaded operator.
        """
        return self.timestamp < other.timestamp

    def __str__(self):
        """String formatting.
        """
        return '%s @ %.3f (%s)' % (self.type_, self.timestamp, self.kwargs)



class ContolChangeParameter:

    """Most common contol change codes, see
    https://professionalcomposers.com/midi-cc-list/
    """

    MODULATION_WHEEL = 1
    BREATH_CONTROL = 2
    VOLUME = 7
    PAN = 10
    EXPRESSION = 11
    SUSTAIN_PEDAL = 64
    PORTAMENTO = 65
    RESONANCE = 71
    FREQUENCY_CUTOFF = 74



class ProgramChangePrograms:

    """Program change programs, see
    https://www.recordingblogs.com/wiki/midi-program-change-message
    """

    ACOUSTIC_GRAND_PIANO = 0
    BRIGHT_ACOUSTIC_GRAND_PIANO = 1
    ELECTRIC_GRAND_PIANO = 2
    HONKY_TONK_PIANO = 3
    ELECTRIC_PIANO_1 = 4
    ELECTRIC_PIANO_2 = 5
    HARPSICORD = 6
    CLAVINET = 7
    CELESTA = 8
    GLOKENSPIEL = 9
    MUSIC_BOX = 10
    VIBRAPHONE = 11
    MARIMBA = 12
    XYLOPHONE = 13
    TUBOLAR_BELL = 14
    DULCIMER = 15
    HAMMOND_ORGAN = 16
    PERCUSSIVE_ORGAN = 17
    ROCK_ORGAN = 18
    CHURCH_ORGAN = 19
    REED_ORGAN = 20
    ACCORDION = 21
    HARMONICA = 22
    TANGO_ACCORDION = 23
    NYLON_STRING_ACOUSTIC_GUITAR = 24
    STEEL_STRING_ACOUSTIC_GUITAR = 25
    JAZZ_ELECTRIC_GUITAR = 26
    CLEAN_ELECTRIC_GUITAR = 27
    MUTED_ELECTRIC_GUITAR = 28
    OVERDRIVEN_GUITAR = 29
    DISTORTION_GUITAR = 30
    GUITAR_HARMONICS = 31
    ACOUSTIC_BASS = 32
    FINGERED_ELECTRIC_BASS = 33
    PICKED_ELECTRIC_BASS = 34
    FRETLESS_BASS = 35
    SLAP_BASS_1 = 36
    SLAP_BASS_2 = 37
    SYNTH_BASS_1 = 38
    SYNTH_BASS_2 = 39
    VIOLIN = 40
    VIOLA = 41
    CELLO = 42
    CONTRABASS = 43
    TREMOLO_STRINGS = 44
    PIZZICATO_STRINGS = 45
    ORCHESTRAL_STRINGS = 46
    TIMPANI = 47
    STRING_ENSEMBLE_1 = 48
    STRING_ENSEMBLE_2 = 49
    SYNTH_STRINGS_1 = 50
    SYNTH_STRINGS_2 = 51
    CHOIR_AAHS = 52
    VOICE_OOHS = 53
    SYNTH_CHOIR = 54
    ORCHESTRA_HIT = 55
    TRUMPET = 56
    TROMBONE = 57
    TUBA = 58
    MUTED_TRUMPET = 59
    FRENCH_HORN = 60
    BRASS_ENSEMBLE = 61
    SYNTH_BRASS_1 = 62
    SYNTH_BRASS_2 = 63
    SOPRANO_SAX = 64
    ALTO_SAX = 65
    TENOR_SAX = 66
    BARITONE_SAX = 67
    OBOE = 68
    ENGLISH_HORN = 69
    BASSOON = 70
    CLARINET = 71
    PICCOLO = 72
    FLUTE = 73
    RECORDER = 74
    PAN_FLUTE = 75
    BOTTLE_BLOW = 76
    SHAKUHACHI = 77
    WHISTLE = 78
    OCARINA = 79
    SYNTH_SQUARE_WAVE = 80
    SYNTH_SAW_WAVE = 81
    SYNTH_CALLIOPE = 82
    SYNTH_CHIFF = 83
    SYNTH_CHARANG = 84
    SYNTH_VOICE = 85
    SYNTH_FIFTHS_SAW = 86
    SYNTH_BRASS = 87
    FANTASIA = 88
    WARM_PAD = 89
    POLYSYNTH = 90
    SPACE_VOX = 91
    BOWED_GLASS = 92
    METAL_PAD = 93
    HALO_PAD = 94
    SWEEP_PAD = 95
    ICE_RAIN = 96
    SOUNDTRACK = 97
    CRYSTAL = 98
    ATMOSPHERE = 99
    BRIGHTNESS = 100
    GOBLINS = 101
    ECHO_DROPS = 102
    SCI_FI = 103
    SITAR = 104
    BANJO = 105
    SHAMISEN = 106
    KOTO = 107
    KALIMBA = 108
    BAG_PIPE = 109
    FIDDLE = 110
    SHANAI = 111
    TINKLE_BELL = 112
    AGOGO = 113
    STEEL_DRUMS = 114
    WOODBLOCK = 115
    TAIKO_DRUM = 116
    MELODIC_TOM = 117
    SYNTH_DRUM = 118
    REVERSE_CYMBAL = 119
    GUITAR_FRET_NOISE = 120
    BREATH_NOISE = 121
    SEASHORE = 122
    BIRD_TWEET = 123
    TELEPHONE_RING = 124
    HELICOPTER = 125
    APPLAUSE = 126
    GUNSHOT = 127



class xMidiFile(MidiFile):

    """Small wrapper around the midio.MidiFile class.
    """

    def __init__(self, ticks_per_beat=480):
        """Constructor.
        """
        MidiFile.__init__(self, ticks_per_beat=ticks_per_beat)

    @staticmethod
    def _load_event_data(file_path, max_num_events, speed=1., roi=None):
        """Load the event data from the input level 2 FITS file.

        Note the FITS file is opened here and automatically closed when the
        the file object goes out of scope.
        """
        input_file = xEventFile(file_path)
        max_num_events = min(max_num_events, len(input_file.event_data))
        cols = [
            input_file.energy_data(mc=True),
            input_file.energy_data(mc=False),
            (input_file.time_data() - input_file.min_good_time()) / speed,
            input_file.phi_data()
        ]
        num_events = len(cols[0])
        logger.info('%d event(s) loaded.', num_events)
        if roi is not None:
            t = cols[2]
            ra, dec = input_file.sky_position_data(mc=False)
            ra0, dec0 = input_file.wcs_reference()
            mask = roi.event_mask(t, ra, dec, ra0, dec0) * (t <= roi.tmax)
            logger.info('Applying event mask....')
            cols = [col[mask] for col in cols]
            num_events = len(cols[0])
            logger.info('%d event(s) remaining.', num_events)
        if num_events > max_num_events:
            cols = [col[:max_num_events] for col in cols]
        return cols

    @staticmethod
    def compute_note_number(energy, mode='Ionian', key='C', pivot=4., dynamic_scaling=1.):
        """Convert an energy (in keV) into the corresponding note number.

        This is accomplished assigninig the A4 note number to the pivot energy, and
        scaling the energies in octave space. Since the extended IXPE bandpass
        over which we calculate the response functions (1--15 keV) corresponds to
        roughly 4 octaves, this is not a terrible match for a sonification project.
        The nominal IXPE energy band (2--8 keV) is only 2 octaves, and in real life
        we might need to tweak things a little bit, which is the purpose of the
        dynamic_scaling argument (see the comments in the code for a detailed
        explanation of the algorithm).

        Args
        ----
        energy : array_like
            The photon energy.

        mode : str
            The name of the musical scale to which the energy values should be snapped.

        key : str
            The key for the aforementioned musical scale.

        pivot : float
            The pivot energy (mapped to A4).

        dynamic_scaling : float
            Empirical parameter for scaling the energy in octave spaces.
        """
        # Convert the energy to the relative position in octave space. Note at this
        # point the pivot is mapped to (floating point) 0.
        notes = NOTES_PER_OCTAVE * numpy.log2(energy / pivot)
        # Multiplicative dynamic scaling. Note at this point the pivot is *still*
        # mapped to (floating point) 0.
        notes *= dynamic_scaling
        # Add the offset so that the pivot is mapped to A4.
        notes += xMidiNote.CENTRAL_A_NUMBER
        # Snap to proper musical scale.
        return xMusicalScale(mode, key).snap(notes)

    @staticmethod
    def compute_velocity(mc_energy, rec_energy, pivot=64., dynamic_scaling=1.):
        """Convert the energy measurement to note velocity.

        The note velocity is determined by the ratio between the reconstructed and
        true energy, i.e., events in the left tail of the energy dispersion are
        rendered with a lower velocity.

        Args
        ----
        mc_energy : array_like
            The true energy.

        rec_energy : array_like
            The reconstructed energy.

        pivot_velocity : float
            The pivot value for the note velocity, assigned when the reconstructed
            energy is equal to the true energy.

        dynamic_scaling : float
            Empirical parameter to enhance the veocity dynamics.
        """
        delta = dynamic_scaling * (rec_energy / mc_energy - 1.)
        velocity = pivot * (delta + 1.)
        # Clip the output values to the valid range and cast to an integer.
        velocity = numpy.clip(numpy.rint(velocity), 0., 127.).astype(int)
        return velocity

    @staticmethod
    def compute_pan(phi, dynamic_scaling=1.):
        """Compute the pan.

        This is simply mapping the photoelectron angle into the 0--127 phisical
        range.

        Args
        ----
        phi : array_like
            The photoelectron angle.

        dynamic_scaling : float
            Empirical parameter to enhance the veocity dynamics.
        """
        # Map phi from [-pi, pi] into [-0.5, 0.5]
        pan = phi / (2. * numpy.pi)
        # Apply the dynamic scaling.
        pan = numpy.clip(pan * dynamic_scaling, -0.5, 0.5)
        # Convert from [-0.5, 0.5] to [0., 127.]
        pan = 127 * (pan + 0.5)
        # Cast to integer.
        pan = numpy.rint(pan).astype(int)
        return pan

    @staticmethod
    def track_name_meta_message(name, time=0):
        """Return a track_name MIDI meta-message.
        """
        return MetaMessage('track_name', name=name, time=time)

    @staticmethod
    def set_tempo_meta_message(bpm, time=0):
        """Return a set_tempo MIDI meta-message.
        """
        return MetaMessage('set_tempo', tempo=bpm2tempo(bpm), time=time)

    @staticmethod
    def program_change_message(program, channel=0):
        """Return a program change MIDI message.
        """
        return Message('program_change', program=program, channel=channel)

    @staticmethod
    def control_change_message(control, value, channel=0):
        """Return a control change message.
        """
        return Message('control_change', control=control, value=value, channel=channel)

    @staticmethod
    def pan_message(value, channel=0):
        """Return a pan message.
        """
        return xMidiFile.control_change_message(ContolChangeParameter.PAN, value, channel)

    def add_midi_track(self, track_name=None, bpm=None):
        """Add a track to the MIDI file.
        """
        track = MidiTrack()
        if track_name is not None:
            track.append(self.track_name_meta_message(track_name))
        if bpm is not None:
            track.append(self.set_tempo_meta_message(bpm))
        self.tracks.append(track)
        return track

    def fill(self, file_path, roi=None, **kwargs):
        """Fill the MIDI file with the event data in the proper format.

        This is the main function where most of the action actually happens.
        """
        event_list = []
        # Read the actual event data.
        args = file_path, kwargs.get('maxnotes'), kwargs.get('speed'), roi
        mc_energy, rec_energy, timestamp, phi = self._load_event_data(*args)
        # Calculate the notes.
        note_kwargs = dict(mode=kwargs.get('scale'), key=kwargs.get('key'),
            pivot=kwargs.get('epivot'), dynamic_scaling=kwargs.get('escaling'))
        logger.info('Note settings: %s', note_kwargs)
        note = self.compute_note_number(mc_energy, **note_kwargs)
        # Calculate the velocity.
        velocity_kwargs = dict(pivot=kwargs.get('vpivot'), dynamic_scaling=kwargs.get('vscaling'))
        logger.info('Velocity settings: %s', velocity_kwargs)
        velocity = self.compute_velocity(mc_energy, rec_energy, **velocity_kwargs)
        # Calculate the pan.
        pan_kwargs = dict(dynamic_scaling=kwargs.get('pscaling'))
        logger.info('Pan settings: %s', pan_kwargs)
        pan = self.compute_pan(phi, **pan_kwargs)
        # Calculate the note duration.
        duration_kwargs = dict(duration=kwargs.get('duration'))
        logger.info('Duration settings: %s', duration_kwargs)
        duration = numpy.full(phi.shape, duration_kwargs.get('duration'))
        for t, n, v, d, p in zip(timestamp, note, velocity, duration, pan):
            _pan = xMidiEvent('control_change', t, control=ContolChangeParameter.PAN, value=p)
            _note_on = xMidiEvent('note_on', t, note=n, velocity=v)
            _note_off = xMidiEvent('note_off', t + d, note=n, velocity=v)
            event_list += [_pan, _note_on, _note_off]
        event_list.sort()
        # Add the main track.
        bpm = kwargs.get('bpm')
        track = self.add_midi_track('source_name', bpm)
        track.append(self.program_change_message(kwargs.get('program')))
        timestamp = numpy.array([event.timestamp for event in event_list])
        delta_time = numpy.diff(timestamp, prepend=0.)
        delta_time = second2tick(delta_time, self.ticks_per_beat, bpm2tempo(bpm))
        delta_time = numpy.rint(delta_time).astype(int)
        for event, dt in zip(event_list, delta_time):
            track.append(Message(event.type_, time=dt, **event.kwargs))
        # If needed, add the pad track.
        if kwargs.get('pad'):
            ch = kwargs.get('padchannel')
            nn = kwargs.get('padnote')
            vel = kwargs.get('padvelocity')
            track = self.add_midi_track('background', bpm)
            track.append(self.program_change_message(95, channel=ch))
            _note_on = xMidiEvent('note_on', 0, note=nn, velocity=vel, channel=ch)
            _note_on = xMidiEvent('note_on', timestamp.max() - 10., note=nn, velocity=vel, channel=ch)
            track.append(Message(_note_on.type_, time=0, **_note_on.kwargs))
        xMidiFile.plot_diagnostics(timestamp, note, velocity, pan, duration)

    @staticmethod
    def _plot_diagnostics_base(title, values, binning=None):
        """Create a simple diagnostic plot.
        """
        plt.figure(title)
        if binning is None:
            binning = numpy.linspace(-0.5, 127.5, 129)
        xHistogram1d(binning, xlabel=title).fill(values).plot()

    @staticmethod
    def plot_diagnostics(timestamp, note, velocity, pan, duration):
        """Create all the diagnostics plots.
        """
        xMidiFile._plot_diagnostics_base('Timestamp', timestamp, numpy.linspace(0., timestamp.max(), 100))
        xMidiFile._plot_diagnostics_base('Note number', note)
        xMidiFile._plot_diagnostics_base('Velocity', velocity)
        xMidiFile._plot_diagnostics_base('Pan', pan)
        xMidiFile._plot_diagnostics_base('Duration', duration, numpy.linspace(0., 2., 100))



def _fluidsynth_call(midi_file_path, sound_font, arg_list):
    """Generic fluidsynth call to operate on MIDI streams.

    This is just a plain subprocess call that can be specialized, e.g., to
    convert a MIDI file to wave or to play a MIDI file.
    """
    assert midi_file_path.endswith('.mid')
    args = ['fluidsynth'] + arg_list
    logger.info('About to execute %s', args)
    try:
        subprocess.call(args)
    except FileNotFoundError as e:
        logger.error(e)
        logger.info('You need fluidsynth to convert from MIDI to wave.')
        logger.info('See https://www.fluidsynth.org/download/ for more details')
        abort()


def midi_to_wav(midi_file_path, sound_font=''):
    """Convert a MIDI file to audio using fluidsynth.
    """
    logger.info('Converting %s to audio...', midi_file_path)
    wav_file_path = midi_file_path.replace('.mid', '.wav')
    arg_list = ['-ni', sound_font, midi_file_path, '-F', wav_file_path]
    _fluidsynth_call(midi_file_path, sound_font, arg_list)
    return wav_file_path


def play_midi(midi_file_path, sound_font=''):
    """Play a MIDI file file through fluidsynth.
    """
    logger.info('Playing MIDI file %s...', midi_file_path)
    arg_list = ['-i', sound_font, midi_file_path]
    _fluidsynth_call(midi_file_path, sound_font, arg_list)


def stereo_to_mono(file_path):
    """Convert a stereo wave file to a mono file containing the average of the
    two channels.
    """
    assert file_path.endswith('.wav')
    logger.info('Converting %s from stereo to mono...', file_path)
    input_stream = wave.open(file_path, 'rb')
    if input_stream.getnchannels() != 2:
        logger.warning('Not a stereo file, returning...')
        return
    logger.info('Reading input data...')
    num_frames = input_stream.getnframes()
    data = numpy.frombuffer(input_stream.readframes(num_frames), dtype=numpy.int16)
    logger.info('Averaging...')
    data = data.astype(numpy.int32).reshape((num_frames, 2))
    average = data.sum(axis=1) / 2
    data[:, 0] = average
    data[:, 1] = average
    data = data.astype(numpy.int16)
    output_file_path = file_path.replace('.wav', '_mono.wav')
    logger.info('Writing output file %s...', output_file_path)
    output_stream = wave.open(output_file_path, 'wb')
    output_stream.setparams(input_stream.getparams())
    output_stream.writeframes(data.tobytes())
    output_stream.close()
