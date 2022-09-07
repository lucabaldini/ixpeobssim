#!/usr/bin/env python
#
# Copyright (C) 2021 the ixpeobssim team.
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

from __future__ import print_function, division

import os

from matplotlib.animation import FFMpegWriter, PillowWriter

from ixpeobssim import IXPEOBSSIM_CONFIG_FITS
from ixpeobssim.evt.animate import xMovingCircle, xSkyAnimation
from ixpeobssim.evt.sonify import xMidiNote, xMidiFile, xMusicalScale,\
    midi_to_wav, play_midi, stereo_to_mono, MIDO_INSTALLED
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger, abort
from ixpeobssim.utils.matplotlib_ import plt
from ixpeobssim.utils.system_ import import_module


__description__ = \
"""Sonify an existing level-2 file.

This is complex application in the beta-testing stage that allows to convert
a photon list into notes, with optional support for animated imagery.

At the fundamental level the application converts a list of photons into a MIDI
file where each note corresponds to an event, based on its energy, position and
polarization. The command-line switches allow to customize the conversion in
many different ways.

There are options to convert the MIDI file into an actual .wav audio file, in
both stereo and mono flavors.

Additionally, the application takes an optional input FIT image and a
time-dependent selection ROI to create an animation to accompany the music.
(Note that you will have to merge the two manually after the fact.)

See the ixpeobssim documentation for installation and usage.
"""

PARSER = xArgumentParser(description=__description__)
PARSER.add_file()
PARSER.add_outfile()
PARSER.add_argument('--maxnotes', type=int, default=1000,
    help='maximum number of notes in the output MIDI file')
PARSER.add_argument('--speed', type=float, default=1.,
    help='scaling factor for the event timing in the output file')
PARSER.add_argument('--scale', type=str, default='Ionian', choices=xMusicalScale.SCALE_DICT,
    help='the musical scale to be used')
PARSER.add_argument('--key', type=str, default='C', choices=xMidiNote.NOTE_NAMES,
    help='the key to be used')
PARSER.add_argument('--epivot', type=float, default=4.,
    help='energy pivot in keV')
PARSER.add_argument('--escaling', type=float, default=1.,
    help='energy dynamic scaling')
PARSER.add_argument('--vpivot', type=float, default=64.,
    help='velocity pivot')
PARSER.add_argument('--vscaling', type=float, default=1.,
    help='velocity dynamic scaling')
PARSER.add_argument('--pscaling', type=float, default=4.,
    help='pan dynamic scaling')
PARSER.add_argument('--duration', type=float, default=0.25,
    help='note duration')
PARSER.add_argument('--bpm', type=float, default=60.,
    help='track tempo in beat per minutes')
PARSER.add_argument('--program', type=int, default=99,
    help='track instrument program to be used')
PARSER.add_boolean('--pad', default=True,
    help='enable the background pad')
PARSER.add_argument('--padprogram', type=int, default=95,
    help='program for the the background pad')
PARSER.add_argument('--padnote', type=int, default=48,
    help='note for the the background pad')
PARSER.add_argument('--padvelocity', type=int, default=36,
    help='the velocity for the the background pad')
PARSER.add_argument('--padchannel', type=int, default=1,
    help='the MIDI channel for the background pad')
PARSER.add_argument('--roiconfig', type=str, default=None,
    help='path to the ROI configuration file')
PARSER.add_boolean('--wav', default=True,
    help='convert the midi file to wav after the fact')
PARSER.add_boolean('--mono', default=False,
    help='create a mono copy of the output wave file')
PARSER.add_argument('--animimage', type=str, default=None,
    help='path to the still FITS image for the animation')
PARSER.add_argument('--animtype', type=str, default='mp4', choices=('mp4', 'avi', 'mov'),
    help='the file type for the animation output file')
PARSER.add_argument('--animfps', type=int, default=25,
    help='animation speed in frame per seconds')
PARSER.add_boolean('--playback', default=False,
    help='play the MIDI file')
PARSER.add_boolean('--verbose', default=False,
    help='print the midi stream to the terminal')
PARSER.add_argument('--soundfont', type=str, default='',
    help='path to the soundfont for the audio conversion')
PARSER.add_boolean('--diagnostics', default=False,
    help='show the diagnostic plots')



def xpsonify(**kwargs):
    """Sonify.
    """
    # Cache the necessary keyword arguments.
    input_file_path = kwargs.get('file')
    output_file_path = kwargs.get('outfile')
    config_file_path = kwargs.get('roiconfig')
    img_file_path = kwargs.get('animimage')
    if output_file_path is None:
        output_file_path = input_file_path.replace('.fits', '.mid')
    # Prepare the (moving) ROI, if needed.
    if config_file_path is not None:
        logger.info('Loading ROI from %s...', config_file_path)
        roi = import_module(config_file_path).ROI
    else:
        roi = None
    # Create the MIDI file.
    mid = xMidiFile()
    mid.fill(input_file_path, roi, **kwargs)
    if kwargs.get('verbose'):
        print(mid)
    logger.info('Saving output MIDI file to %s...', output_file_path)
    mid.save(output_file_path)
    # Convert to wave,
    if kwargs.get('wav'):
        wav_file_path = midi_to_wav(output_file_path, kwargs.get('soundfont'))
        # And, possibly, in mono too.
        if kwargs.get('mono'):
            stereo_to_mono(wav_file_path)
    # Make the animation.
    if img_file_path:
        file_type = kwargs.get('animtype')
        fps = kwargs.get('animfps')
        logger.info('Animating %s...', img_file_path)
        img = xSkyAnimation(img_file_path)
        interval = int(round(1000 / fps))
        logger.debug('Animation interval set to %d ms', interval)
        anim = img.run(roi, interval)
        video_file_path = output_file_path.replace('.mid', '.%s' % file_type)
        logger.info('Saving video to %s...', video_file_path)
        anim.save(video_file_path, FFMpegWriter(fps))
    if kwargs.get('playback'):
        play_midi(output_file_path, kwargs.get('soundfont'))
    logger.info('All done, bye :-)')
    if kwargs.get('diagnostics'):
        plt.show()


def main():
    """
    """
    if not MIDO_INSTALLED:
        abort()
    xpsonify(**PARSER.parse_args().__dict__)


if __name__ == '__main__':
    main()
