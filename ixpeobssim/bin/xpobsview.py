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

from __future__ import print_function, division

import numpy

from ixpeobssim.evt.event import xEventFile
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, last_line_color


# pylint: disable=invalid-name

__description__ = \
"""Observation timeline and spacecraft data quick look.

This application retrieves the TIMELINE, GTI and OCTI data from a series of
event lists and plot them as strip charts.
"""

PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()


def _interleave_arrays(a1, a2=None):
    """Interleave two numpy arrays for the purpose of plotting epoch flags.

    See https://stackoverflow.com/questions/5347065/interweaving-two-numpy-arrays
    """
    if a2 is None:
        a2 = a1
    return numpy.vstack((a1, a2)).reshape((-1,), order='F')


def _epochs_to_plottable(data, start_met=None, stop_met=None):
    """Convenience function to turn a series of epocs, in the form of
    START and STOP vectors, into a series of (x, y) values that can be plotted
    as a strip chart.
    """
    met = numpy.sort(numpy.append(data['START'], data['STOP']))
    met = _interleave_arrays(met)
    val = numpy.zeros(met.shape)
    val[1::4] = 1.
    val[2::4] = 1.
    if start_met is not None:
        met = numpy.append(start_met, met)
        val = numpy.append(0, val)
    if stop_met is not None:
        met = numpy.append(met, stop_met)
        val = numpy.append(val, 0)
    return met, val


def _epoch_labels(data, y):
    """Draw the duration of a series of epochs on a strip chart.
    """
    met = 0.5 * (data['START'] + data['STOP'])
    duration = (data['STOP'] - data['START'])
    kwargs = dict(rotation=90., color=last_line_color(), ha='center', va='bottom', size='small')
    for t, dt in zip(met, duration):
        plt.text(t, y, '%.1f s' % dt, **kwargs)


def xpobsview(**kwargs):
    """Main function.
    """
    # Vertical spacing between adjacent strip charts.
    spacing = 0.2
    # Horrible hack to keep track of how many strip charts we have plotted, and
    # at which y coordinate we should put the next one.
    global index
    index = 0

    def yoffset():
        """Nested function calculating the current y offset.
        """
        global index
        val = index + spacing * (index + 1)
        index += 1
        return val

    plt.figure('Observation timeline')
    for i, file_path in enumerate(kwargs.get('filelist')):
        event_file = xEventFile(file_path)
        start_met = event_file.start_met()
        stop_met = event_file.stop_met()
        # For the first file only we plot the stuff that' s identical across
        # differnt DUs, i.e., the timeline and the GTIs.
        if i == 0:
            # Timeline data...
            timeline_data = event_file.timeline_data()
            if timeline_data is not None:
                met = _interleave_arrays(timeline_data['START'], timeline_data['STOP'])
                in_saa = _interleave_arrays(timeline_data['IN_SAA'])
                occult = _interleave_arrays(timeline_data['TARGET_OCCULT'])
                plt.plot(met, in_saa + yoffset(), label='SAA')
                for t in timeline_data['START']:
                    plt.axvline(t, ls='dashed', color='gray', lw=0.75)
                plt.plot(met, occult + yoffset(), label='Target occulted by Earth')
            # GTI data...
            gti_data = event_file.gti_data()
            met, gti = _epochs_to_plottable(gti_data, start_met, stop_met)
            yoff = yoffset()
            plt.plot(met, gti + yoff, label='Good time intervals')
            _epoch_labels(gti_data, yoff)
        # OCTI data...
        octi_data = event_file.octi_data()
        if octi_data is not None:
            met, octi = _epochs_to_plottable(octi_data, start_met, stop_met)
            du_id = event_file.du_id()
            yoff = yoffset()
            label = 'On-orbit calibration time intervals for DU %d' % du_id
            plt.plot(met, octi + yoff, label=label)
            _epoch_labels(octi_data, yoff)

        ymax = index * (1. + 2 * spacing) + 3.
        setup_gca(xmin=start_met, xmax=stop_met, legend=True, ymax=ymax,
                  xlabel='MET [s]', ylabel='Timeline flags')
        plt.yticks([])


def main():
    xpobsview(**PARSER.parse_args().__dict__)
    plt.show()



if __name__ == '__main__':
    main()
