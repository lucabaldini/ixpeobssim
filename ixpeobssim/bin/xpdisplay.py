#!/usr/bin/env python
#
# Copyright (C) 2022, the ixpeobssim team.
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

"""xpdisplay app.
"""

from __future__ import print_function, division

from astropy.io import fits
import numpy

from ixpeobssim.evt.display import xL1EventFile, xXpolGrid
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, xStatBox


__description__ = \
"""IXPE single event display.

This application provides visualization support for track images contained in
Level-1 IXPE files.
"""

PARSER = xArgumentParser(description=__description__)
PARSER.add_file()
PARSER.add_argument('--evtlist', type=str,
    help='path to the auxiliary (Level-2 file) event list')



def _run_display(file_path, **kwargs):
    """Run the event display.
    """
    grid = xXpolGrid()
    event_file = xL1EventFile(file_path)
    threshold = event_file.zero_sup_threshold()
    logger.info('Zero suppression threshold: %d', threshold)
    if kwargs.get('evtlist'):
        event_list = xEventFile(kwargs.get('evtlist'))
        event_data = event_list.time_data(), event_list.energy_data(),\
            *event_list.sky_position_data(), *event_list.stokes_data()
        for met, energy, ra, dec, q, u in zip(*event_data):
            event = event_file.bisect_met(met)
            plt.figure('IXPE single event display', figsize=(10, 10))
            grid.draw_event(event, zero_sup_threshold=threshold, padding=False, values=True)
            if abs(event.timestamp - met) <= 1.e-6:
                box = xStatBox()
                box.add_entry('MET = %.6f s' % met)
                box.add_entry('Energy = %.3f keV' % energy)
                box.add_entry('(R.A., Dec) = (%.3f, %.3f) degrees' % (ra, dec))
                box.add_entry('(Q, U) = (%.3f, %.3f)' % (q, u))
                #box.x0 = 0.025
                #box.y0 = 0.975
                #box.halign = 'left'
                #box.valign = 'top'
                box.plot(transform=plt.gcf().transFigure)
            grid.show_display()


def xpdisplay(**kwargs):
    """Application entry point.
    """
    return _run_display(kwargs.get('file'), **kwargs)


def main():
    """main() entry point.
    """
    xpdisplay(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()
