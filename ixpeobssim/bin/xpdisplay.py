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

import os

from astropy.io import fits
import numpy

from ixpeobssim import IXPEOBSSIM_DATA
from ixpeobssim.core.hist import xHistogram1d
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
PARSER.add_argument('--timestamp', type=float, default=None,
    help='timestamp of a single specific event to be displayed')
PARSER.add_argument('--evtlist', type=str,
    help='path to the auxiliary (Level-2 file) event list')
PARSER.add_argument('--resample', type=float, default=None,
    help='the power-law index for resampling events in energy')
PARSER.add_boolean('--absorption', True,
    help='draw the reconstructed absorption_point')
PARSER.add_boolean('--barycenter', True,
    help='draw the reconstructed barycenter')
PARSER.add_boolean('--direction', True,
    help='draw the reconstructed track direction')
PARSER.add_boolean('--pixpha', False,
    help='indicate the pixel PHA values')
PARSER.add_boolean('--indices', False,
    help='draw the row and column indices of the readout matrix')
PARSER.add_argument('--cmap', type=str, default='Reds',
    help='the color map for the pixel values')
PARSER.add_argument('--cmapoffset', type=int, default=10,
    help='the PHA offset for the color map')
PARSER.add_argument('--axside', type=float, default=None,
    help='the axis side for the event display')
PARSER.add_boolean('--autosave', False,
    help='save the event displays automatically')
PARSER.add_argument('--imgformat', type=str, default='png',
    help='the image format for the output files when autosave is True')




def event_box(met, energy, ra, dec, q, u):
    """Draw a text box with the event information.
    """
    box = xStatBox()
    box.add_entry('MET = %.6f s' % met)
    box.add_entry('Energy = %.3f keV' % energy)
    box.add_entry('(R.A., Dec) = (%.3f, %.3f) degrees' % (ra, dec))
    box.add_entry('(Q, U) = (%.3f, %.3f)' % (q, u))
    box.x0 = 0.025
    box.y0 = 0.975
    box.plot(transform=plt.gcf().transFigure)


def load_level_2_data(file_path, resample_index=None, pivot_energy=8., interactive=False):
    """Load the event data from the Level-2 event list.
    """
    event_list = xEventFile(file_path)
    met = event_list.time_data()
    energy = event_list.energy_data()
    ra, dec = event_list.sky_position_data()
    q, u = event_list.stokes_data()
    if resample_index is not None:
        logger.info('Resampling input level-2 data with index %.3f', resample_index)
        mask = numpy.random.uniform(size=len(energy)) <= (energy /  pivot_energy)**resample_index
        logger.info('%d event(s) out of %s remaining.', mask.sum(), len(mask))
        met, energy, ra, dec, q, u = [item[mask] for item in (met, energy, ra, dec, q, u)]
    if interactive:
        # Debug plot for the input energy spectrum.
        plt.figure('Input energy spectrum')
        h = xHistogram1d(numpy.linspace(2., 8., 20)).fill(energy)
        h.plot()
    return met, energy, ra, dec, q, u


def display_event(event, grid, threshold, base_file_name=None, box_info=None,
    padding=False, **kwargs):
    """Single-stop event display.
    """
    draw_kwargs = dict(values=kwargs.get('pixpha'), indices=kwargs.get('indices'),
        canvas_side=kwargs.get('axside'), zero_sup_threshold=threshold, padding=padding)
    plt.figure('IXPE single event display', figsize=(9., 10.))
    logger.info('Drawing event @ MET %.6f', event.timestamp)
    # Draw the bare event...
    grid.draw_event(event, **draw_kwargs)
    # ... then the reconstruction elements.
    if kwargs.get('absorption'):
        event.recon.draw_absorption_point()
    if kwargs.get('barycenter'):
        event.recon.draw_barycenter()
    if kwargs.get('direction'):
        event.recon.draw_track_direction()
    if box_info is not None:
        event_box(*box_info)
    if kwargs.get('autosave'):
        file_name = '%s_%.6f.%s' % (base_file_name, event.timestamp, kwargs.get('imgformat'))
        file_path = os.path.join(IXPEOBSSIM_DATA, file_name)
    else:
        file_path = None
    grid.show_display(file_path)


def run_display(file_path, **kwargs):
    """Run the event display.
    """
    base_file_name = os.path.basename(file_path).replace('.fits', '')
    grid = xXpolGrid(cmap_name=kwargs.get('cmap'), cmap_offset=kwargs.get('cmapoffset'))
    event_file = xL1EventFile(file_path)
    threshold = event_file.zero_sup_threshold()
    logger.info('Zero suppression threshold: %d', threshold)
    #draw_kwargs = dict(values=kwargs.get('pixpha'), indices=kwargs.get('indices'),
    #    zero_sup_threshold=threshold, padding=False, canvas_side=kwargs.get('axside'))
    # If we are targeting a specific event, we show it and exit immediately.
    # Note in this case we're not drawing the info box---shall we make the extra effort?
    if kwargs.get('timestamp'):
        event = event_file.bisect_met(kwargs.get('timestamp'))
        display_event(event, grid, threshold, base_file_name, **kwargs)
        return
    # If we are passing an event list, loop over that one. Note that in this case
    # we do have all the final, calibrated information at hand, and we can
    # display it in a dedicated box.
    if kwargs.get('evtlist'):
        l2_data = load_level_2_data(kwargs.get('evtlist'), kwargs.get('resample'))
        for met, energy, ra, dec, q, u in zip(*l2_data):
            event = event_file.bisect_met(met)
            # We only show the calibrated information if the timestamps in the
            # level 1 and level 2 files agree to within 2 mus.
            if abs(event.timestamp - met) <= 2.e-6:
                box_info = (met, energy, ra, dec, q, u)
            else:
                box_info = None
            display_event(event, grid, threshold, base_file_name, box_info, **kwargs)
    # Finally, handle the case where we're looking at a raw level-1 file.
    else:
        for event in event_file:
            display_event(event, grid, threshold, base_file_name, **kwargs)


def xpdisplay(**kwargs):
    """Application entry point.
    """
    return run_display(kwargs.get('file'), **kwargs)


def main():
    """main() entry point.
    """
    xpdisplay(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()
