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

"""xpevtdisplay app.
"""

from __future__ import print_function, division

import os

from astropy.io import fits
import numpy

from ixpeobssim import IXPEOBSSIM_DATA
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.evt.clustering import DBscan
from ixpeobssim.evt.display import xL1EventFile, xXpolGrid, xDisplayArgumentParser,\
    load_event_list, display_event, xDisplayCard
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, xStatBox


__description__ = \
"""IXPE single event display.

This application provides visualization support for track images contained in
Level-1 IXPE files.
"""

PARSER = xDisplayArgumentParser(__description__)
PARSER.add_argument('--timestamp', type=float, default=None,
    help='timestamp of a single specific event to be displayed')



def composite_figure(obs_name=None, figsize=(12., 9.)):
    """Create a composite figure hosting all the graphical elements for the observation
    display.

    This is largely mutuated from the corresponding function in xpobsdisplay
    """
    width, height = figsize
    display_width = height / width
    text_width = 0.5 * (1. - display_width)
    # Create the top-level figure.
    fig = plt.figure('Event display (%s)' % obs_name, figsize=figsize)
    # The axes for the event display---this is guaranteed to be square, and
    # placed right on the left of the figure, spanning its entire height.
    ax_display = fig.add_axes((0., 0., display_width, 1.))
    # The axes object for the text card.
    rect = (display_width, 0., text_width, 0.75)
    ax_text = fig.add_axes(rect)
    plt.axis('off')
    return fig, ax_display, ax_text



def xpevtdisplay(**kwargs):
    """Run the event display.
    """
    # Set the random seed, if necessary.
    random_seed = kwargs.get('seed')
    if random_seed is not None:
        logger.info('Setting random seed to %d...', random_seed)
        numpy.random.seed(random_seed)

    file_path = kwargs.get('file')
    base_file_name = os.path.basename(file_path).replace('.fits', '')
    grid = xXpolGrid(cmap_name=kwargs.get('cmap'), cmap_offset=kwargs.get('cmapoffset'))
    l1_file = xL1EventFile(file_path)
    threshold = l1_file.zero_sup_threshold()
    logger.info('Zero suppression threshold: %d', threshold)
    dbscan = DBscan(threshold, min_density_points=kwargs.get('clumindensity'),
        min_cluster_size=kwargs.get('cluminsize'))

    card = xDisplayCard(kwargs.get('targetname'), l1_file.hdu_list['EVENTS'].header)
    # If we are targeting a specific event, we show it and exit immediately.
    # Note in this case we're not drawing the info box---shall we make the extra effort?
    if kwargs.get('timestamp'):
        event = l1_file.bisect_met(kwargs.get('timestamp'))
        display_event(event, grid, threshold, dbscan, base_file_name, **kwargs)
        return
    # If we are passing an event list, loop over that one. Note that in this case
    # we do have all the final, calibrated information at hand, and we can
    # display it in a dedicated box.
    # Note the autostop mechanism is implemented in the load_event_list() hook,
    # so that we get events nicely spaces across the entire observation.
    if kwargs.get('evtlist'):
        event_list = load_event_list(kwargs.get('evtlist'), **kwargs)
        for i, (met, energy, ra, dec, q, u) in enumerate(zip(*event_list)):
            event = l1_file.bisect_met(met)
            assert abs(event.timestamp - met) <= 2.e-6
            fig, ax_display, ax_text = composite_figure()
            # Draw the text card.
            plt.sca(ax_text)
            card.set_event_data(met, energy, ra, dec, q, u)
            card.draw(x0=0., y0=0.95, line_spacing=0.09)
            # Draw the actual event display.
            plt.sca(ax_display)
            display_event(event, grid, threshold, dbscan, base_file_name, **kwargs)
            plt.close(fig)
    # Finally, handle the case where we're looking at a raw level-1 file.
    # Note that we have to handle the autostop by hand in this branch of the code.
    else:
        for i, event in enumerate(l1_file):
            display_event(event, grid, threshold, dbscan, base_file_name, **kwargs)
            if i + 1 == kwargs.get('autostop'):
                break


def main():
    """main() entry point.
    """
    xpevtdisplay(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()
