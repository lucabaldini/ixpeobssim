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

"""xpobsdisplay app.
"""

from __future__ import print_function, division

import os

from astropy.io import fits
from astropy.visualization.wcsaxes import WCSAxes
import matplotlib
import numpy

from ixpeobssim import IXPEOBSSIM_DATA
from ixpeobssim.binning.base import xEventBinningBase
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.evt.clustering import DBscan
from ixpeobssim.evt.display import xL1EventFile, xXpolGrid, xDisplayArgumentParser,\
    event_box, load_level_2_data, display_event
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, setup_gca_stokes


__description__ = \
"""IXPE single event display.

This application provides visualization support for track images contained in
Level-1 IXPE files.
"""

PARSER = xDisplayArgumentParser(__description__)



def composite_figure(wcs, obs_name=None, figsize=(18., 9.), left_pad=0.06,
                     right_pad=0.01, bot_pad=0.1):
    """Create a composite figure hosting all the graphical elements for the observation
    display.

    I have been struggling with this for many hours, trying with subplots and
    gridspecs, and finally came to the realization that the only sensible way to
    achieve a reasonable positioning for the various pads was to use
    matplotlib.figure.add_axes() and calculate all the coordinates by hand.
    """
    width, height = figsize
    aspect_ratio = width / height
    display_width = height / width
    plot_width = 0.5 * (1. - display_width) - left_pad - right_pad
    # Create the top-level figure.
    fig = plt.figure('Observation display (%s)' % obs_name, figsize=figsize)
    # The axes for the event display---this is guaranteed to be square, and
    # placed right on the left of the figure, spanning its entire height.
    ax_display = fig.add_axes((0., 0., display_width, 1.))
    # Now the skymap and the polarization axes on the top left of the figure.
    # Note these have a square aspect ratio, and will naturally place themselves
    # in the right place vertically.
    rect = (display_width + left_pad, 0.5, plot_width, 0.5)
    ax_skymap = fig.add_axes(rect, projection=wcs)
    ax_skymap.set_aspect('equal')
    rect = (display_width + plot_width + 2. * left_pad, 0.5, plot_width, 0.5)
    ax_polarization = fig.add_axes(rect)
    setup_gca_stokes()
    # The spectrum axes on the bottom-right corner.
    rect = (display_width + plot_width + 2. * left_pad, bot_pad, plot_width, 0.5 - bot_pad)
    ax_spectrum = fig.add_axes(rect)
    rect = (display_width + left_pad, 0., plot_width, 0.5)
    ax_text = fig.add_axes(rect)
    plt.axis('off')
    return fig, ax_display, ax_skymap, ax_spectrum, ax_polarization


def xpobsdisplay(**kwargs):
    """Run the observation event display.
    """
    file_path = kwargs.get('file')
    base_file_name = os.path.basename(file_path).replace('.fits', '')
    grid = xXpolGrid(cmap_name=kwargs.get('cmap'), cmap_offset=kwargs.get('cmapoffset'))
    l1_file = xL1EventFile(file_path)
    threshold = l1_file.zero_sup_threshold()
    logger.info('Zero suppression threshold: %d', threshold)
    dbscan = DBscan(threshold, min_density_points=kwargs.get('clumindensity'),
        min_cluster_size=kwargs.get('cluminsize'))

    l2_file = xEventFile(kwargs.get('evtlist'))
    energy = l2_file.energy_data()

    # Build the glorious, composite panel.
    wcs = xEventBinningBase._build_image_wcs()
    # Temporary stuff...
    # Note the number of bins is relevant, here. 2--8 keV in steps of 40 eV is 150 bins.
    hist_spec = xHistogram1d(numpy.linspace(2., 8., 151), xlabel='Energy [keV]',
        ylabel='Counts')
    hist_spec.fill(energy)
    # ... end of hack

    # We do need an event list, here...
    if not kwargs.get('evtlist'):
        raise RuntimeError('Please provide an event list...')
    l2_data = load_level_2_data(kwargs.get('evtlist'), **kwargs)
    for met, energy, ra, dec, q, u in zip(*l2_data):
        event = l1_file.bisect_met(met)
        # We only show the calibrated information if the timestamps in the
        # level 1 and level 2 files agree to within 2 mus.
        #if abs(event.timestamp - met) <= 2.e-6:
        #    box_info = (met, energy, ra, dec, q, u)
        #else:
        #box_info = None

        fig, ax_display, ax_skymap, ax_spectrum, ax_polarization = composite_figure(wcs)

        plt.sca(ax_spectrum)
        hist_spec.plot()
        setup_gca(logy=True, grids=True)

        # Event display is blocking, and goes last.
        plt.sca(ax_display)
        display_event(event, grid, threshold, dbscan, base_file_name, **kwargs)



        plt.close(fig)



def main():
    """main() entry point.
    """
    xpobsdisplay(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()
    # ra, dec = 0., 0.
    # x, y, binning = xEventBinningBase._pixelize_skycoords(ra, dec, wcs_)
    # data, _, _ = numpy.histogram2d(x, y, bins=binning)
    # ax_hists[0].imshow(data)
    # ax_hists[0].set_xlabel('Right Ascension')
    # ax_hists[0].set_ylabel('Declination')
