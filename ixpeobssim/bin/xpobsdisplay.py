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
from ixpeobssim.irf import load_arf, load_modf
from ixpeobssim.evt.clustering import DBscan
from ixpeobssim.evt.display import xL1EventFile, xXpolGrid, xDisplayArgumentParser,\
    event_box, load_level_2_data, display_event
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.evt.kislat2015 import xStokesAnalysis
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, setup_gca_stokes,\
    xTextCard, plot_ellipse, DEFAULT_COLORS


__description__ = \
"""IXPE single event display.

This application provides visualization support for track images contained in
Level-1 IXPE files.
"""

PARSER = xDisplayArgumentParser(__description__)
PARSER.add_irfname()



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
    #setup_gca_stokes()
    # The spectrum axes on the bottom-right corner.
    rect = (display_width + plot_width + 2. * left_pad, bot_pad, plot_width, 0.5 - bot_pad)
    ax_spectrum = fig.add_axes(rect)
    rect = (display_width + left_pad, 0., plot_width, 0.5)
    ax_text = fig.add_axes(rect)
    plt.axis('off')
    return fig, ax_display, ax_skymap, ax_polarization, ax_spectrum, ax_text



def xpobsdisplay(**kwargs):
    """Run the observation event display.
    """
    file_path = kwargs.get('file')
    base_file_name = os.path.basename(file_path).replace('.fits', '')
    grid = xXpolGrid(cmap_name=kwargs.get('cmap'), cmap_offset=kwargs.get('cmapoffset'))
    # Open the Level-1 file and retrieve the necessary information.
    l1_file = xL1EventFile(file_path)
    threshold = l1_file.zero_sup_threshold()
    logger.info('Zero suppression threshold: %d', threshold)
    # Setup the DBscan
    dbscan = DBscan(threshold, min_density_points=kwargs.get('clumindensity'),
        min_cluster_size=kwargs.get('cluminsize'))
    # Open the Level-2 file and retrieve the necessary info,
    l2_file = xEventFile(kwargs.get('evtlist'))
    xref, yref = l2_file.wcs_reference()
    wcs_kwargs = dict(xref=xref, yref=yref, npix=200)
    wcs_ = xEventBinningBase._build_image_wcs(default_img_side=10., **wcs_kwargs)
    energy_data = l2_file.energy_data()
    ra_data, dec_data = l2_file.sky_position_data()
    q_data, u_data = l2_file.stokes_data()
    # Prepate the text card with the observation-related info that does not
    # change within the event loop.
    header = l2_file.hdu_list['EVENTS'].header
    card = xTextCard()
    card.set_line('Target Name', header['OBJECT'])
    card.set_line('Observation Start', header['DATE-OBS'])
    card.set_line('Observation End', header['DATE-END'])
    card.set_line('Detector Unit', '%s (%s)' % (header['DETNAM'], header['DET_ID']))
    card.set_line('Spacer', None)

    # Temporary stuff: for the time being we fill the plots all at once at the
    # beginning, but this should happen in the event loop.
    # Note the number of bins is relevant, here. 2--8 keV in steps of 40 eV is 150 bins.
    hist_spec = xHistogram1d(numpy.linspace(2., 8., 151), xlabel='Energy [keV]', ylabel='Counts')
    hist_spec.fill(energy_data)
    x, y, binning = xEventBinningBase._pixelize_skycoords(ra_data, dec_data, wcs_)
    hist_skymap, _, _ = numpy.histogram2d(x, y, bins=binning)
    aeff = load_arf(kwargs.get('irfname'), l2_file.du_id())
    modf = load_modf(kwargs.get('irfname'), l2_file.du_id())
    analysis = xStokesAnalysis(q_data, u_data, energy_data, modf, aeff, l2_file.livetime())
    # ... end of hack

    # We do need an event list, here...
    if not kwargs.get('evtlist'):
        raise RuntimeError('Please provide an event list...')
    l2_data = load_level_2_data(kwargs.get('evtlist'), **kwargs)
    for t, E, ra, dec, q, u in zip(*l2_data):
        event = l1_file.bisect_met(t)
        assert abs(event.timestamp - t) <= 2.e-6
        # Create the composite panel---I can't seem to be able to understand why
        # the event display is not refreshed if I don't create and delete the
        # damned thing within the event loop and destroy it at each event.
        # I am sure this is pointing at something fundamentally wrong in the
        # code and I should look at it in details...
        fig, ax_display, ax_skymap, ax_polarization, ax_spectrum, ax_text = composite_figure(wcs_)
        # Update the skymap.
        ax_skymap.imshow(hist_skymap)
        ax_skymap.set_xlabel('Right Ascension')
        ax_skymap.set_ylabel('Declination')
        # Update the polarization plot.
        plt.sca(ax_polarization)
        # Need to refine this!
        table = analysis.polarization_table(numpy.linspace(2., 8., 2))
        for sigma in (1., 2., 3):
            q, u, dq, du = [table[key] for key in ('QN', 'UN', 'QN_ERR', 'UN_ERR')]
            plot_ellipse((q, u), 2. * sigma * dq, 2. * sigma * du, zorder=10,
                color=DEFAULT_COLORS[0])
        setup_gca_stokes(side=0.12, pd_grid=numpy.linspace(0.05, 0.1, 2))
        # Update the spectrum.
        plt.sca(ax_spectrum)
        hist_spec.plot()
        setup_gca(logy=True, grids=True)
        # Update the text card.
        plt.sca(ax_text)
        card.set_line('Mission elsapsed time', t, '%.6f', 's')
        card.set_line('Energy', E, '%.2f', 'keV')
        card.set_line('Right ascention', ra, '%.3f', 'decimal degrees')
        card.set_line('Declination', dec, '%.3f', 'decimal degrees')
        card.set_line('Stokes parameters', '(%.4f, %.4f)' % (q, u))
        card.draw(x0=0., y0=0.95, line_spacing=0.09)
        # And, finally, the actual event display---since this is blocking,
        # it needs to go last.
        plt.sca(ax_display)
        display_event(event, grid, threshold, dbscan, base_file_name, **kwargs)
        # See the initial remark about the need to destroy the figure.
        plt.close(fig)


def main():
    """main() entry point.
    """
    xpobsdisplay(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()
