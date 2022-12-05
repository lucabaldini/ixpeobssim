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
    event_box, load_event_list, display_event, xDisplayCard
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.evt.kislat2015 import xStokesAnalysis
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, setup_gca_stokes,\
    plot_ellipse, DEFAULT_COLORS


__description__ = \
"""IXPE composite carousel to display observations.
"""

PARSER = xDisplayArgumentParser(__description__)
PARSER.add_irfname()
PARSER.add_ebounds()
PARSER.add_argument('--npix', type=int, default=200,
    help='number of pixels per side for the count map in sky coordinates')



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
    # Now the count map and the polarization axes on the top left of the figure.
    # Note these have a square aspect ratio, and will naturally place themselves
    # in the right place vertically.
    rect = (display_width + left_pad, 0.5, plot_width, 0.5)
    ax_cmap = fig.add_axes(rect, projection=wcs)
    ax_cmap.set_aspect('equal')
    # This additional, ghost axes object is to hold the color bar count map,
    # so that we get it displayed without the count map axes being resized---this
    # a horrible hack, but I could not figure out a better way to do it.
    rect = (display_width + left_pad, 0.525, plot_width, 0.5)
    ax_cmap_colorbar = fig.add_axes(rect)
    plt.axis('off')
    rect = (display_width + plot_width + 2. * left_pad, 0.5, plot_width, 0.5)
    ax_polarization = fig.add_axes(rect)
    #setup_gca_stokes()
    # The spectrum axes on the bottom-right corner.
    rect = (display_width + plot_width + 2. * left_pad, bot_pad, plot_width, 0.475 - bot_pad)
    ax_spectrum = fig.add_axes(rect)
    rect = (display_width + left_pad, 0., plot_width, 0.5)
    ax_text = fig.add_axes(rect)
    plt.axis('off')
    return fig, ax_display, ax_cmap, ax_cmap_colorbar, ax_polarization, ax_spectrum, ax_text


def polarization_analysis(q, u, energy, modf, aeff, mask):
    """
    """
    analysis = xStokesAnalysis(q, u, energy, modf, aeff, None)
    I, Q, U = analysis._sum_stokes_parameters(mask)
    mu = analysis._effective_mu(mask)
    W2 = analysis.W2(mask)
    QN, UN, dI, dQ, dU, dQN, dUN, cov, pval, conf, sig = analysis.calculate_stokes_errors(I, Q, U, mu, W2)
    #pd, pd_err, pa, pa_err = analysis.calculate_polarization(I, Q, U, mu, W2)
    return QN, UN, dQN, dUN


def xpobsdisplay(**kwargs):
    """Run the observation event display.

    TODO:

    * add significance to Stokes plot
    * add colorbar to the display
    * add sigma levels on Stokes plot
    """
    # We do need an event list, here...
    if not kwargs.get('evtlist'):
        raise RuntimeError('Please provide an event list...')

    # Cache the global settings...
    file_path = kwargs.get('file')
    emin, emax = kwargs.get('emin'), kwargs.get('emax')
    npix = kwargs.get('npix')
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
    wcs_kwargs = dict(xref=xref, yref=yref, npix=npix)
    wcs_ = xEventBinningBase._build_image_wcs(default_img_side=10., **wcs_kwargs)
    time_data = l2_file.time_data()
    energy_data = l2_file.energy_data()
    ra_data, dec_data = l2_file.sky_position_data()
    q_data, u_data = l2_file.stokes_data()
    energy_mask = numpy.logical_and(energy_data >= emin, energy_data < emax)

    # Setup all the binned data products.
    card = xDisplayCard(l2_file.hdu_list['EVENTS'].header)
    energy_binning = numpy.arange(0., 12.02, 0.04)
    hist_spec = xHistogram1d(energy_binning, xlabel='Energy [keV]', ylabel='Counts')
    cmap_data = numpy.zeros((npix, npix), dtype=float)
    aeff = load_arf(kwargs.get('irfname'), l2_file.du_id())
    modf = load_modf(kwargs.get('irfname'), l2_file.du_id())

    # Load the event data from the event list.
    event_list = load_event_list(kwargs.get('evtlist'), **kwargs)
    previous_met = time_data[0]

    # Start the loop over the event list.
    for met, energy, ra, dec, q, u in zip(*event_list):
        # Retrieve the actual event from the underlying level-1 file.
        event = l1_file.bisect_met(met)
        assert abs(event.timestamp - met) <= 2.e-6
        # Create the mask for the current chunck of data.
        logger.info('Filtering events between MET %.6f and %.6f...', previous_met, met)
        mask = numpy.logical_and(time_data >= previous_met, time_data < met)
        logger.info('Done, %d event(s) left before the energy cut.', mask.sum())

        # Update all the binned products: the energy spcetrum...
        hist_spec.fill(energy_data[mask])
        # ... the count map...
        mask *= energy_mask
        x, y, binning = xEventBinningBase._pixelize_skycoords(ra_data[mask], dec_data[mask], wcs_)
        counts, _, _ = numpy.histogram2d(x, y, bins=binning)
        cmap_data += counts
        # ... and the polarization analysis. Note that, instead of keeping track
        # of all the necessary stuff to accumulate the polarization analysis in
        # cuncks, we repeat the entire calculation feeding in all the events up to
        # the current met---this is slightly suboptimal, but simpler and less
        # error-prone.
        mask = time_data < met
        mask *= energy_mask
        qn, un, dqn, dun = polarization_analysis(q_data, u_data, energy_data, modf, aeff, mask)

        # Create the composite panel---I can't seem to be able to understand why
        # the event display is not refreshed if I don't create and delete the
        # damned thing within the event loop and destroy it at each event.
        # I am sure this is pointing at something fundamentally wrong in the
        # code and I should look at it in details...
        fig, ax_display, ax_cmap, ax_cmap_colorbar, ax_polarization,\
            ax_spectrum, ax_text = composite_figure(wcs_)
        # Update the count map.
        plt.sca(ax_cmap)
        im = ax_cmap.imshow(cmap_data)
        ax_cmap.set_xlabel('Right Ascension')
        ax_cmap.set_ylabel('Declination')
        plt.grid()
        plt.colorbar(im, ax=ax_cmap_colorbar, location='top')
        # Update the polarization plot.
        plt.sca(ax_polarization)
        for sigma in (1., 2., 3):
            plot_ellipse((qn, un), 2. * sigma * dqn, 2. * sigma * dun, zorder=10,
                color=DEFAULT_COLORS[0])
        setup_gca_stokes(side=0.12, pd_grid=numpy.linspace(0.05, 0.1, 2))
        # Update the spectrum.
        plt.sca(ax_spectrum)
        hist_spec.plot()
        setup_gca(logy=True, grids=True, xticks=numpy.arange(0., 12.1, 2.))
        for x in emin, emax:
            plt.axvline(x, color='gray')
        plt.axvspan(0., emin, alpha=0.25, color='gray')
        plt.axvspan(emax, 12., alpha=0.25, color='gray')
        # Update the text card.
        plt.sca(ax_text)
        card.set_event_data(met, energy, ra, dec, q, u)
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