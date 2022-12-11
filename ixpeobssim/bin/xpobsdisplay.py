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
    load_event_list, display_event, xDisplayCard
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.evt.kislat2015 import xStokesAnalysis
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, setup_gca_stokes,\
    plot_ellipse, DEFAULT_COLORS


__description__ = \
"""IXPE composite carousel to display observations.

This application extends xpevtdisplay to include the most relevant cumulative
distributions for a given observation.
"""

PARSER = xDisplayArgumentParser(__description__)
PARSER.add_irfname()
PARSER.add_argument('--pdmax', type=float, default=0.2,
    help='maximum polarization degree for the Stokes plot')
PARSER.add_argument('--pdstep', type=float, default=0.05,
    help='polarization degree step for the Stokes plot grid')
PARSER.add_argument('--xref', type=float, default=None,
                    help='the horizontal position of the image center')
PARSER.add_argument('--yref', type=float, default=None,
                    help='the vertical position of the image center')
PARSER.add_argument('--npix', type=int, default=200,
    help='number of pixels per side for the count map in sky coordinates')
PARSER.add_argument('--pixsize', type=float, default=None,
    help='pixel size in arcseconds for the count map in sky coordinates')
PARSER.add_argument('--subtitle', type=str, default=None,
    help='subtitle for the animation')



def composite_figure(wcs, obs_name=None, figsize=(18., 9.), left_pad=0.06,
                     right_pad=0.01, bot_pad=0.1, title_height=0.09):
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
    display_height = 1. - title_height
    plot_height = 0.5 * display_height
    # Create the top-level figure.
    fig = plt.figure('Observation display (%s)' % obs_name, figsize=figsize)
    # The axis for the title
    ax_title = fig.add_axes((0., display_height, display_width, title_height))
    plt.axis('off')
    # The axes for the event display---this is guaranteed to be square, and
    # placed right on the left of the figure, spanning its entire height.
    ax_display = fig.add_axes((0., 0., display_width, display_height))
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
    rect = (display_width, 0., plot_width * width / height, 0.5)
    ax_text = fig.add_axes(rect)
    plt.axis('off')
    ax_text.set_aspect('equal')
    return fig, ax_title, ax_display, ax_cmap, ax_cmap_colorbar, ax_polarization, ax_spectrum, ax_text


def polarization_analysis(q, u, energy, modf, aeff, mask):
    """
    """
    analysis = xStokesAnalysis(q, u, energy, modf, aeff, None)
    I, Q, U = analysis._sum_stokes_parameters(mask)
    mu = analysis._effective_mu(mask)
    W2 = analysis.W2(mask)
    QN, UN, dI, dQ, dU, dQN, dUN, cov, pval, conf, sig = analysis.calculate_stokes_errors(I, Q, U, mu, W2)
    #pd, pd_err, pa, pa_err = analysis.calculate_polarization(I, Q, U, mu, W2)
    return QN, UN, dQN, dUN, sig


def met_span(event_file):
    """Return the times for the first and the last event in a given (Level-1 or Level-2)
    file.
    """
    met = event_file.hdu_list['EVENTS'].data['TIME']
    first_met, last_met = met[0], met[-1]
    span = (last_met - first_met) / 1000.
    logger.info('Wall-clock file span: %.6f--%.6f s (%.3f ks)', first_met, last_met, span)
    return first_met, last_met


def xpobsdisplay(**kwargs):
    """Run the observation event display.
    """
    # We do need an event list, here...
    if not kwargs.get('evtlist'):
        raise RuntimeError('Please provide an event list...')

    # Set the random seed, if necessary.
    random_seed = kwargs.get('seed')
    if random_seed is not None:
        logger.info('Setting random seed to %d...', random_seed)
        numpy.random.seed(random_seed)

    # Cache the global settings...
    file_path = kwargs.get('file')
    emin, emax = kwargs.get('emin'), kwargs.get('emax')
    npix = kwargs.get('npix')
    pixsize = kwargs.get('pixsize')
    pdmax = kwargs.get('pdmax')
    pdstep = kwargs.get('pdstep')
    pd_grid = numpy.arange(pdstep, pdmax, pdstep)
    sig_color = 'black'
    sig_arrowprops=dict(arrowstyle='->', connectionstyle='angle3', color=sig_color)
    sig_kwargs = dict(xycoords='data', textcoords='axes fraction',
        arrowprops=sig_arrowprops, backgroundcolor='white', color=sig_color, ha='center')
    base_file_name = os.path.basename(file_path).replace('.fits', '')
    grid = xXpolGrid(cmap_name=kwargs.get('cmap'), cmap_offset=kwargs.get('cmapoffset'))

    # Open the Level-1 file and retrieve the necessary information.
    l1_file = xL1EventFile(file_path)
    l1_first_met, l1_last_met = met_span(l1_file)

    threshold = l1_file.zero_sup_threshold()
    logger.info('Zero suppression threshold: %d', threshold)
    # Setup the DBscan
    dbscan = DBscan(threshold, min_density_points=kwargs.get('clumindensity'),
        min_cluster_size=kwargs.get('cluminsize'))

    # Open the Level-2 file and retrieve the necessary info,
    l2_file = xEventFile(kwargs.get('evtlist'))
    l2_first_met, l2_last_met = met_span(l2_file)
    xref, yref = l2_file.wcs_reference()
    if kwargs.get('xref') is not None:
        xref = kwargs.get('xref')
    if kwargs.get('yref') is not None:
        yref = kwargs.get('yref')
    wcs_kwargs = dict(xref=xref, yref=yref, npix=npix, pixsize=pixsize)
    wcs_ = xEventBinningBase._build_image_wcs(default_img_side=10., **wcs_kwargs)
    time_data = l2_file.time_data()
    energy_data = l2_file.energy_data()
    ra_data, dec_data = l2_file.sky_position_data()
    q_data, u_data = l2_file.stokes_data()
    energy_mask = numpy.logical_and(energy_data >= emin, energy_data < emax)
    total_events = energy_mask.sum()

    # Setup all the binned data products.
    card = xDisplayCard(kwargs.get('targetname'), l2_file.hdu_list['EVENTS'].header)
    energy_binning = numpy.arange(0., 12.02, 0.04)
    hist_spec = xHistogram1d(energy_binning, xlabel='Energy [keV]', ylabel='Counts')
    cmap_data = numpy.zeros((npix, npix), dtype=float)
    aeff = load_arf(kwargs.get('irfname'), l2_file.du_id())
    modf = load_modf(kwargs.get('irfname'), l2_file.du_id())

    # Load the event data from the event list.
    event_list = load_event_list(kwargs.get('evtlist'), **kwargs)
    previous_met = time_data[0]

    # Start the loop over the event list.
    for i, (met, energy, ra, dec, q, u) in enumerate(zip(*event_list)):
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
        qn, un, dqn, dun, sig = polarization_analysis(q_data, u_data, energy_data, modf, aeff, mask)

        # Since we are at it, we take advantage of the fact that the polarization
        # analysis is re-done from the beginning every time to cache the event
        # statistics.
        num_events = mask.sum()
        elapsed_time = met - time_data[0]

        # Create the composite panel---I can't seem to be able to understand why
        # the event display is not refreshed if I don't create and delete the
        # damned thing within the event loop and destroy it at each event.
        # I am sure this is pointing at something fundamentally wrong in the
        # code and I should look at it in details...
        fig, ax_title, ax_display, ax_cmap, ax_cmap_colorbar, ax_polarization,\
            ax_spectrum, ax_text = composite_figure(wcs_)
        # Set the title.
        plt.sca(ax_title)
        title = 'Replay of a sample of events obtained by one of IXPE\'s three detectors'
        subtitle = kwargs.get('subtitle')
        plt.text(0.05, 0.7, title, size='x-large', va='center', ha='left')
        if subtitle is not None:
            plt.text(0.05, 0.3, '(%s)' % subtitle, size='large', va='center', ha='left')
        # Update the count map.
        plt.sca(ax_cmap)
        im = ax_cmap.imshow(cmap_data)
        ax_cmap.set_xlabel('Right Ascension')
        ax_cmap.set_ylabel('Declination')
        plt.grid()
        plt.colorbar(im, ax=ax_cmap_colorbar, location='top')
        # Update the polarization plot.
        plt.sca(ax_polarization)
        color = DEFAULT_COLORS[0]
        plt.plot(qn, un, 'o', color=color)
        for sigma in (1., 2., 3):
            plot_ellipse((qn, un), 2. * sigma * dqn, 2. * sigma * dun, zorder=10, color=color)
        delta = 0.5 * (dqn + dun) * sigma * numpy.sqrt(0.5)
        x0 = qn - delta * numpy.sign(qn)
        y0 = un - delta * numpy.sign(un)
        plt.text(x0, y0, '%.d$\\sigma$' % sigma, color=color, backgroundcolor='white',
            ha='center', va='center', zorder=11, clip_on=True,
            bbox=dict(boxstyle='square,pad=0.', fc='white', ec='none'))
        if sig > 3.:
            text = 'Polarization significance: %.2f $\\sigma$' % sig
            plt.gca().annotate(text, xy=(qn, un), xytext=(0.5, 1.1), **sig_kwargs)
        setup_gca_stokes(side=pdmax, pd_grid=pd_grid)
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
        # Update the cumulative statistics.
        card.update_cumulative_statistics(num_events, emin, emax)
        # Update the event data.
        card.set_event_data(met, energy, ra, dec, q, u)
        card.draw(x0=0.02, y0=0.99, line_spacing=0.08)
        # Draw the small progress bar.
        frac = num_events / total_events
        radius = 0.085
        pos = (radius + 0.025, 0.575)
        plt.gca().pie([1. - frac, frac], wedgeprops={'width': 0.025}, startangle=90,
            colors=['lightgray', color], center=pos, radius=radius)
        plt.text(*pos, '%.1f%%' % (100. * frac), size='small', ha='center',
            va='center', color='black')
        # I am not sure why, but we do have to reset the canvas size to get the display right.
        plt.gca().set_xlim(0., 1.)
        plt.gca().set_ylim(0., 1.)
        # And, finally, the actual event display---since this is blocking,
        # it needs to go last.
        plt.sca(ax_display)
        file_name = '%s_%04d.%s' % (base_file_name, i, kwargs.get('imgformat'))
        display_event(event, grid, threshold, dbscan, file_name, **kwargs)
        # See the initial remark about the need to destroy the figure.
        plt.close(fig)


def main():
    """main() entry point.
    """
    xpobsdisplay(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()
