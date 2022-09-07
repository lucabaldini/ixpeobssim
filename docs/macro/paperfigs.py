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

from __future__ import print_function, division

import matplotlib
import numpy
from astropy.visualization.wcsaxes import WCSAxes

from ixpeobssim.binning.misc import xBinnedMap
from ixpeobssim.binning.polarization import xBinnedPolarizationMapCube
from ixpeobssim.core.fitsio import xFITSImageBase
from ixpeobssim.core.hist import xHistogram2d
import ixpeobssim.core.pipeline as pipeline
import ixpeobssim.evt.xspec_ as xspec_
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.instrument.sc import dithering_pattern
from ixpeobssim.irf import load_arf, load_psf
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, _set_rc_param
from ixpeobssim.utils.units_ import degrees_to_arcmin

from ixpeobssim.config.crab_pulsar import ROI_MODEL as CRAB_PRS_ROI


# Make the font slightly bigger.
_set_rc_param('font.size', 14.0)



def plot_count_spectrum(energy=6., phase=0.6):
    """
    """
    psr = CRAB_PRS_ROI['Crab pulsar']
    aeff = load_arf()
    time_grid = numpy.linspace(0., 1., 100)
    count_spectrum = psr.create_count_spectrum(aeff, time_grid)
    count_spectrum.zlabel = 'Count spectrum [s$^{-1}$ keV$^{-1}$]'
    plt.figure('Crab count spectrum', figsize=(8, 7.2))
    plt.subplots_adjust(top=0.80, bottom=0.1)
    count_spectrum.plot(logz=True, vmin=1.1e-6, num_contours=100)
    fmt = matplotlib.ticker.LogFormatterMathtext()
    fmt.create_dummy_axis()
    count_spectrum.plot_contours(logz=True, cfmt=fmt)
    setup_gca(ylabel='Pulse phase')

    y = 0.78
    dx = 0.42
    dy = 0.40

    axins1 = plt.gca().inset_axes([-0.08, y, dx, dy])
    slice1 = count_spectrum.hslice(phase)
    axins1.plot(slice1.x, slice1.y, color='black')
    axins1.set_yscale('log')
    axins1.set_xlim(1., 12.)
    axins1.set_xlabel('Energy [keV]')
    axins1.set_xticks(numpy.linspace(2., 10., 5))
    axins1.xaxis.tick_top()
    axins1.xaxis.set_label_position('top')
    axins1.grid(which='major')
    plt.axhline(phase, color='white', ls='dashed')

    plt.gca().annotate('', xy=(2., phase), xycoords='data', xytext=(0, 65),
        textcoords='offset points', size=20,
        arrowprops=dict(arrowstyle="wedge,tail_width=0.7", fc="white", ec="white"))

    axins2 = plt.gca().inset_axes([0.80, y, dx, dy])
    slice2 = count_spectrum.vslice(energy)
    axins2.plot(slice2.x, slice2.y, color='black')
    axins2.set_xlim(0., 1.)
    axins2.set_xlabel('Pulse phase')
    axins2.xaxis.tick_top()
    axins2.xaxis.set_label_position('top')
    axins2.yaxis.tick_right()
    axins2.yaxis.set_label_position('right')
    axins2.grid(which='major')
    plt.axvline(energy, color='white', ls='dashed')

    plt.gca().annotate('', xy=(energy, 0.9), xycoords='data', xytext=(130, 0),
        textcoords='offset points', size=20,
        arrowprops=dict(arrowstyle="wedge,tail_width=0.7", fc="white", ec="white"))


    plt.figure('Crab ppf', figsize=(8, 7.2))
    logz = True
    ppf = count_spectrum.ppf
    ppf.zlabel = 'Percent point function [keV]'
    plt.subplots_adjust(top=0.80, bottom=0.1)
    ppf.plot(logz=True, vmin=1., vmax=12., num_contours=100)
    ppf.plot_contours(num_contours=[2, 3, 4, 5])
    setup_gca(ylabel='Pulse phase', xlabel='$\\xi$')

    axins = plt.gca().inset_axes([-0.08, y, dx, dy])
    slice = ppf.hslice(phase)
    axins.plot(slice.x, slice.y, color='black')
    axins.set_ylim(1., 12.)
    axins.set_yscale('log')
    axins.set_xlabel('$\\xi$')
    axins.set_xticks(numpy.linspace(0., 1., 6))
    yticks = [1, 2, 3, 4, 5, 10]
    axins.set_yticks(yticks)
    axins.set_yticklabels(yticks)
    axins.xaxis.tick_top()
    axins.xaxis.set_label_position('top')
    axins.grid(which='major')
    plt.axhline(phase, color='white', ls='dashed')

    plt.gca().annotate('', xy=(0.2, phase), xycoords='data', xytext=(0, 65),
        textcoords='offset points', size=20,
        arrowprops=dict(arrowstyle="wedge,tail_width=0.7", fc="white", ec="white"))


def plot_dithering_path(start_met, duration, step=1, ylabel='y [arcmin]'):
    """
    """
    dithering = dithering_pattern()
    t = numpy.arange(start_met, duration + 0.5 * step, step)
    x, y = dithering(t)
    plt.plot(x, y, color='black')
    plt.gca().set_aspect('equal')
    setup_gca(xmin=-2., xmax=2., ymin=-2., ymax=2., grids=True,
        xlabel='x [arcmin]', ylabel=ylabel)
    plt.text(1.8, 1.6, '%d s' % t[-1], fontsize=12, ha='right')


def plot_dithering():
    """
    """
    t = numpy.linspace(0., 50000., 10000000)
    x, y = dithering_pattern()(t)
    psf = load_psf(du_id=1)
    dx, dy = psf.delta(len(t))
    x += degrees_to_arcmin(dx)
    y += degrees_to_arcmin(dy)
    binning = numpy.linspace(-2., 2., 200)
    h = xHistogram2d(binning, binning, zlabel='Scaled counts [a. u.]').fill(x, y)
    h.content /= h.content.max()

    fig = plt.figure('Dithering', (8, 10))
    ax = fig.add_gridspec(3, 3, hspace=0.35, bottom=0.075, top=0.975)
    ax1 = fig.add_subplot(ax[1:3, 0:3])
    h.plot()
    ax1.set_aspect('equal')
    setup_gca(xmin=-2., xmax=2., ymin=-2., ymax=2., grids=True,
        xlabel='x [arcmin]', ylabel='y [arcmin]')

    ax0 = fig.add_subplot(ax[0, 0])
    plot_dithering_path(0, 500)
    ax1 = fig.add_subplot(ax[0, 1])
    plot_dithering_path(0, 5000, ylabel=None)
    ax2 = fig.add_subplot(ax[0, 2])
    plot_dithering_path(0, 50000, ylabel=None)


def plot_xspec(duration=100000., rebin=2):
    """
    """
    pipeline.reset('toy_point_source')
    pipeline.xpobssim(duration=duration, saa=False, occult=False)
    for algorithm in ['PHA1', 'PHA1Q', 'PHA1U']:
        file_list = pipeline.xpbin(*pipeline.file_list(), algorithm=algorithm)
        pipeline.xpgrppha(*file_list, comm='GROUP 0 275 %d' % rebin)
    file_list = pipeline.file_list('pha1*', 'grppha')
    fit_output = pipeline.xpxspec(*file_list, model='pollin * powerlaw', plot=False)
    #xspec_.plot()
    fig, axs = plt.subplots(6, 1, figsize=(8, 10), sharex=True,
        gridspec_kw=dict(bottom=0.06, top=0.98, height_ratios=[1., 0.4, 1., 0.4, 1., 0.4]))
    # Remove horizontal space between axes
    fig.subplots_adjust(hspace=0.15)
    label_dict = {'PHA1': 'I', 'PHA1Q': 'Q', 'PHA1U': 'U'}
    for i, alg in enumerate(['PHA1', 'PHA1Q', 'PHA1U']):
        for du_id in DU_IDS:
            data = xspec_.retrieve_plot_data(alg, du_id)
            plt.sca(axs[i * 2])
            xspec_.plot_normalized_counts(data, du_id)
            label = '%s [s$^{-1}$ keV$^{-1}$]' % label_dict[alg]
            setup_gca(grids=True, logy=True, ylabel=label)
            if i == 0:
                fit_data = xspec_.current_fit_output()
                #fit_data.stat_box(position='lower left').plot()
                plt.legend()
            plt.sca(axs[i * 2 + 1])
            xspec_.plot_residuals(data)
            if i != 2:
                plt.gca().set_xlabel('')
    fig.align_ylabels()


def plot_polarization_maps(duration=2000000.):
    """
    """
    pipeline.reset('casa')
    pipeline.xpobssim(duration=duration, saa=False, occult=False)
    pipeline.xpbin(*pipeline.file_list(), algorithm='CMAP', npix=250, pixsize=2)
    pipeline.xpbin(*pipeline.file_list(), algorithm='PMAP', npix=50, pixsize=10)
    cmap = xBinnedMap.from_file_list(pipeline.file_list('cmap'))
    pmap = xBinnedPolarizationMapCube.from_file_list(pipeline.file_list('pmap'))

    fig = plt.figure('Polarization map', (8, 10))
    ax = fig.add_gridspec(2, 2, bottom=0.075, top=0.975, height_ratios=(1, 2.45), hspace=0.025)
    ax1 = fig.add_subplot(ax[1, 0:2], axes_class=WCSAxes)
    ax1.reset_wcs(cmap.fits_image.wcs)
    ax1.imshow(cmap.fits_image.data)
    ax1.grid(color='gray')
    ax1.set_xlabel('Right Ascension (J2000)')
    ax1.set_ylabel('Declination (J2000)')
    mask = pmap.calculate_significance_mask(3., 0.)
    pmap._overlay_arrows(0, mask)
    plt.text(0.05, 0.92, 'Stokes I', color='white', transform = ax1.transAxes)

    ax2 = fig.add_subplot(ax[0, 0], axes_class=WCSAxes)
    ax2.reset_wcs(pmap.wcs, slices=('x', 'y', 0))
    ax2.imshow(pmap.Q[0])
    for axis in ('x', 'y'):
        ax2.tick_params(axis, labelsize=0.)
    ax2.set_xlabel(' ')
    ax2.set_ylabel(' ')
    ax2.grid(color='gray')
    plt.text(0.1, 0.85, 'Stokes Q', color='white', transform = ax2.transAxes)

    ax3 = fig.add_subplot(ax[0, 1], axes_class=WCSAxes)
    ax3.reset_wcs(pmap.wcs, slices=('x', 'y', 0))
    ax3.imshow(pmap.U[0])
    for axis in ('x', 'y'):
        ax3.tick_params(axis, labelsize=0.)
    ax3.set_xlabel(' ')
    ax3.set_ylabel(' ')
    ax3.grid(color='gray')
    plt.text(0.1, 0.85, 'Stokes U', color='white', transform = ax3.transAxes)



if __name__ == '__main__':
    plot_count_spectrum()
    #plot_dithering()
    #plot_xspec()
    #plot_polarization_maps()
    plt.show()
