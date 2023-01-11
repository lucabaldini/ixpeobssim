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

"""xpradialprofile app.
"""

from __future__ import print_function, division

import os

import numpy
from scipy.optimize import curve_fit

from ixpeobssim.core.fitting import fit_histogram
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.core.modeling import xFitModelBase, xConstant, xLorentzian
from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.irf import load_psf
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.astro import angular_separation
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca, last_line_color, residual_plot
from ixpeobssim.utils.os_ import check_input_file
from ixpeobssim.utils.units_ import degrees_to_arcmin, arcmin_to_degrees, arcmin_to_arcsec


__description__ = \
"""Create a radial profile plot, in sky coordinates, for a given observation.

This program takes a level-2 file as an input and creates a radial histogram of
the event counts, weighted with 1. / r to correct for the solid angle. This is
a useful diagnostics to gauge the background level in a given observation.

In order to provide more context to the plots, the PSF radial profile for the
relevant DU is overlayied, and the estimated fraction of source events outside
a circle of radius r is annotated over the PSF line for a discrete set of values
of r.
"""


PARSER = xArgumentParser(description=__description__)
PARSER.add_filelist()
PARSER.add_argument('--ra', type=float, default=None,
                    help='RA of the center of the radial plot in decimal degrees')
PARSER.add_argument('--dec', type=float, default=None,
                    help='Dec of the center of the radial plot in decimal degrees')
PARSER.add_argument('--rmin', type=float, default=0.,
                    help='minimum radial distance in arcminutes')
PARSER.add_argument('--rmax', type=float, default=10.,
                    help='maximum radial distance in arcminutes')
PARSER.add_argument('--rcore', type=float, default=0.5,
                    help='core radius for the psf fit in arcminutes')
PARSER.add_argument('--rbins', type=int, default=200,
                    help='number of bins for the radial profile')
PARSER.add_boolean('--autocenter', default=True,
                   help='recenter the image based on the count map')
PARSER.add_boolean('--psf', default=True,
                   help='overimpose the PSF radial profile')
PARSER.add_boolean('--residuals', default=False,
                   help='plot the estimated residual background')
PARSER.add_boolean('--interactive', default=False,
                   help='plot some diagnostic plot')
PARSER.add_irfname()



def fit_offset(delta_ra, delta_dec, rmax=2., num_bins=100, interactive=False):
    """Fit the angular separation between the sky positions of the events and the
    nominal source position in the two ra, dec projections with a Lorentzian model,
    and return the optimal offset values to recenter the count map.

    Arguments
    ---------
    delta_ra : array_like
        The R. A. difference between the sky position and the nominal source position in
        decimal degrees.

    delta_dec : array_like
        The Dec. difference between the sky position and the nominal source position in
        decimal degrees.

    r_max : float
        The maximum angular separation for the underlying histograms to be fitted.

    num_bins : int
        The number of bins for the underlying histograms to be fitted.

    interactive : bool
        If True, show some diagnostic plots.
    """
    # Convert the coordinates in arcmin.
    delta_ra = degrees_to_arcmin(delta_ra)
    delta_dec = degrees_to_arcmin(delta_dec)
    # Create the histograms in right ascension and declination to fit the offset.
    binning = numpy.linspace(-rmax, rmax, num_bins)
    model_ra = xConstant() + xLorentzian()
    model_dec = xConstant() + xLorentzian()
    hist_ra = xHistogram1d(binning, xlabel='$\\Delta$ R. A. [arcmin]').fill(delta_ra)
    hist_dec = xHistogram1d(binning, xlabel='$\\Delta$ Dec. [arcmin]').fill(delta_dec)
    # Fit the histograms with a model.
    p0 = (0., 1000., 0., 0.1)
    try:
        fit_histogram(model_ra, hist_ra, p0=p0)
    except RuntimeError as e:
        logger.error('Cannot fit delta-ra histogram: %s', e)
    try:
        fit_histogram(model_dec, hist_dec, p0=p0)
    except RuntimeError as e:
        logger.error('Cannot fit delta-dec histogram: %s', e)
    # Convert back the fitted offset into arcmin.
    offset = arcmin_to_degrees(model_ra.Peak), arcmin_to_degrees(model_dec.Peak)
    logger.info('Offset best estimate: %s', offset)
    # If necessary, show the diagnostic plots.
    if interactive:
        plt.figure()
        hist_ra.plot()
        model_ra.plot()
        model_ra.stat_box()
        setup_gca(grids=True)
        plt.figure()
        hist_dec.plot()
        model_dec.plot()
        model_dec.stat_box()
        setup_gca(grids=True)
    return offset


def fit_psf_model(irf_name, du_id, radial_hist, rcore):
    """Calculate the psf radial model for a given IRF name and DU, fitter to the core
    of a given radial histogram.

    This is essentially loading the PSF for a given IRF name and Detector Unit,
    and fitting the corresponding radial profile

    Arguments
    ---------
    irf_name : str
        The IRF name.

    du_id : int
        The Detector Unit identifier.

    radial_hist : xHistogram1d instance
        The radial histogram containing the counts corrected for the solid angle.

    rmax : float
        The maximum angular separation for the fit.
    """
    psf = load_psf(irf_name, du_id)
    fit_model = lambda r, norm: norm * psf(arcmin_to_arcsec(r))
    r = radial_hist.bin_centers()
    mask = r < rcore
    r = r[mask]
    counts = radial_hist.content[mask]
    # Fit the PSF profile ro the inner core of the radial profile.
    popt, pcov = curve_fit(fit_model, r, counts, sigma=numpy.sqrt(counts))
    return lambda r: fit_model(r, popt)


def xpradialprofile(**kwargs):
    """Create and show the radial profile plot.
    """
    file_list = kwargs.get('filelist')
    for file_path in file_list:
        check_input_file(file_path, 'fits')
        file_name = os.path.basename(file_path)
        # Open the event file and retrieve the sky-position data.
        event_file = xEventFile(file_path)
        ra, dec = event_file.sky_position_data()
        # Retrieve the nominal center for the radial profile---this is by
        # default the reference pixel in the WCS of the input file, and the user
        # can ovverride this via the--ra and --dec command-line switches.
        obj_ra, obj_dec = event_file.wcs_reference()
        ra0, dec0 = kwargs.get('ra'), kwargs.get('dec')
        if ra0 is None:
            ra0 = obj_ra
        if dec0 is None:
            dec0 = obj_dec
        logger.info('Nominal center for radial plot set to (%.4f, %.4f)', ra0, dec0)
        # If necessary, optimize the center of the count map based on the
        # count map into the event file.
        if kwargs.get('autocenter'):
            logger.info('Recentering the map...')
            offset_ra, offset_dec = fit_offset(ra - ra0, dec - dec0,
                interactive=kwargs.get('interactive'))
            ra0 += offset_ra
            dec0 += offset_dec
            logger.info('Final center for radial profile set to (%.4f, %.4f)', ra0, dec0)
        # Calculate the angular separation and create the weighted radial histogram.
        angsep = degrees_to_arcmin(angular_separation(ra, dec, ra0, dec0))
        binning = numpy.linspace(kwargs.get('rmin'), kwargs.get('rmax'), kwargs.get('rbins'))
        hist = xHistogram1d(binning, xlabel='Radial distance [arcmin]')
        hist.fill(angsep, weights=1. / angsep)
        figure_name = '%s radial profile' % file_name
        if kwargs.get('residuals'):
            ax1, ax2 = residual_plot(figure_name)
        else:
            plt.figure(figure_name)
        hist.plot()
        # If necessary, overlay the PSF profile.
        if kwargs.get('psf'):
            irf_name = kwargs.get('irfname')
            du_id = event_file.du_id()
            rcore = kwargs.get('rcore')
            # Retrieve the PSF model fitted to the core of the count ditribution.
            psf_model = fit_psf_model(irf_name, du_id, hist, rcore)
            r = hist.bin_centers()
            plt.plot(r, psf_model(r), label='PSF %s (DU %s)' % (irf_name, du_id))
            # If necessary, plot the residual background.
            if kwargs.get('residuals'):
                plt.sca(ax2)
                res = (hist.content - psf_model(r))
                sigma_res = hist.errors()
                plt.errorbar(r, res, sigma_res, fmt='o')
                setup_gca(xlabel=hist.labels[0], ylabel='Residuals', grids=True,
                    ymin=0., ymax=3. * res.mean())
                plt.sca(ax1)
            # Now calculate the cumulative source counts as a function of the
            # angular separation based on the PSF model. Note that it is not
            # obvious that we can use the fitted PSF model directly, as
            # the normalization was fitted in the weighted (that is, corrected
            # for the solid angle) setting. We therefore resort to create a spline
            # with r * psf(r), normalize it so that the integral between 0 and
            # rcore matches the total number of total counts within rcore, and
            # than use it to integrate outside any desired radius.
            psf_spline = xInterpolatedUnivariateSpline(r, r * psf_model(r))
            norm = (angsep <= rcore).sum() / psf_spline.integral(0., rcore)
            psf_spline = psf_spline.scale(norm)
            # At this point we have a model for the radial cumulative source
            # counts that we can use to label the plot.
            for r0 in numpy.linspace(1., 5., 5):
                tot_counts = (angsep >= r0).sum()
                src_counts = psf_spline.integral(r0, 10.)
                src_ratio = src_counts / tot_counts
                print(r0, tot_counts, src_counts, src_counts/tot_counts)
                y = psf_model(r0)
                plt.plot(r0, y, 'o', color=last_line_color())
                plt.text(r0, y, ' %.2f%%' % (src_ratio * 100), color=last_line_color())
        setup_gca(logy=True, grids=True, legend=True)
    plt.show()


def main():
    """main() entry point.
    """
    xpradialprofile(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()
