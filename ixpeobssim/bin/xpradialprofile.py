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
from ixpeobssim.core.modeling import xFitModelBase, xConstant, xGaussian
from ixpeobssim.evt.event import xEventFile
from ixpeobssim.irf import load_psf
from ixpeobssim.utils.argparse_ import xArgumentParser
from ixpeobssim.utils.astro import angular_separation
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.utils.os_ import check_input_file
from ixpeobssim.utils.units_ import degrees_to_arcmin, arcmin_to_degrees, arcmin_to_arcsec


__description__ = \
"""Create a radial profile plot, in sky coordinated, for a given observation.

This program takes a level-2 file as an input
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
PARSER.add_argument('--rbins', type=int, default=200,
                    help='number of bins for the radial profile')
PARSER.add_boolean('--recenter', default=True,
                   help='recenter the image based on the count map')
PARSER.add_boolean('--psf', default=True,
                   help='overimpose the PSF radial profile')
PARSER.add_irfname()



def calculate_center(delta_ra, delta_dec, rmax=2., num_bins=100, interactive=False):
    """
    """
    logger.info('Calculating the center of the count map in sky coordinates...')
    delta_ra = degrees_to_arcmin(delta_ra)
    delta_dec = degrees_to_arcmin(delta_dec)
    binning = numpy.linspace(-rmax, rmax, num_bins)
    model_ra = xConstant() + xGaussian()
    model_dec = xConstant() + xGaussian()
    hist_ra = xHistogram1d(binning, xlabel='$\\Delta$ R. A. [arcmin]').fill(delta_ra)
    hist_dec = xHistogram1d(binning, xlabel='$\\Delta$ Dec. [arcmin]').fill(delta_dec)
    fit_histogram(model_ra, hist_ra)
    fit_histogram(model_dec, hist_dec)
    offset = arcmin_to_degrees(model_ra.Peak), arcmin_to_degrees(model_dec.Peak)
    logger.info('Offset best estimate: %s', offset)
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


def xpradialprofile(**kwargs):
    """
    """
    file_list = kwargs.get('filelist')
    for file_path in file_list:
        check_input_file(file_path, 'fits')
        file_name = os.path.basename(file_path)
        event_file = xEventFile(file_path)
        ra, dec = event_file.sky_position_data()



        # Retrieve the center for the radial profile.
        obj_ra, obj_dec = event_file.wcs_reference()
        ra0, dec0 = kwargs.get('ra'), kwargs.get('dec')
        if ra0 is None:
            ra0 = obj_ra
        if dec0 is None:
            dec0 = obj_dec
        logger.info('Center for radial plot set to (%.3f, %.3f)', ra0, dec0)
        if kwargs.get('recenter'):
            logger.info('Recentering the map...')
            offset_ra, offset_dec = calculate_center(ra - ra0, dec - dec0)
            ra0 += offset_ra
            dec0 += offset_dec
        angsep = degrees_to_arcmin(angular_separation(ra, dec, ra0, dec0))
        binning = numpy.linspace(kwargs.get('rmin'), kwargs.get('rmax'), kwargs.get('rbins'))
        hist = xHistogram1d(binning, xlabel='Radial distance [arcmin]')
        hist.fill(angsep, weights=1. / angsep)
        plt.figure('%s radial profile' % file_name)
        hist.plot()
        if kwargs.get('psf'):
            du_id = event_file.du_id()
            psf = load_psf(kwargs.get('irfname'), du_id)
            r = hist.bin_centers()
            model = lambda r, norm: norm * psf(arcmin_to_arcsec(r))
            mask = r < 0.75
            counts = hist.content[mask]
            popt, pcov = curve_fit(model, r[mask], counts, sigma=numpy.sqrt(counts))
            label = 'PSF %s (DU %s)' % (kwargs.get('irfname'), du_id)
            plt.plot(r, popt * psf(arcmin_to_arcsec(r)), label=label)
        setup_gca(logy=True, grids=True, legend=True)
    plt.show()



def main():
    """main() entry point.
    """
    xpradialprofile(**PARSER.parse_args().__dict__)



if __name__ == '__main__':
    main()
