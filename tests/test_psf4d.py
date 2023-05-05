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

from __future__ import print_function, division


"""Unit test for the 2d part of the irf.psf module.
"""

import os
import sys
import unittest

import numpy
from scipy.stats import gmean

from ixpeobssim import IXPEOBSSIM_IRFGEN_DATA
from ixpeobssim.core.hist import xHistogram1d, xHistogram2d
from ixpeobssim.instrument import DU_IDS
from ixpeobssim.instrument.mma import FOCAL_LENGTH
from ixpeobssim.irf.caldb import irf_file_path
from ixpeobssim.irf.psf import xPointSpreadFunction, xPointSpreadFunction2d,\
    xPointSpreadFunction4d
from ixpeobssim.irf import load_psf, DEFAULT_IRF_NAME
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
from ixpeobssim.utils.units_ import degrees_to_arcmin
from ixpeobssim.utils.logging_ import logger

if sys.flags.interactive:
    plt.ion()


class TestIxpePsf4d(unittest.TestCase):

    """Unit test for the IXPE point-spread function in the 2d version.
    """

    def test_1d(self):
        """Compare the underlying 1-dimensional splines for version 12 and 14 of
        the PSF.
        """
        for du_id in DU_IDS:
            psf12 = xPointSpreadFunction(irf_file_path('ixpe:obssim:v12', du_id, 'psf'))
            psf14 = xPointSpreadFunction(irf_file_path('ixpe:obssim:v14', du_id, 'psf'))
            self.assertTrue(numpy.allclose(psf12.x, psf14.x))
            self.assertTrue(numpy.allclose(psf12.y, psf14.y))

    def test_loading_error(self, irf_name='ixpe:obssim:v12'):
        """Make sure that for the IRF v12 the 2d and 4d PSF cannot be loaded.
        """
        for du_id in DU_IDS:
            with self.assertRaises(KeyError) as context:
                psf = xPointSpreadFunction2d(irf_file_path(irf_name, du_id, 'psf'))
            logger.info(context.exception)
            with self.assertRaises(KeyError) as context:
                psf = xPointSpreadFunction4d(irf_file_path(irf_name, du_id, 'psf'))
            logger.info(context.exception)

    def test_2d(self, irf_name='ixpe:obssim:v14', num_samples=1000000, half_size=0.1):
        """Test the PSF 2D.

        This is plotting the underlying image and sampling the corresponding
        distribution.
        """
        for du_id in DU_IDS:
            psf = xPointSpreadFunction2d(irf_file_path(irf_name, du_id, 'psf'))
            plt.figure('PSF image %s DU %d' % (irf_name, du_id))
            psf.psf_image.plot(stretch='log', zlabel='PSF', vmin=1.e-7, vmax=1.e-2)
            x, y = psf.delta(num_samples)
            binning = numpy.linspace(-half_size, half_size, 250)
            h = xHistogram2d(binning, binning, xlabel='Delta R. A.', ylabel='Delta Dec.').fill(x, y)
            plt.figure('PSF sampling %s DU %d' % (irf_name, du_id))
            h.plot(logz=True)
            plt.gca().set_aspect('equal')

    def test_4d_rscaling(self, irf_name='ixpe:obssim:v14'):
        """Test the PSF 4D.
        """
        r = numpy.linspace(0., 400., 401)
        for du_id in DU_IDS:
            psf = xPointSpreadFunction4d(irf_file_path(irf_name, du_id, 'psf'))
            plt.figure('PSF 4d radial scaling %s DU %d' % (irf_name, du_id))
            for theta in (0., 2., 5.):
                for energy in (1., 2., 8., 12.):
                    scale = psf.rscale_interpolator((theta, energy, r))
                    label = '%.1f arcmin, %.2f keV' % (theta, energy)
                    plt.plot(r, scale, label=label)
            setup_gca(grids=True, legend=True, xlabel='radius [arcsec]',
                      ylabel='EEF scaling')

    def test_4d_sampling(self, irf_name='ixpe:obssim:v14', num_samples=1000000,
                         half_size=8.029536666666667):
        """Sample the 4-dimensional PSF.
        """
        x_binning = numpy.linspace(-half_size, half_size, 722)
        r_bins = numpy.logspace(-1.3, 1.1, 21)

        for du_id in DU_IDS:
            psf = xPointSpreadFunction4d(irf_file_path(irf_name, du_id, 'psf'))

            def make_2d_hist(energy, theta):
                x, y = psf.delta(len(energy), energy, theta, num_samples)
                h = xHistogram2d(x_binning, x_binning,
                                 xlabel='Delta R. A.', ylabel='Delta Dec.')
                h.fill(degrees_to_arcmin(x), degrees_to_arcmin(y))
                h_norm = h * (1. / h.sum())
                return h_norm

            def radial_profile(h):
                x, y = h.bin_centers(axis=0), h.bin_centers(axis=1)
                vals = []
                radii = []
                errs = []
                for r_low, r_up in zip(r_bins[:-1], r_bins[1:]):
                    mask = (x**2 + y**2 >= r_low**2) * (x**2 + y**2 < r_up**2)
                    radii.append(gmean([r_low, r_up]))
                    counts = numpy.sum(h.entries[mask])
                    norm_counts = numpy.sum(h.content[mask])
                    vals.append(norm_counts)
                    errs.append(numpy.sqrt(counts) * norm_counts / counts)
                return numpy.array(radii), numpy.array(vals), numpy.array(errs)

            ref_h = make_2d_hist(2.29, 0.)
            _, reference_profile, _ = radial_profile(ref_h)
            for theta, energy in ((0., 2.29), (0., 4.51), (0., 6.4), (4., 6.4)):
                label = '%s DU %d, %.1f arcsec, %.2f keV' %\
                        (irf_name, du_id, theta, energy)
                r, vals, errs = radial_profile(make_2d_hist(energy,theta))
                vals /= reference_profile
                errs /= reference_profile
                plt.figure('radial profile MMA %d' % du_id)
                plt.errorbar(r, vals, yerr=errs, fmt='o--', label=label)
                setup_gca(grids=True, logx=True, legend=True,
                          xlabel='radius [arcmin]',
                          ylabel='Radial profile / reference')


if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)
