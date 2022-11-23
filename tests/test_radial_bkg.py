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

"""Unit tests for the implementation of the radial background.
"""


from __future__ import print_function, division

import sys
import unittest

import numpy
from scipy.optimize import curve_fit

from ixpeobssim.core.hist import xHistogram1d, xHistogram2d
from ixpeobssim.core.modeling import xLine, xConstant
from ixpeobssim.core.fitting import fit_histogram, fit
from ixpeobssim.core.rand import xUnivariateGenerator
from ixpeobssim.instrument.gpd import GPD_PHYSICAL_MAX_RADIUS, gpd_map_binning,\
    GPD_PHYSICAL_HALF_SIDE_X, GPD_PHYSICAL_HALF_SIDE_Y, within_gpd_physical_area
from ixpeobssim.srcmodel.bkg import xRadialBackgroundGenerator
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca


if sys.flags.interactive:
    plt.ion()



class TestRadialBackground(unittest.TestCase):

    """Unit test for radial background.
    """

    def test_oversample(self, num_events=500000):
        """Small snippet to gauge the heuristic for the necessary oversample factor.
        """
        slope, oversample = [], []
        for s in numpy.linspace(-1., 2., 11):
            gen = xRadialBackgroundGenerator(s)
            r = gen.rvs(num_events)
            phi = 2. * numpy.pi * numpy.random.random(num_events)
            x, y = gen.polar_to_cartesian(r, phi)
            mask = within_gpd_physical_area(x, y)
            x = x[mask]
            y = y[mask]
            ovrsmpl = len(r) / len(x)
            logger.info('Oversampling for %.3f radial slope: %.5f', s, ovrsmpl)
            slope.append(s)
            oversample.append(ovrsmpl)
        plt.figure('Radial background oversampling')
        plt.plot(slope, oversample, 'o')

        def pol2(x, c0, c1, c2):
            return c0 + c1 * x + c2 * x**2.

        popt, pcov = curve_fit(pol2, slope, oversample)
        logger.info('Oversample best fit quadratic: %s', popt)
        xgrid = numpy.linspace(-1., 2., 100)
        plt.plot(xgrid, pol2(xgrid, *popt))
        setup_gca(xlabel='Radial slope', ylabel='Oversampling factor', grids=True)

    def test_sample(self, num_events=1000000, slope=0.2):
        """Convenience function for the direct transform (radial -> xy).
        """
        half_size = xRadialBackgroundGenerator.HALF_SIDE
        radius = numpy.sqrt(2.) * half_size
        x, y = xRadialBackgroundGenerator(slope).rvs_xy(num_events)
        r = numpy.sqrt(x**2. + y**2.)
        plt.figure('Radial background r')
        binning = numpy.linspace(0., radius, 100)
        rmin = half_size / 10.
        mask = r > rmin
        hist = xHistogram1d(binning, xlabel='r').fill(r[mask], weights=1. / r[mask])
        model = fit_histogram(xLine(), hist, xmax=half_size)
        model.plot()
        model.stat_box()
        hist.plot()
        slope_hat = model.Slope * half_size / model(0.5 * half_size)
        logger.info('Best-fit slope: %.3f (target %.3f)', slope_hat, slope)
        plt.figure('Radial background xy')
        xbinning, ybinning = gpd_map_binning(GPD_PHYSICAL_HALF_SIDE_X,
            GPD_PHYSICAL_HALF_SIDE_Y, 25)
        hist = xHistogram2d(xbinning, ybinning, xlabel='x', ylabel='y').fill(x, y)
        hist.plot()

    def test_sample_constant(self, num_events=1000000, num_bins=100):
        """Test the radial generator with zero slope.
        """
        half_size = xRadialBackgroundGenerator.HALF_SIDE
        xbinning, ybinning = gpd_map_binning(GPD_PHYSICAL_HALF_SIDE_X,
            GPD_PHYSICAL_HALF_SIDE_Y, num_bins)
        x, y = xRadialBackgroundGenerator(0.).rvs_xy(num_events)
        plt.figure('Constant radial background x')
        hist = xHistogram1d(xbinning, xlabel='x').fill(x)
        model = fit_histogram(xConstant(), hist)
        hist.plot()
        model.plot()
        model.stat_box()
        self.assertTrue(model.chisq < model.ndof + 5. * numpy.sqrt(2. * model.ndof))
        plt.figure('Constant radial background y')
        hist = xHistogram1d(ybinning, xlabel='y').fill(y)
        model = fit_histogram(xConstant(), hist)
        hist.plot()
        model.plot()
        model.stat_box()
        self.assertTrue(model.chisq < model.ndof + 5. * numpy.sqrt(2. * model.ndof))



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)
