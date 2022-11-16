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
from ixpeobssim.core.modeling import xLine
from ixpeobssim.core.fitting import fit_histogram, fit
from ixpeobssim.core.rand import xUnivariateGenerator
from ixpeobssim.srcmodel.bkg import xRadialBackgroundGenerator
from ixpeobssim.utils.logging_ import logger
from ixpeobssim.utils.matplotlib_ import plt, setup_gca


if sys.flags.interactive:
    plt.ion()






class TestRadialBackground(unittest.TestCase):

    """Unit test for radial background.
    """

    def test_oversample(self, half_size=7., num_events=500000):
        """
        """
        slope, oversample = [], []
        for s in numpy.linspace(-1., 2., 11):
            gen = xRadialBackgroundGenerator(half_size, half_size, s)
            r = gen.rvs(num_events)
            phi = 2. * numpy.pi * numpy.random.random(num_events)
            x, y = gen.trim(*gen.polar_to_cartesian(r, phi))
            ovrsmpl = len(r) / len(x)
            logger.info('Oversampling for %.3f radial slope: %.5f', s, ovrsmpl)
            slope.append(s)
            oversample.append(ovrsmpl)
        plt.figure('Radial background oversampling')
        plt.plot(slope, oversample, 'o')
        setup_gca(xlabel='Radial slope', ylabel='Oversampling factor', grids=True)

    def _test_direct(self, half_size=7., num_events=10000000, alpha=0.2):
        """Convenience function for the direct transform (radial -> xy).
        """
        radius = numpy.sqrt(2.) * half_size
        x, y = xRadialBackgroundGenerator(half_size, half_size, alpha).rvs_xy(num_events)
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
        alpha_hat = model.Slope * half_size / model(0.5 * half_size)
        logger.info('Best-fit slope: %.3f (target %.3f)', alpha_hat, alpha)

        def pdf(x, N, b):
            return N * (1. + b * abs(x)**1.5)

        plt.figure('Radial background x')
        binning = numpy.linspace(-half_size, half_size, 100)
        hist = xHistogram1d(binning, xlabel='x').fill(x)
        hist.plot()
        _x = hist.bin_centers()
        popt, pcov = curve_fit(pdf, _x, hist.content)
        print(popt)
        plt.plot(_x, pdf(_x, *popt))
        plt.figure('Radial background y')
        hist = xHistogram1d(binning, xlabel='y').fill(y)
        hist.plot()
        _x = hist.bin_centers()
        popt, pcov = curve_fit(pdf, _x, hist.content)
        print(popt)
        plt.plot(_x, pdf(_x, *popt))
        plt.figure('Radial background xy')
        hist = xHistogram2d(binning, binning, xlabel='x', ylabel='y').fill(x, y)
        hist.plot()



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)
