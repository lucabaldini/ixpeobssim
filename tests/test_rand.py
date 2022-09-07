#!/usr/bin/env python
#
# Copyright (C) 2019, the ixpeobssim team.
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

"""Unit test for the core.rand module.
"""

from __future__ import print_function, division

import sys
import unittest

import numpy
numpy.random.seed(667)

from ixpeobssim.core.rand import xUnivariateGenerator, xUnivariateAuxGenerator
from ixpeobssim.core.hist import xHistogram1d
from ixpeobssim.utils.matplotlib_ import plt, last_line_color
from ixpeobssim.utils.logging_ import logger

if sys.flags.interactive:
    plt.ion()


class TestRandUnivariate(unittest.TestCase):

    """Unit test for univariate random generator.
    """

    def _basic_test(self, name, x, y, size=1000000, num_bins=100):
        """Convenience method to streamline the basic test with arbitrary data
        points.

        This is creating an instance of a xUnivariateGenerator, generating a
        random sample of the appropriate size, and calculating a chisquare
        with a completely determined model.
        """
        logger.info('Random test: %s', name)
        # Create the random generator based on the input arrays.
        g = xUnivariateGenerator(x, y)
        # Extract the random numbers and build a histogram.
        r = g.rvs(size)
        h = xHistogram1d(numpy.linspace(x.min(), x.max(), num_bins + 1)).fill(r)
        # Poor's man proxy for the chisquare (note we are not integrating the
        # model, just calculating it at the bin centers).
        _xc = h.bin_centers(0)
        exp = g(_xc) * size * (x.max() - x.min()) / num_bins
        obs = h.content
        chisq = ((obs - exp)**2. / exp).sum()
        delta = (chisq - num_bins) / numpy.sqrt(2. * num_bins)
        logger.info('Chisq = %.3f / %d dof (%.3f sigma)', chisq, num_bins, delta)
        self.assertTrue(abs(delta) < 4.)
        plt.figure('Test: %s' % name)
        h.plot()
        plt.plot(_xc, exp)

    def test_uniform(self):
        """Test a uniform distribution.
        """
        n = 100
        _x = numpy.linspace(0., 1., n)
        _y = numpy.full(n, 1.)
        self._basic_test('uniform', _x, _y)

    def test_triangle(self):
        """Test a triangular distribution.
        """
        _x = numpy.linspace(0., 1., 100)
        _y = 2 * _x
        self._basic_test('triangle', _x, _y)

    def test_gauss(self):
        """Test a gaussian distribution.
        """
        _x = numpy.linspace(-5., 5., 1000)
        _y = 1. / numpy.sqrt(2. * numpy.pi) * numpy.exp(-_x**2. / 2.)
        self._basic_test('gauss', _x, _y)

    def test_negative(self):
        """Test a distribution with negative weights.
        """
        _x = numpy.linspace(0., 1., 100)
        _y = 2 * _x - 1.
        with self.assertRaises(SystemExit):
            self._basic_test('negative triangle', _x, _y)



class TestRandAuxUnivariateGaussian(unittest.TestCase):

    """Unit test for univariate random generators with an auxiliary variable.
    """

    @classmethod
    def setUpClass(cls):
        """Setup: create a univariate generator depending on an auxiliary
        variable.
        """
        cls.rv = numpy.linspace(-5., 5., 500)
        cls.aux = numpy.linspace(0., 1., 100)
        cls.generator = xUnivariateAuxGenerator(cls.rv, cls.aux, cls._aux_gauss)

    @staticmethod
    def _aux_gauss(rv, aux):
        """Small utility function: a gauassian where the mean and sigma depend
        on an auxiliary variable (and in a fairly nonsense fashion, too).

        We assume that aux lies in the [0, 1] interval.
        """
        mu = aux
        sigma = 1. - 0.25 * aux
        norm = 1. / numpy.sqrt(2 * numpy.pi) / sigma
        return norm * numpy.exp(-(rv - mu)**2. / (2 * sigma**2.))

    def test_negative(self):
        """Test a negative pdf.
        """
        rv = numpy.linspace(-5., 5., 500)
        aux = numpy.linspace(0., 1., 100)
        z = lambda x, y: numpy.full((rv.size, aux.size), -1.)
        with self.assertRaises(SystemExit):
            generator = xUnivariateAuxGenerator(rv, aux, z)

    def test_pdf(self):
        """Plot the 2d pdf and a few slices.
        """
        plt.figure('Gaussian aux generator')
        self.generator.plot()
        plt.figure('Gaussian aux generator slices')
        for _aux in (1., 0.5, 0.):
            s = self.generator.slice(_aux)
            s.plot(label='Slice @ aux=%.3f' % _aux)
            delta = abs(s(self.rv) - self._aux_gauss(self.rv, _aux))
            self.assertTrue(delta.max() < 1.e-9)
        plt.legend()

    def test_ppf(self):
        """Plot the 2d ppf.
        """
        ppf = self.generator.ppf
        plt.figure('Gaussian aux ppf')
        ppf.plot()

    def test_rvs(self, size=1000000, num_bins=200, min_entries=5):
        """Test throwing out random numbers.
        """
        logger.info('Testing throwing random numbers...')
        binning = numpy.linspace(self.rv.min(), self.rv.max(), num_bins + 1)
        plt.figure('Gaussian aux rvs')
        for _aux in (1.0, 0.5, 0.):
            logger.info('Auxiliary variable: %.5f', _aux)
            r = self.generator.rvs(numpy.full((size,), _aux))
            h = xHistogram1d(binning).fill(r)
            s = self.generator.slice(_aux)
            h.plot()
            scale = size / num_bins * (self.rv.max() - self.rv.min())
            s.plot(label='Slice @ aux=%.3f' % _aux, scale=scale)
            # Calculate the chisquare.
            _xc = h.bin_centers(0)
            exp = s(_xc) * scale
            # Mask the bins with too few entries, as a single bin in the tail
            # with one entry can screw up the chisquare completely.
            mask = exp >= min_entries
            _xc = _xc[mask]
            exp = exp[mask]
            obs = h.content[mask]

            def debug_residuals():
                """Small convenience function to debug issue #244.
                """
                res = (obs - exp) / numpy.sqrt(exp)
                plt.figure('Residuals at aux %.1f' % _aux)
                plt.plot(_xc, res)
                plt.xlabel('rv')
                plt.ylabel('Relative residuals')

            #debug_residuals()

            nu = len(obs)
            chisq = ((obs - exp)**2. / exp).sum()
            delta = (chisq - nu) / numpy.sqrt(2. * nu)
            logger.info('Chisq = %.3f / %d dof (%.3f sigma)', chisq, num_bins, delta)
            self.assertTrue(abs(delta) < 5., 'delta = %.3f' % delta)
        plt.legend()



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)
