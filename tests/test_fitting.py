#!/usr/bin/env python
#
# Copyright (C) 2021, the ixpeobssim team.
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

"""Unit tests for the core.fitting module.
"""


from __future__ import print_function, division

import unittest
import sys

import numpy

from ixpeobssim.core.fitting import linear_analytical_fit, power_law_analytical_fit
from ixpeobssim.utils.matplotlib_ import plt, setup_gca
if sys.flags.interactive:
    plt.ion()


class TestFitting(unittest.TestCase):

    """Unit test for the core.fitting module.
    """

    @classmethod
    def setUpClass(cls):
        """Test setup.
        """
        cls.xmin = 1.
        cls.xmax = 10.
        cls.n = 10
        cls.npars = 5
        cls.x = numpy.linspace(cls.xmin, cls.xmax, cls.n)
        cls.x2 = numpy.random.uniform(cls.xmin, cls.xmax, size=(cls.n, cls.npars))
        cls.xgrid = numpy.linspace(cls.xmin, cls.xmax, 100)
        cls.ngrid = len(cls.xgrid)

    @staticmethod
    def _setup_gca(**kwargs):
        """
        """
        setup_gca(xlabel='x [a.u.]', ylabel='y [a.u.]', grids=True, **kwargs)

    @staticmethod
    def _broadcast(func, x, *params):
        """
        """
        params = [numpy.tile(param, (len(x), 1)).T for param in params]
        return func(x, *params).T

    def test_linear1(self):
        """Simple linear fit, with one dimensional x and y.
        """
        m = 2.
        q = 1.
        y = m * self.x + q
        m0, q0 = linear_analytical_fit(self.x, y)
        self.assertAlmostEqual(m0, m)
        self.assertAlmostEqual(q0, q)
        plt.figure('Linear 1')
        plt.plot(self.x, y, 'o')
        plt.plot(self.xgrid, m0 * self.xgrid + q0)
        self._setup_gca()

    def test_linear2(self):
        """Linear fit with a common x and multiple y vectors.
        """
        m = numpy.linspace(1., 5., self.npars)
        q = numpy.full(self.npars, 0.)
        f = lambda x, m, q: m * x + q
        y = self._broadcast(f, self.x, m, q)
        m0, q0 = linear_analytical_fit(self.x, y)
        self.assertTrue(numpy.allclose(m0, m))
        self.assertTrue(numpy.allclose(q0, q))
        plt.figure('Linear 2')
        plt.plot(self.x, y, 'o')
        plt.plot(self.xgrid, self._broadcast(f, self.xgrid, m0, q0))
        self._setup_gca()

    def test_linear3(self):
        """Linear fit with multiple x and multiple y vectors.
        """
        m = numpy.linspace(1., 5., self.npars)
        q = numpy.full(self.npars, 0.)
        f = lambda x, m, q: m * x + q
        y = f(self.x2, m, q)
        m0, q0 = linear_analytical_fit(self.x2, y)
        self.assertTrue(numpy.allclose(m0, m))
        self.assertTrue(numpy.allclose(q0, q))
        plt.figure('Linear 3')
        plt.plot(self.x2, y, 'o')
        plt.plot(self.xgrid, self._broadcast(f, self.xgrid, m0, q0))
        self._setup_gca()

    def test_power_law1(self):
        """Test a simple power-law fit.
        """
        norm = 10.
        index = -2.
        y = norm * self.x**index
        norm0, index0 = power_law_analytical_fit(self.x, y)
        self.assertAlmostEqual(norm0, norm)
        self.assertAlmostEqual(index0, index)
        plt.figure('Power law 1')
        plt.plot(self.x, y, 'o')
        plt.plot(self.xgrid, norm0 * self.xgrid**index0)
        self._setup_gca(logx=True, logy=True)

    def test_power_law2(self):
        """Power-law fit with single x and multiple y vectors.
        """
        norm = numpy.full(self.npars, 1.)
        index = numpy.linspace(-2., -1., self.npars)
        f = lambda x, norm, index: norm * x**index
        y = self._broadcast(f, self.x, norm, index)
        norm0, index0 = power_law_analytical_fit(self.x, y)
        self.assertTrue(numpy.allclose(norm0, norm))
        self.assertTrue(numpy.allclose(index0, index))
        plt.figure('Power law 2')
        plt.plot(self.x, y, 'o')
        plt.plot(self.xgrid, self._broadcast(f, self.xgrid, norm0, index0))
        self._setup_gca(logx=True, logy=True)

    def test_power_law3(self):
        """Power-law fit with multiple x and multiple y vectors.
        """
        norm = numpy.full(self.npars, 1.)
        index = numpy.linspace(-2., -1., self.npars)
        f = lambda x, norm, index: norm * x**index
        y = f(self.x2, norm, index)
        norm0, index0 = power_law_analytical_fit(self.x2, y)
        self.assertTrue(numpy.allclose(norm0, norm))
        self.assertTrue(numpy.allclose(index0, index))
        plt.figure('Power law 3')
        plt.plot(self.x2, y, 'o')
        plt.plot(self.xgrid, self._broadcast(f, self.xgrid, norm0, index0))
        self._setup_gca(logx=True, logy=True)



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)
