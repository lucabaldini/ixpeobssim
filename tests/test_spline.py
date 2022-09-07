#!/usr/bin/env python
#
# Copyright (C) 2015--2019, the ixpeobssim team.
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


"""Unit test for the core.spline module.
"""

import sys

import unittest
import scipy.special
import numpy

from ixpeobssim.core.spline import xInterpolatedUnivariateSpline
from ixpeobssim.core.spline import xInterpolatedUnivariateSplineLinear
from ixpeobssim.core.spline import xInterpolatedUnivariateLogSpline
from ixpeobssim.core.spline import xInterpolatedUnivariateLogSplineLinear
from ixpeobssim.core.spline import xInterpolatedBivariateSpline

from ixpeobssim.utils.matplotlib_ import plt, last_line_color

if sys.flags.interactive:
    plt.ion()


class TestUnivariateSpline(unittest.TestCase):

    """Unit test for univariate splines.
    """

    @classmethod
    def setUpClass(cls):
        """Setup.

        Create a few objects to be used for testing.
        """
        cls.num_points = 100
        cls.x1 = numpy.linspace(0., 2. * numpy.pi, cls.num_points)
        cls.y1 = numpy.sin(cls.x1)
        cls.x2 = numpy.linspace(0., numpy.pi, cls.num_points)
        cls.y2 = numpy.sin(cls.x2)
        cls.x3 = numpy.linspace(0., 10., 100)
        cls.y3 = 3. * cls.x3
        fmt = dict(xlabel='x [a. u.]', ylabel='y [a. u.]')
        cls.s1 = xInterpolatedUnivariateSplineLinear(cls.x1, cls.y1, **fmt)
        cls.s2 = xInterpolatedUnivariateSplineLinear(cls.x2, cls.y2, **fmt)
        cls.s3 = xInterpolatedUnivariateSplineLinear(cls.x3, cls.y3, **fmt)
        cls.s4 = xInterpolatedUnivariateSpline(cls.x1, cls.y1, **fmt)

    def test_len(self):
        """Test the basic object instantiation.
        """
        # Check we get the number of points right.
        self.assertEqual(len(self.s1), self.num_points)

    def test_evaluation(self):
        """Test the object evaluation.
        """
        # This is a linear interpolator, so the interpolated values must
        # be identical, within rounding errors, to the original grid of
        # values.
        _delta = abs(self.s1(self.x1) - self.y1)
        self.assertTrue(_delta.max() < 1e-9, 'max. diff. %.9f' % _delta.max())

        # s1 and s2 are built with different sets of points, but with the same
        # underlying function, so they should be fairly close at any given
        # point.
        _x = numpy.linspace(0, numpy.pi, 10)
        _delta = abs(self.s1(_x) - self.s2(_x))
        self.assertTrue(_delta.max() < 1e-3, 'max. diff. %.9f' % _delta.max())

        # And now test with a spline with k = 3.
        _delta = abs(self.s4(self.x1) - self.y1)
        self.assertTrue(_delta.max() < 1e-9, 'max. diff. %.9f' % _delta.max())

    def test_multiplication(self):
        """Test the interpolator multiplication.
        """
        # Evaluate s1*s2 in x2 should give the same answer than multiplying
        # s1(x2)*y2.
        _m = self.s1 * self.s2
        _delta = abs(_m(self.x2) - self.s1(self.x2)*self.y2)
        self.assertTrue(_delta.max() < 1e-9, 'max. diff. %.9f' % _delta.max())

        # And the result of the multiplication should be an instance of
        # the original operand class.
        self.assertTrue(isinstance(_m, xInterpolatedUnivariateSplineLinear))

        # Make sure the labels are right.
        self.assertTrue(_m.xlabel == self.s1.xlabel)
        self.assertTrue(_m.ylabel == None)

    def test_division(self):
        """Test the interpolator division
        """
        # Evaluate s1/s2 in x2 should give the same answer than dividing
        # s1(x2) / y2.
        _d = self.s1 / self.s2
        # Small workaround to avoid zero-division errors.
        mask = self.y2 != 0.
        _delta = abs(_d(self.x2)[mask] - self.s1(self.x2)[mask] / self.y2[mask])
        self.assertTrue(_delta.max() < 1e-9, 'max. diff. %.9f' %_delta.max())

        # And the result of the division should be an instance of
        # the original operand class.
        self.assertTrue(isinstance(_d, xInterpolatedUnivariateSplineLinear))

        # Make sure the labels are right.
        self.assertTrue(_d.xlabel == self.s1.xlabel)
        self.assertTrue(_d.ylabel == None)

    def test_sum(self):
        """Test the interpolator sum.
        """
        # Evaluate s1 + s2 in x2 should give the same answer than adding
        # s1(x2) + y2.
        _s = self.s1 + self.s2
        _delta = abs(_s(self.x2) - (self.s1(self.x2) + self.y2))
        self.assertTrue(_delta.max() < 1e-9, 'max. diff. %.9f' % _delta.max())

        # And the result of the multiplication should be an instance of
        # the original operand class.
        self.assertTrue(isinstance(_s, xInterpolatedUnivariateSplineLinear))

        # Make sure the labels are right.
        self.assertTrue(_s.xlabel == self.s1.xlabel)
        self.assertTrue(_s.ylabel == self.s1.ylabel)

    def test_subtraction(self):
        """Test the interpolator subtraction.
        """
        _s = self.s1 - self.s2
        _delta = abs(_s(self.x2) - (self.s1(self.x2) - self.y2))
        self.assertTrue(_delta.max() < 1e-9, 'max. diff. %.9f' % _delta.max())

        # And the result of the multiplication should be an instance of
        # the original operand class.
        self.assertTrue(isinstance(_s, xInterpolatedUnivariateSplineLinear))

        # Make sure the labels are right.
        self.assertTrue(_s.xlabel == self.s1.xlabel)
        self.assertTrue(_s.ylabel == self.s1.ylabel)

    def test_extrapolation(self):
        """Test interpolator extrapolation.
        """
        # Calculate one extrapolated value by hand and compare it to the
        # value from the interpolator.
        _xa = self.x1[-2]
        _xb = self.x1[-1]
        _ya = self.y1[-2]
        _yb = self.y1[-1]
        _x = _xb + 0.2
        _y = _ya + (_yb - _ya)/(_xb - _xa)*(_x - _xa)
        _delta = abs(self.s1(_x) - _y)
        self.assertTrue(_delta < 1e-9, 'max. diff. %.9f' % _delta)

    def test_scale(self, scale_factor=2.):
        """Test the scaling.
        """
        _s = self.s1.scale(scale_factor)
        _delta = abs(_s.y - scale_factor * self.s1.y)
        self.assertTrue(_delta.max() < 1e-9, 'max. diff. %.9f' % _delta.max())

        # Make sure the labels are right.
        self.assertTrue(_s.xlabel == self.s1.xlabel)
        self.assertTrue(_s.ylabel == self.s1.ylabel)

    def test_norm(self):
        """Test the normalization calculation.
        """
        _delta = abs(self.s3.norm() - 100.*3./2)
        self.assertTrue(_delta < 1e6, 'norm. diff. %.9f' % _delta)

    def test_cdf(self):
        """ The cdf must be 0 at xmin and 1 at xmax.
        """
        cdf = self.s3.build_cdf()
        _delta = abs(cdf(self.s3.xmin()))
        self.assertTrue(_delta < 1e-3, 'cdff(xmin) %.9f' % _delta)
        _delta = abs(cdf(self.s3.xmax()) - 1.)
        self.assertTrue(_delta < 1e-3, 'cdf(xmax) - 1 %.9f' % _delta)
        plt.figure('Cumulative distribution function')
        cdf.plot()

    def test_ppf(self):
        """ The ppf must be defined between 0 and 1 (where is equal to the
        xmin and xmax values of the original spline).
        """
        ppf = self.s3.build_ppf()
        _delta = abs(ppf.xmin())
        self.assertTrue(_delta < 1e-3, 'ppf xmin %.9f' % _delta)
        _delta = abs(ppf.xmax() - 1.)
        self.assertTrue(_delta < 1e-3, 'ppf (xmax - 1) %.9f' % _delta)
        _delta = abs(ppf(0) - self.s3.xmin())
        self.assertTrue(_delta < 1e-3, 'ppf(0) - xmin %.9f' % _delta)
        _delta = abs(ppf(1) - self.s3.xmax())
        self.assertTrue(_delta < 1e-3, 'ppf(1) - xmax %.9f' % _delta)
        plt.figure('Percent-point function')
        ppf.plot()

    def test_cdf_erf(self):
        """Test the cdf for a gaussian function.
        """
        _x = numpy.linspace(-5, 5, 100)
        _y = 1. / numpy.sqrt(2. * numpy.pi) *numpy.exp(-0.5 * _x**2)
        pdf = xInterpolatedUnivariateSplineLinear(_x, _y, xlabel='x [a. u.]')
        cdf = pdf.build_cdf()
        delta = abs(cdf(_x) - 0.5*(1. + scipy.special.erf(_x/numpy.sqrt(2.))))
        max_delta = delta.max()
        err_msg = 'maximum absolute delta %.4e' % max_delta
        self.assertTrue(max_delta < 5e-4, err_msg)
        plt.figure('Gaussian test')
        pdf.plot(label='pdf', scale=numpy.sqrt(2. * numpy.pi))
        cdf.plot(label='cdf')
        plt.axis([-3.5, 3.5, 0., 1.1])
        plt.legend()

    def test_sort(self):
        """Test the automatic sorting functionality.
        """
        _x = numpy.random.sample(100)
        _y = _x**2
        s = xInterpolatedUnivariateSplineLinear(_x, _y)
        _x.sort()
        self.assertTrue((s.x == _x).all())
        self.assertTrue((s.y == _x**2).all())

    def test_non_unique(self):
        """The spline constructor must fail when non-unique values are passed.
        """
        _x = numpy.array([1, 1, 2, 3, 4])
        _y = _x**2
        with self.assertRaises(AssertionError):
            s = xInterpolatedUnivariateSplineLinear(_x, _y)
            del(s)

    def test_inverse(self):
        """Test the spline inversion.
        """
        x = numpy.linspace(0., 1., 100)
        y = x**2.
        s1 = xInterpolatedUnivariateSpline(x, y, xlabel='x', ylabel='y')
        s2 = s1.inverse()
        self.assertTrue(numpy.allclose(s2.y, numpy.sqrt(y)))
        s3 = s1.inverse(0.2, 0.8)
        plt.figure('Inverse spline')
        s2.plot()
        s3.plot()

    def test_plot(self):
        """Test the plotting method.
        """
        plt.figure('Basic plotting')
        self.s1.plot(color='red')
        self.s3.plot(overlay=True)



class TestUnivariateSplineLog(unittest.TestCase):

    """Unit test for univariate splines in log space.
    """

    def power_law_integral(self, norm, index, xmin, xmax):
        """Convenience function for the test.
        """
        return norm / (1. - index) * (xmax**(1. -index) - xmin**(1. - index))

    def test_power_law(self):
        """Basic powe-law test
        """
        norm = 1.
        index = 2.
        emin = 1.
        emax = 15.
        num_points = 5
        _x = numpy.logspace(numpy.log10(emin), numpy.log10(emax), num_points)
        _y = norm * _x**(-index)
        slin = xInterpolatedUnivariateSplineLinear(_x, _y)
        slog = xInterpolatedUnivariateLogSplineLinear(_x, _y)
        # Test normalization.
        target_norm = self.power_law_integral(norm, index, emin, emax)
        log_norm = slog.norm()
        delta = abs(target_norm - log_norm)/target_norm
        msg = 'delta = %.3e' % delta
        self.assertTrue(delta < 0.01, msg)
        # Basic plotting
        plt.figure('Simple power law')
        slin.plot(logx=True, logy=True, overlay=True, label='Linear')
        slog.plot(logx=True, logy=True, label='Logarithmic')
        # Test extrapolation.
        _x = numpy.logspace(numpy.log10(emax), numpy.log10(2 * emax), 10)
        _y = slog(_x)
        plt.plot(_x, _y, color=last_line_color(), ls='dashed')
        plt.legend()

    def test_power_law_cutoff(self):
        """Test a power-law with exponential cutoff.
        """
        norm = 1.
        index = 2.
        emin = 1.
        emax = 15.
        cutoff = 5.
        num_points = 10
        _x = numpy.logspace(numpy.log10(emin), numpy.log10(emax), num_points)
        _y = norm * _x**(-index) * numpy.exp(-_x / cutoff)
        slin = xInterpolatedUnivariateSplineLinear(_x, _y)
        slog = xInterpolatedUnivariateLogSplineLinear(_x, _y)
        # Basic plotting
        plt.figure('Power law with exponential cutoff')
        slin.plot(logx=True, logy=True, overlay=True, label='Linear')
        slog.plot(logx=True, logy=True, label='Logarithmic')
         # Test extrapolation.
        _x = numpy.logspace(numpy.log10(emax), numpy.log10(2 * emax), 10)
        _y = slog(_x)
        plt.plot(_x, _y, color=last_line_color(), ls='dashed')
        plt.legend()



class TestBivariateSpline(unittest.TestCase):

    """Unit test for bivariate splines.
    """

    @classmethod
    def setUpClass(cls):
        """Setup.

        Create a few objects to be used for testing.
        """
        cls.x = numpy.linspace(0., 1., 100)
        cls.y = numpy.linspace(0., 2., 200)
        cls.test_points = ((0.1, 0.2), (0.5, 1.), (0.9, 1.9))
        cls.fmt = dict(xlabel='x [a.u.]', ylabel='y [a.u.]', zlabel='z [a.u.]')

    def basic_test(self, name, z, plot_countours=True):
        """Small utility function to test bivariate splines.
        """
        s = xInterpolatedBivariateSpline(self.x, self.y, z, **self.fmt)
        plt.figure('%s 2d test' % name)
        s.plot()
        if plot_countours:
            s.plot_contours()
        _x = 0.5
        plt.figure('%s 2d vslice @ x = %.3f test' % (name, _x))
        s.vslice(_x).plot()
        _y = 1.
        plt.figure('%s 2d hslice @ y = %.3f test' % (name, _y))
        s.hslice(_y).plot()
        plt.figure('%s 2d ppf test' % name)
        ppf = s.build_horizontal_ppf()
        ppf.plot()
        ppf.plot_contours()
        return s

    def test_constant(self, value=10.):
        """Test a constant spline.
        """
        z = numpy.full((self.x.size, self.y.size), value)
        s = self.basic_test('Constant', z, False)
        for _x, _y in self.test_points:
            self.assertAlmostEqual(s(_x, _y), value)

    @staticmethod
    def _bilinear_func(x, y):
        """Test bilinear function.
        """
        return x + 2. * y

    def test_bilinear(self):
        """Test a bilinear array.

        Note that here we construct the z array by hand, without using the
        numpy meshgrid facility, in order to remove the ambiguity between the
        matlab and the scipy indexing conventions.
        """
        z = numpy.zeros((self.x.size, self.y.size), dtype=float)
        for i, _x in enumerate(self.x):
            for j, _y in enumerate(self.y):
                z[i, j] = self._bilinear_func(_x, _y)
        s = self.basic_test('Bilinear', z)
        for _x, _y in self.test_points:
            self.assertAlmostEqual(s(_x, _y), self._bilinear_func(_x, _y))

    def test_bilinear_func(self):
        """Test a bilinear spline.
        """
        s = self.basic_test('Bilinear function', self._bilinear_func)
        for _x, _y in self.test_points:
            self.assertAlmostEqual(s(_x, _y), self._bilinear_func(_x, _y))



if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)
